#!/bin/bash

#'
#' Script to create files for baypass.
#'
#' Before running this on the cluster, combine the necessary sample info
#' files into a tar file:
#'
#' cd THIS_REPO \
#'     && cd _data \
#'     && export SAMP_DIR=tany-sample-info \
#'     && mkdir ${SAMP_DIR} \
#'     && cp full-sample-info.csv tany-abundances.csv ./${SAMP_DIR}/ \
#'     && tar -cz --no-xattrs --exclude ".*" -f ${SAMP_DIR}.tar.gz ${SAMP_DIR} \
#'     && scp ${SAMP_DIR}.tar.gz lnell@transfer.chtc.wisc.edu:/staging/lnell/
#'

. /app/.bashrc
conda activate baypass-env

export THREADS=$(count_threads)

# Edit these for your system:
export DATA_TAR_FULL_PATH="/staging/lnell/tany-sample-info.tar.gz"
export SAMPLE_INFO="full-sample-info.csv"
export ABUNDANCES="tany-abundances.csv"
export COUNTS_TAR_FULL_PATH="/staging/lnell/dna/snape/count_masked_time.tar"
# All sites passing `04-combine-syncs` (masked pools allowed; see below):
export COUNTS_IN="count_masked_time.sync.gz"
export NAMES_IN="count_masked_time.names.gz"
export PARENT_DIR_OUT="/staging/lnell/dna/baypass"

#' Minimum proportion of pools with data (coverage > 0 after masking) for a
#' SNP to be kept. BayPass treats a pool with 0 reads at a SNP as missing
#' data (manual v3.1, section 2.2.1).
#' PROVISIONAL value; see TODO.md.
export MIN_PROP_POOLS=0.75






#' ========================================================================
#' Inputs
#' ========================================================================


require_files ${DATA_TAR_FULL_PATH} ${COUNTS_TAR_FULL_PATH}


export OUT_DIR="baypass-inputs"
# For sample info files:
export SAMPLE_FILE_PREFIX="tany"


mkdir ${OUT_DIR}
cd ${OUT_DIR}


# Files containing sample and abundance info:
tar -xzf ${DATA_TAR_FULL_PATH} -C ./
check_exit_status "cp, extract sample info tar" $?
cd $(basename ${DATA_TAR_FULL_PATH%.tar.gz})
mv ${SAMPLE_INFO} ${ABUNDANCES} ../
cd ..
rm -r $(basename ${DATA_TAR_FULL_PATH%.tar.gz})


# Files containing counts and info:
tar -xf ${COUNTS_TAR_FULL_PATH} -C ./
check_exit_status "cp, extract counts tar" $?
cd $(basename ${COUNTS_TAR_FULL_PATH%.tar})
mv ${COUNTS_IN} ${NAMES_IN} ../
cd ..
rm -r $(basename ${COUNTS_TAR_FULL_PATH%.tar})

#' Masked pools have ".:.:.:.:.:." at a site; poolfstat needs counts, so
#' write them as 0 counts (missing data in BayPass).
gunzip -c ${COUNTS_IN} \
    | sed 's/\.:\.:\.:\.:\.:\./0:0:0:0:0:0/g' \
    | gzip > counts_missing0.sync.gz
check_exit_status "convert masked sites to 0 counts" $?
rm ${COUNTS_IN}
export COUNTS_IN="counts_missing0.sync.gz"



#' ========================================================================
#' Write baypass input files
#' ========================================================================


cat << EOF > ../baypass-inputs.R
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages({
    library(tidyverse)
    library(poolfstat)
})

suppressWarnings(n_subs <- as.integer(args[[1]]))
if (is.na(n_subs)) stop("Improper input to R script")

.n_threads <- ${THREADS}
min_prop_pools <- ${MIN_PROP_POOLS}
in_files <- list(samp_info = "${SAMPLE_INFO}",
                 names = "${NAMES_IN}",
                 counts = "${COUNTS_IN}",
                 abunds = "${ABUNDANCES}")
out_prefix <- "${SAMPLE_FILE_PREFIX}"

samp_df <- read_csv(in_files[["samp_info"]], col_types = "cDccidd")

pools <- read_lines(in_files[["names"]])
n_haps <- pools |>
    set_names() |>
    map_int(\(x) 2L * samp_df[["n_adults"]][samp_df[["biotech_id"]] == x]) |>
    as.list()

cat("Reading pool data...\n")
#' Remove non-biallelic SNPs and those with pooled minor-allele frequency
#' < 0.01 (default filters), read into pooldata object.
#' Pools with 0 coverage at a site are kept as missing data
#' (default min.cov.per.pool = -1).
pool_dat <- popsync2pooldata(sync.file = in_files[["counts"]],
                             poolsizes = n_haps,
                             poolnames = pools,
                             nthreads = .n_threads)

#' Remove SNPs near fixation because
#' "... nearly fixed markers (e.g., π < 0.05 or π > 0.95) are not expected to be
#' associated with any covariate, even if they are adjacent to strongly
#' associated SNPs, which can severely degrade the detection of clustered
#' signals." (BayPass v1.31 manual)
#'
cat("Filtering pool data...\n")
#' Means are over pools with data at that SNP. Also requires data in at
#' least `min_prop_pools` of pools.
# Create a logical vector for which SNPs aren't near fixation:
has_data <- pool_dat@readcoverage > 0
ref_probs <- pool_dat@refallele.readcount / pool_dat@readcoverage
ref_probs[!has_data] <- NA
avg_ref_probs <- rowMeans(ref_probs, na.rm = TRUE)
enough_pools <- rowMeans(has_data) >= min_prop_pools
cat(sprintf("SNPs read: %i; with data in >= %.0f%% of pools: %i\n",
            pool_dat@nsnp, 100 * min_prop_pools, sum(enough_pools)))
poly_snps <- enough_pools & avg_ref_probs > 0.1 & avg_ref_probs < 0.9
cat(sprintf("SNPs kept: %i\n", sum(poly_snps)))
# Now updating fields:
pool_dat@nsnp <- sum(poly_snps)
pool_dat@refallele.readcount <- pool_dat@refallele.readcount[poly_snps,]
pool_dat@readcoverage <- pool_dat@readcoverage[poly_snps,]
pool_dat@snp.info <- pool_dat@snp.info[poly_snps,]


#' Now convert to format for baypass.
#' If 'n_subs' is > 1, then I use 'subsamplesize = pool_dat@nsnp %/% n_subs'
#' below to generate 'n_subs' pseudo-independent datasets.
#' This is because baypass doesn't scale linearly with number of
#' threads, so including more threads doesn't help with larger datasets
#' (> 8 isn't generally worthwhile).
#' Using multiple datasets is more efficient.
#' Also, this:
#' > the sub-sampling approach allows comparing results across
#' > (pseudo-independent) sub-sets, each of which additionally has a lower
#' > level of background LD (if a thinning approach has been performed
#' > to generate the sub-sets).
#' (BayPass v2.31 manual)

cat("Writing pool data...\n")
sub_size <- ifelse(n_subs > 1L, pool_dat@nsnp %/% n_subs, -1L)
pooldata2genobaypass(pool_dat, prefix = out_prefix,
                     subsamplesize = sub_size)


#' Lastly, write a file for the covariate, log(density + 1):
cat("Writing covariate data...\n")
pop_df <- read_csv(in_files[["abunds"]], col_types = "ifdd") |>
    mutate(gen = as.integer(gen))
pools |>
    str_remove_all("^SN-|_S.?.?\$") |>
    str_split("-") |>
    map(as.integer) |>
    set_names(pools) |>
    imap(\(x, p) {
        yr <- ifelse(x[2] < 77, 2000L + x[2], 1900L + x[2])
        gen <- x[1]
        lgl <- pop_df[["gen"]] == gen & pop_df[["year"]] == yr
        tibble(pool = p, logN = pop_df[["logN"]][lgl])
    }) |>
    list_rbind() |>
    # Make it the same order as 'pools':
    arrange(match(pool, pools)) |>
    # Now extract just 'logN' and write to file:
    getElement("logN") |>
    round(4) |>
    paste(collapse = " ") |>
    write_lines(file = paste0(out_prefix, "-log_n.txt"))

EOF


chmod +x ../baypass-inputs.R

../baypass-inputs.R 3

mv ${SAMPLE_INFO} ${ABUNDANCES} ${COUNTS_IN} ${NAMES_IN} ../
cd ..
tar -czf ${OUT_DIR}_3.tar.gz ${OUT_DIR}


mv ${SAMPLE_INFO} ${ABUNDANCES} ${COUNTS_IN} ${NAMES_IN} ./${OUT_DIR}/
cd ${OUT_DIR}
rm ${SAMPLE_FILE_PREFIX}.genobaypass* ${SAMPLE_FILE_PREFIX}.poolsize \
    ${SAMPLE_FILE_PREFIX}.snpdet* ${SAMPLE_FILE_PREFIX}-log_n.txt
../baypass-inputs.R 1



#' ========================================================================
#' Deal with outputs
#' ========================================================================

# No longer needed:
rm ${SAMPLE_INFO} ${ABUNDANCES} ${COUNTS_IN} ${NAMES_IN}


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

cp ${OUT_DIR}_3.tar.gz ${OUT_DIR}.tar.gz ${PARENT_DIR_OUT}/

rm -r ${OUT_DIR} baypass-inputs.R


exit 0


# tany.genobaypass.sub1
# tany.genobaypass.sub2
# tany.genobaypass.sub3
# tany-log_n.txt
# tany.poolsize
# tany.snpdet.sub1
# tany.snpdet.sub2
# tany.snpdet.sub3
