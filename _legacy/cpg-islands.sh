
#'
#' This file creates a bed file of CpG islands in the Tanytarsus gracilentus
#' reference assembly.
#'

#' I first clone `lucasnell/TaJoCGI`, then cd to its directory.
#' Then I install the faster Cython version:

python3 cySetup.py build_ext --inplace

#' Change the following on your system for inputs, outputs, and threads:
FULL_REFERENCE="${HOME}/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz"
REFERENCE=$(basename ${FULL_REFERENCE} | sed 's/.gz//g')
OUT_BED="${HOME}/_data/Tgraci_cpg.bed"
N_THREADS=6

#' I had some issues using a gzipped input file for some reason
gunzip -c ${FULL_REFERENCE} \
    > ${REFERENCE}


#' Then I run it on the assembly.
./TaJoCGI.py -c ${N_THREADS} -o ${OUT_BED} ${REFERENCE}


# Not sure why, but it's only identifying 7 CpG islands...
