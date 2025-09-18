

library(tidyverse)
library(lme4)

library(parallel)
options(mc.cores = max(1L, parallel::detectCores()-2L))


# =============================================================================*
# Data ----
# =============================================================================*

pop_df <- read_csv("_data/tany-abundances.csv", col_types = "ifdd") |>
    mutate(date = year + ifelse(gen == "a", 0, 0.5),
           crash = N <= 30) |>
    arrange(date) |>
    # generations and years since crash:
    mutate(gsc = map_int(1:n(),
                         \(i) {
                             idx <- which(crash[1:i])
                             if (length(idx) == 0) return(NA)
                             # Because we have abunds for all generations,
                             # we can just use indices:
                             return(i - tail(idx,1))
                         }),
           ysc = map_int(1:n(),
                         \(i) {
                             idx <- which(crash[1:i])
                             if (length(idx) == 0) return(NA)
                             return(year[i] - year[tail(idx,1)])
                         }))

# To more quickly access rows below:
pop_list <- pop_df |> select(-1:-2) |> split(row(pop_df[,1]))


# # To visualize what crashes signify:
# pop_df |>
#     ggplot(aes(date, logN)) +
#     geom_vline(xintercept = c(1992, 1997, 1999, 2000, 2002, 2007, 2008, 2009,
#                               2010, 2013, 2014, 2015), color = "gray90") +
#     geom_line() +
#     geom_point(data = pop_df |> filter(N <= 30)) +
#     theme_classic()


if (!file.exists("_data/allele-freqs-time-crash.rds")) {

    # Takes ~12 sec w 6 threads
    freq_df <- read_rds("_data/allele-freqs-time.rds") |>
        mutate(n_eff_int = n_eff |> round() |> as.integer(),
               info = str_split(pool, "-|_"),
               year = map_int(info, \(x) as.integer(x[3])),
               year = ifelse(year < 70L, 2000L + year, 1900L + year),
               gen = map_chr(info, \(x) letters[as.integer(x[2])]) |>
                   factor(levels = letters[1:2]),
               info = 2 * ((year + (as.integer(gen) - 1) / 2) - 1977) + 1,
               info = mclapply(info, \(idx) pop_list[[idx]])) |>
        unnest(info) |>
        select(year, gen, everything())

    write_rds(freq_df, "_data/allele-freqs-time-crash.rds", compress = "gz")

} else {

    freq_df <- read_rds("_data/allele-freqs-time-crash.rds")

}


mod_freq_df <- freq_df |>
    mutate(locus = interaction(contig, pos, drop = TRUE, sep = ":"),
           pool = interaction(year, gen, drop = TRUE, sep = ":")) |>
    select(year, gen, pool, locus, freq, n_eff, n_eff_int, N, logN, date,
           crash, gsc,  ysc)



# =============================================================================*
# Model fitting ----
# =============================================================================*

forms <- list(freq ~ (1 | locus) + (1 | pool),
              freq ~ gsc + (1 + gsc | locus) + (1 | pool),
              freq ~ ysc + (1 + ysc | locus) + (1 | pool),
              freq ~ logN + (1 + logN | locus) + (1 | pool))

# # Not using `crash` as independent variable because it causes a rank
# # deficient model matrix.
# # This is presumably because levels of `crash` are nested within levels
# # of `pool`.
# for (i in 1:length(forms)) {
for (i in 4) {
    t0 <- Sys.time()
    mod <- glmer(forms[[i]], family = binomial, weights = n_eff,
                  data = mod_freq_df, nAGQ = 0)
    write_rds(mod, sprintf("_data/crash-mod%i.rds", i-1L), compress = "gz")
    t1 <- Sys.time()
    dt <- as.numeric(difftime(t1, t0, units = "min"))
    cat(sprintf("Model %i took %.2f min\n", i-1L, dt))
}; rm(t0, t1, mod, dt)
# # Model 0 took 6.23 min
# # Model 1 took 104.31 min
# # Model 2 took 128.70 min
# # Model 3 took 17.93 min
#

map(0:3, \(i) {
    .aic <- AIC(read_rds(sprintf("_data/crash-mod%i.rds", i)))
    tibble(mod = paste0("mod", i), AIC = .aic)
}) |>
    list_rbind() |>
    mutate(dAIC = AIC - min(AIC)) |>
    arrange(AIC)

# # A tibble: 4 × 3
#   mod        AIC  dAIC
#   <chr>    <dbl> <dbl>
# 1 mod3  8666310.    0
# 2 mod1  8671939. 5629.
# 3 mod2  8672240. 5930.
# 4 mod0  8674539. 8228.


mod3 <- read_rds("_data/crash-mod3.rds")


# # https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-significance-of-random-effects
# # Takes 8.3 min
# cV <- ranef(mod3, condVar = TRUE)
# write_rds(cV, "_data/crash-mod3-cV.rds", compress = "gz")

cV <- read_rds("_data/crash-mod3-cV.rds")

# Getting sample sizes per group:
locus_n_map <- mod_freq_df |>
    group_by(locus) |>
    # summarize(n = mean(n_eff)) |>
    summarize(n = min(n_eff)) |>
    (\(x) {
        m <- as.list(x$n)
        names(m) <- x$locus
        return(m)
    })()

# Rough estimate of p-values using estimate +/- (SE* 1.96)
pval_df <- cV |>
    as.data.frame() |>
    filter(grpvar == "locus", term == "logN") |>
    select(-grpvar, -term) |>
    rename(locus = grp, val = condval, stdev = condsd) |>
    as_tibble() |>
    mutate(n = map_dbl(locus, \(l) locus_n_map[[l]]),
           se = stdev / sqrt(n),
           z = abs(val) / se,
           pval = 2 * pnorm(-abs(z)),
           phred = -log10(pval))


pval_df |>
    arrange(pval)

pval_df$phred |> hist()

# Candidate loci:
candi_df <- pval_df |>
    filter(pval < 1e-5) |>
    arrange(pval) |>
    mutate(locus = locus |> paste() |> str_split(":"),
           contig = map_chr(locus, \(x) x[[1]]) |>
               factor(levels = levels(freq_df$contig)),
           pos = map_chr(locus, \(x) x[[2]]) |> as.integer()) |>
    select(contig, pos, everything(), -locus)




# =============================================================================*
# Matching loci to genes ----
# =============================================================================*


#' This removes transcript id that's an extension separated by
#' either the last '-' or '.'
#' This converts the transcript id to a gene id.
trans2genes <- function(trans) {
    locs <- str_locate_all(trans, "\\.|-") |>
        map_int(\(x) {
            if (nrow(x) == 0) return(NA_integer_)
            x[nrow(x),"end"]
        })
    out <- character(length(trans))
    out[is.na(locs)] <- trans[is.na(locs)]
    out[!is.na(locs)] <- str_sub(trans[!is.na(locs)], 1,
                                 locs[!is.na(locs)]-1)
    return(out)
}



# Lists to lookup HOGs for each gene and GO terms for HOG:
hog_go_list <- "~/_data/midgenomes/ortho-analyses/orthofinder-extraction/All_HOG_GO/N0-GO-by-HOG.tsv" |>
    read_tsv(col_types = "cc") |>
    mutate(go = str_split(go, ";")) |>
    (\(x) {
        z <- x$go
        names(z) <- x$hog
        return(z)
    })()
hog_gene_list <- "~/_data/midgenomes/ortho-analyses/orthofinder-extraction/All_HOG_GO/N0-GO-by-species-genes.tsv" |>
    read_tsv(col_types = "cccc") |>
    filter(species == "Tgraci") |>
    group_by(gene) |>
    summarize(hog = list(unique(hog))) |>
    (\(x) {
        z <- x$hog
        names(z) <- x$gene
        return(z)
    })()



gene_df <- "~/_data/Tgraci_genes.gff3" |>
    read_tsv(col_names = c("contig", "source", "type", "start", "end",
                           "score", "strand", "phase", "attributes"),
             comment = "#", col_types = "ccciicccc") |>
    mutate(contig = factor(contig, levels = levels(freq_df$contig)),
           gene = attributes |>
               str_remove_all("^ID|geneID|=") |>
               str_split(";") |>
               map_chr(\(x) x[[1]]),
           gene = trans2genes(gene)) |>
    select(contig, strand, start, end, gene) |>
    arrange(contig, strand, start)


# Functional annotation for all genes.
# (Have to use `read_lines` because it has different numbers of columns by row)
description_list <- "~/_data/midgenomes/go-terms/Tgraci_go-terms.tsv" |>
    read_lines(progress = FALSE) |>
    base::`[`(-1) |>
    str_split("\t") |>
    lapply(\(x) {
        descx <- x[startsWith(x, "description")]
        tibble(gene = x[[1]],
               desc = ifelse(length(descx) == 0, NA_character_,
                           tolower(paste(gsub("^description:", "", descx),
                                         collapse = ";"))))
    }) |>
    list_rbind() |>
    mutate(gene = trans2genes(gene)) |>
    filter(!is.na(desc)) |>
    # To remove duplicate descriptions:
    group_by(gene) |>
    summarize(desc = desc |>
                  str_split(";") |>
                  do.call(what = c) |>
                  unique() |>
                  list()) |>
    # convert to list for quick access:
    (\(x) {
        z <- x[["desc"]]
        names(z) <- x[["gene"]]
        return(z)
    })()



# GO terms for all Tgraci genes:
go_list <- "~/_data/midgenomes/go-terms/Tgraci_go-terms.tsv" |>
    read_lines(progress = FALSE) |>
    base::`[`(-1) |>
    str_split("\t") |>
    lapply(\(x) {
        gox <- x[startsWith(x, "go")]
        tibble(gene = x[[1]],
               go = ifelse(length(gox) == 0, NA_character_,
                           tolower(paste(gox, collapse = ";"))))
    }) |>
    list_rbind() |>
    mutate(gene = trans2genes(gene)) |>
    filter(!is.na(go)) |>
    # To remove duplicate GO terms and genes:
    group_by(gene) |>
    summarize(go = go |>
                  str_split(";") |>
                  do.call(what = c) |>
                  unique() |>
                  toupper() |>
                  list()) |>
    # Now add go terms from hogs:
    mutate(go = map2(gene, go, \(.gene, .go) {
        if (is.list(.go)) .go <- do.call(c, .go)
        .hog <- hog_gene_list[[.gene]]
        if (length(.hog) == 0) return(.go)
        hog_gos <- hog_go_list[[.hog]]
        .go <- unique(toupper(c(.go, hog_gos)))
        return(.go)
    })) |>
    # convert to list for quick access:
    (\(x) {
        z <- x[["go"]]
        names(z) <- x[["gene"]]
        return(z)
    })()




# Match contigs and positions to nearest gene it's inside or upstream of
# Arguments should be vectors, and it returns a list of dataframe to unnest
gene_matcher <- function(contig, pos) {
    stopifnot(length(contig) == length(pos))
    map2(contig, pos, \(cg, ps) {
        lgl <- gene_df$contig == cg &
            ((gene_df$strand == "+" & ps <= gene_df$end) |
                 (gene_df$strand == "-" & ps >= gene_df$start))
        d <- gene_df[lgl,]
        if (nrow(d) == 0) return(tibble(gene = NA_character_, dist = NA_integer_))
        dists <- ifelse(d$strand == "+", d$start - ps, ps - d$end)
        insides <- dists < 0
        if (sum(insides) > 1) {
            gn <- d$gene[insides]
            md <- dists[insides]
        } else {
            md <- min(dists)
            gn <- d$gene[dists == md]
        }
        return(tibble(gene = gn, dist = md))
    })
}




candi_genes_df <- candi_df |>
    mutate(gene_info = gene_matcher(contig, pos)) |>
    unnest(gene_info) |>
    filter(!is.na(gene)) |>
    filter(dist < 5000)

candi_genes_df

#' Extracting sequence for one that didn't match with gene:
library(jackalope)
contig <- candi_df$contig[[1]]
pos <- candi_df$pos[[1]]

ref <- read_fasta("~/_data/midgenomes/assemblies/Tgraci_assembly.fasta.gz")
chrom <- ref$chrom(which(ref$chrom_names() == contig))
# 10kb upstream and belowstream
seq <- str_sub(max(1, pos - 10e3), min(nchar(chrom), pos + 10e3))
write_lines(seq, sprintf("~/Desktop/Tgraci_%s_%i.fa", contig, pos))

# LEFT OFF ----
#' Now blast this sequence...



#' --------
#' Blasting top genes:
#' --------
#'
#' Using BLASTP on AA sequences for each gene, matching to swissprot database.
#'
#' 1. CONTIG41-_anno2.g1857
#'      - Plays an olfactory role that is not restricted to pheromone
#'        sensitivity.
#'        https://www.uniprot.org/uniprotkb/Q9VDD3/entry
#' 2. CONTIG44-_anno1.g2930
#'      - Required for normal morphology and function of ciliated sensory organs.
#'        https://www.uniprot.org/uniprotkb/Q9VPF0/entry
#' 3. CONTIG44+_anno1.g3490
#'      - cGMP-dependent protein kinase 1 (cell signalling)
#'        https://www.uniprot.org/uniprotkb/Q03042/entry
#' 4. CONTIG11+_anno1.g10789
#'      - Plays a role in guidance and morphology of central and peripheral
#'        axons and in synaptic morphology.
#'        https://www.uniprot.org/uniprotkb/Q9VF87/entry
#' 5. CONTIG39-_anno1.g6697
#'      - An integral component of the septate junction. May play a role in
#'        cell-cell interactions that are necessary for proper development.
#'        Vital for embryonic development.
#'        https://www.uniprot.org/uniprotkb/Q9V8R9/entry
#'





candi_funcs_df <- candi_genes_df |>
    mutate(desc = map(gene, \(g) {
        if (! g %in% names(description_list)) return(NA_character_)
        description_list[[g]]
    }),
    go = map(gene, \(g) {
        if (! g %in% names(go_list)) return(NA_character_)
        go_list[[g]]
    })) |>
    filter(!is.na(go), map_int(go, length) > 0)

candi_funcs_df |>
    select(contig, pos, pval, gene:go) |>
    mutate(across(desc:go, \(x) map_chr(x, \(xx) paste(xx, collapse = ";")))) |>
    select(-desc) |>
    print(n = 30)



# library(GO.db)
GOBPOFFSPRING <- GO.db::GOBPOFFSPRING

#' Get a GO term plus all its offspring. Can optionally be done recursively.
#'
get_offs <- function(go, recursive = TRUE) {
    one_go <- function(g) {
        stopifnot(length(g) == 1)
        if (recursive) {
            g1 <- g
            g2 <- rep(NA_character_, 2)
            while (length(g2) > length(g1)) {
                if (all(!is.na(g2))) g1 <- g2
                g2 <- map(g1, \(.g) c(.g, GOBPOFFSPRING[[.g]])) |>
                    do.call(what = c) |>
                    na.exclude() |>
                    paste() |>
                    unique()
            }
        } else {
            g2 <- map(g, \(.g) c(.g, GOBPOFFSPRING[[.g]])) |>
                do.call(what = c) |>
                na.exclude() |>
                paste() |>
                unique()
        }
        return(sort(g2))
    }
    all_gos <- map(go, one_go) |> do.call(what = c)
    return(all_gos)
}

get_term <- function(gos) {
    x <- suppressMessages(
        AnnotationDbi::select(GO.db::GO.db, gos, "TERM")[,"TERM"])
    if (!is.null(names(gos))) names(x) <- names(gos)
    return(x)
}


focal_go_list <- list(anoxia = "GO:0034059" |> get_offs(),
                      hypoxia = "GO:0001666" |> get_offs(),
                      defense = "GO:0098542" |> get_offs(),
                      toxic = "GO:0009636" |> get_offs())



focal_candi_funcs_df <- candi_funcs_df |>
    mutate(anoxia = map_lgl(go, \(g) any(g %in% focal_go_list[["anoxia"]])),
           hypoxia = map_lgl(go, \(g) any(g %in% focal_go_list[["hypoxia"]])),
           defense = map_lgl(go, \(g) any(g %in% focal_go_list[["defense"]])),
           toxic = map_lgl(go, \(g) any(g %in% focal_go_list[["toxic"]]))) |>
    filter(if_any(anoxia:toxic)) |>
    select(contig, pos, pval, gene, dist, anoxia:toxic)



focal_candi_funcs_df$gene[[3]] |>
    (\(x) go_list[[x]])() |>
    keep(\(x) x %in% do.call(c, focal_go_list)) |>
    set_names() |>
    get_term()

#' --------
#' Blasting top genes:
#' --------
#'
#' Using BLASTP on AA sequences for each gene, matching to swissprot database.
#'
#' 1. CONTIG44-_anno1.g3108
#'      -
#' 2. CONTIG25-_anno1.g4894
#'      -
#' 3. CONTIG30-_anno1.g4688
#'      -
#' 4. CONTIG33+_anno1.g14012
#'      -
#' 5. CONTIG43-_anno1.g5259
#'      -
#'

