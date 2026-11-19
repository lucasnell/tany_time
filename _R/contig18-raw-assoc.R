#'
#' How much did mapping exaggerate the CONTIG18 locus's association with
#' hm_logN_6? (Owner, 2026-11-19; exploratory, below threshold.) Compares
#' the haplotype's association using
#'   raw     per-pool read counts of haplotypes A and B from reads anchored
#'           next to the locus (`_R/contig18-flank.R` →
#'           `_data/contig18-flank-counts.csv`), and
#'   mapped  per-pool counts from the first rerun's BayPass inputs (mp050):
#'           B-allele and total reads averaged over the 9 SNPs (B is the
#'           reference allele at 373,124 and the alternate at the others,
#'           in poolfstat's labels; `_data/contig18-haplotype.txt`).
#' 17 pools (SN-1-95 excluded).
#'
#'   1. Pearson r and least-squares slope of B frequency on hm_logN_6
#'      across pools, for each, over pools with data in both.
#'   2. BayPass 3.1 (local), standard covariate model, on the two versions
#'      of the locus embedded among 2,000 random SNPs of mp050 subsample 1
#'      (so the model's hyperparameters are estimated as usual), with the
#'      subsample-1 Ω (SN-1-95 dropped) fixed: importance-sampling Bayes
#'      factor (BF, dB) and effect estimate (Beta_is), and the MCMC effect
#'      (-covmcmc: M_Beta, SD_Beta), 3 seeds each. The auxiliary model of
#'      the main runs isn't used: its Bayes factor depends on the share of
#'      associated SNPs over the whole data set, so a one-locus comparison
#'      wouldn't be meaningful. This compares raw with mapped, not with the
#'      main runs' Bayes factors.
#'
#' Writes `_data/contig18-raw-assoc.txt`.
#'

library(tidyverse)

ARCHIVE_18POOLS <- TRUE
source("_R/data-dir.R")

out_file <- "_data/contig18-raw-assoc.txt"
sink(out_file, split = TRUE)

set.seed(1236521059)
pools18 <- read_lines("_bash/samples/time-pools-18.txt")
pools <- read_lines("_bash/samples/time-pools.txt")
k18 <- match(pools, pools18)
ctg <- "CONTIG18"
lo <- 373124
hi <- 373206



# =============================================================================*
# ---- counts ----
# =============================================================================*


raw <- read_csv("_data/contig18-flank-counts.csv", show_col_types = FALSE) |>
    select(pool, A, B)

tmp <- tempfile()
untar(tany_dir("baypass-inputs", "baypass-inputs_mp050_3.tar.gz"), exdir = tmp)
d <- file.path(tmp, "baypass-inputs_mp050")
geno <- map(1:3, \(k) {
    det <- read_table(file.path(d, paste0("tany.snpdet.sub", k)),
                      col_names = c("contig", "pos", "ref", "alt"),
                      col_types = "cdcc")
    g <- data.table::fread(file.path(d, paste0("tany.genobaypass.sub", k)),
                           header = FALSE) |>
        as.matrix()
    list(det = det, g = g)
})
hap18 <- read_lines(file.path(d, "tany.poolsize")) |> str_split_1(" ") |>
    as.integer()
loc <- map(geno, \(x) {
    i <- which(x$det$contig == ctg & x$det$pos >= lo & x$det$pos <= hi)
    ref <- x$g[i, seq(1, ncol(x$g), 2), drop = FALSE][, k18, drop = FALSE]
    alt <- x$g[i, seq(2, ncol(x$g), 2), drop = FALSE][, k18, drop = FALSE]
    b <- ifelse(x$det$pos[i] == lo, 1, 0)
    list(B = ref * b + alt * (1 - b), N = ref + alt)
})
mB <- do.call(rbind, map(loc, "B"))
mN <- do.call(rbind, map(loc, "N"))
stopifnot(nrow(mB) == 9)
mapped <- tibble(pool = pools,
                 B = round(colMeans(mB)),
                 N = round(colMeans(mN))) |>
    mutate(A = N - B) |>
    select(pool, A, B)

cov <- scan("_bash/03-select-time/covariates/hm_logN_6.txt", quiet = TRUE)
stopifnot(length(cov) == length(pools))
both <- tibble(pool = pools, hm_logN_6 = cov) |>
    left_join(raw |> rename(raw_A = A, raw_B = B), by = "pool") |>
    left_join(mapped |> rename(map_A = A, map_B = B), by = "pool") |>
    mutate(across(c(raw_A, raw_B, map_A, map_B), \(x) replace_na(x, 0L)),
           raw_f = ifelse(raw_A + raw_B > 0, raw_B / (raw_A + raw_B), NA),
           map_f = ifelse(map_A + map_B > 0, map_B / (map_A + map_B), NA))
cat("== Per pool: B reads / all reads, raw (anchored reads) and mapped\n")
both |>
    mutate(across(c(raw_f, map_f, hm_logN_6), \(x) round(x, 2))) |>
    arrange(hm_logN_6) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- 1. correlation and slope ----
# =============================================================================*


ok <- !is.na(both$raw_f) & !is.na(both$map_f)
fit_one <- function(f) {
    m <- lm(f[ok] ~ both$hm_logN_6[ok])
    c(r = cor(f[ok], both$hm_logN_6[ok]), slope = unname(coef(m)[2]),
      slope_se = unname(summary(m)$coefficients[2, 2]))
}
cat("\n== 1. B frequency vs. hm_logN_6 over", sum(ok),
    "pools with data in both\n")
rbind(raw = fit_one(both$raw_f), mapped = fit_one(both$map_f)) |>
    round(4) |>
    print()



# =============================================================================*
# ---- 2. BayPass, standard covariate model ----
# =============================================================================*


#' Background: 2,000 random SNPs of subsample 1 (17 pools).
g1 <- geno[[1]]$g
bg <- sample.int(nrow(g1), 2000)
cols17 <- as.vector(rbind(2 * k18 - 1, 2 * k18))
bg_g <- g1[bg, cols17]
#' The two versions of the locus: rows 1 (raw) and 2 (mapped); counts
#' per pool as (B, A).
loc_rows <- rbind(as.vector(rbind(both$raw_B, both$raw_A)),
                  as.vector(rbind(both$map_B, both$map_A)))
unlink(tmp, recursive = TRUE)

untar(tany_dir("baypass-init", "baypass-init-mp050-sub1.tar.gz"),
      files = "baypass-init-mp050-sub1/tany-init-mp050-sub1_mat_omega.out",
      exdir = tmp)
om <- as.matrix(read.table(file.path(tmp, "baypass-init-mp050-sub1",
                                     "tany-init-mp050-sub1_mat_omega.out")))
om17 <- om[k18, k18]

run_dir <- file.path(tmp, "bp")
dir.create(run_dir)
write.table(rbind(loc_rows, bg_g), file.path(run_dir, "g.txt"),
            row.names = FALSE, col.names = FALSE)
write_lines(paste(hap18[k18], collapse = " "), file.path(run_dir, "ps.txt"))
write_lines(paste(cov, collapse = " "), file.path(run_dir, "cov.txt"))
write.table(om17, file.path(run_dir, "om.txt"), row.names = FALSE,
            col.names = FALSE)

bp <- function(seed, mcmc) {
    pre <- file.path(run_dir, sprintf("r%s_%i", ifelse(mcmc, "mc", "is"), seed))
    args <- c("-gfile", file.path(run_dir, "g.txt"),
              "-poolsizefile", file.path(run_dir, "ps.txt"),
              "-efile", file.path(run_dir, "cov.txt"),
              "-omegafile", file.path(run_dir, "om.txt"),
              "-outprefix", pre, "-seed", seed, "-nthreads", 4)
    if (mcmc) args <- c(args, "-covmcmc")
    st <- system2("baypass", args, stdout = FALSE, stderr = FALSE)
    stopifnot(st == 0)
    f <- paste0(pre, if (mcmc) "_summary_betai.out" else "_summary_betai_reg.out")
    read_table(f, col_types = cols()) |>
        filter(MRK %in% 1:2) |>
        mutate(version = c("raw", "mapped")[MRK], seed = seed)
}
is_res <- map(c(11, 22, 33), \(s) bp(s, FALSE)) |> list_rbind()
mc_res <- map(c(11, 22, 33), \(s) bp(s, TRUE)) |> list_rbind()

cat("\n== 2. BayPass standard covariate model (3 seeds)\n")
cat("Importance sampling:\n")
is_res |>
    select(version, seed, any_of(c("BF(dB)", "Beta_is", "eBPis"))) |>
    as.data.frame() |>
    print()
cat("\nMCMC (-covmcmc):\n")
mc_res |>
    select(version, seed, any_of(c("M_Beta", "SD_Beta", "eBPmc"))) |>
    mutate(z = M_Beta / SD_Beta) |>
    as.data.frame() |>
    print()
cat("\nBackground SNPs (2,000), for scale: importance-sampling BF (dB)",
    "quantiles 50/90/99%, first seed:\n")
f0 <- file.path(run_dir, "ris_11_summary_betai_reg.out")
bgbf <- read_table(f0, col_types = cols()) |> filter(MRK > 2)
print(round(quantile(bgbf[["BF(dB)"]], c(0.5, 0.9, 0.99)), 2))

unlink(tmp, recursive = TRUE)
sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
