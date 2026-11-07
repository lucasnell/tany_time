#'
#' Bayes factor (BF) cutoff from pseudo-observed data (PODs;
#' `_R/baypass-pod-inputs.R` → `03-select-time/04-baypass-pod`). POD SNPs
#' have no association with any covariate, so their BFs from the covariate
#' model are a null distribution for the real runs (same Ω, pool sizes,
#' depths, missing data, filters, seeds, and MCMC settings).
#'
#' Reports, per missing-data setting x covariate:
#'   1. whether the POD mimics the data: FMD between the Ω re-estimated
#'      from the POD (core model) and the Ω it was simulated from, against
#'      FMD between real subsamples of a setting (0.22-0.33,
#'      `_data/baypass-omega-fmd.csv`);
#'   2. null quantiles of BF (dB) from single runs and of the median over
#'      the 10 runs (the statistic used for candidates);
#'   3. the null probability that the median BF reaches each of a grid of
#'      cutoffs, and the expected number of null SNPs passing it among the
#'      SNPs the real analysis tests (expected false positives).
#' Steps 1-4 use no real BF results. Step 5 (`compare_real <- TRUE`)
#' applies the rule fixed by the owner (2026-11-07): empirical p-value of
#' each real SNP's median BF against the POD medians of the same setting
#' and covariate (no pooling), Storey's q-values (`qvalue`), and candidates
#' at q <= 0.05 and q <= 0.10, both reported.
#'
#' Step 4: generalized Pareto fit to the POD tail (sensitivity analysis for
#' p-values below 1 / (n_pod + 1)); step 5 reports its q-values next to
#' the empirical ones.
#'
#' Reads `tany_dir("baypass-pod")` (and `_data/baypass-bf.rds` for step 5).
#' Writes `_data/baypass-pod-threshold.txt`, `_data/baypass-pod-null.csv`
#' (step 3 table), `_data/baypass-pod-gpd.csv` (step 4 fits),
#' `_plots/baypass-pod-gpd-qq.pdf`, and `_data/baypass-qvalues.rds` (step
#' 5: per real SNP, empirical and GPD p and q).
#'

library(tidyverse)
library(geigen) # geigen function inside fmd.dist

source("_R/data-dir.R")

if (!exists("compare_real")) compare_real <- FALSE

out_file <- "_data/baypass-pod-threshold.txt"
sink(out_file, split = TRUE)

pod_dir <- tany_dir("baypass-pod")
mp_tags <- c("mp050", "mp075", "mp100")
covariates <- c("hm_logN_6", "logN", "logN_lag1")
n_runs <- 10
bf_ceiling <- 10 * log10((2 * 5000 - 1) * 0.99 / 0.01)
#' SNPs tested in the real analysis (all 3 subsamples; `00-baypass-inputs`
#' output, `_bash/RERUN-CHECKLIST.md` I1):
n_real <- c(mp050 = 1023857, mp075 = 955748, mp100 = 528372)
cutoffs <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, bf_ceiling)

fmd.dist <- function(mat1, mat2) {
    return(sqrt(sum(log(geigen(mat1,mat2)$values)**2)))
}

#' Extract one member of a job's tar to a temporary folder; returns its path.
untar_one <- function(name, member) {
    tmp <- tempfile()
    untar(file.path(pod_dir, paste0(name, ".tar")),
          files = file.path(name, member), exdir = tmp)
    file.path(tmp, name, member)
}
read_omega <- function(fn) {
    read_table(fn, col_names = FALSE, col_types = cols()) |>
        as.matrix() |>
        unname()
}



# =============================================================================*
# ---- 1. does the POD mimic the data? ----
# =============================================================================*


cat("== 1. FMD: Omega re-estimated from POD vs. Omega simulated from\n")
tibble(tag = mp_tags) |>
    mutate(fmd = map_dbl(tag, \(t) {
        nm <- sprintf("baypass-pod-%s-core-run1", t)
        om_pod <- read_omega(untar_one(nm, sprintf("tany-pod-%s-core-run1_mat_omega.out", t)))
        om_sim <- read_omega(untar_one(nm, sprintf("tany-init-%s-sub1_mat_omega.out", t)))
        fmd.dist(om_pod, om_sim)
    })) |>
    as.data.frame() |>
    print()
cat("(Real subsamples within a setting: FMD 0.22-0.33.)\n")



# =============================================================================*
# ---- 2. null BF distribution ----
# =============================================================================*


job_df <- crossing(tag = mp_tags, covariate = covariates, run = 1:n_runs) |>
    mutate(name = sprintf("baypass-pod-%s-%s-run%i", tag, covariate, run),
           found = file.exists(file.path(pod_dir, paste0(name, ".tar"))))
cat(sprintf("\nPOD covariate-model outputs found: %i of %i\n",
            sum(job_df$found), nrow(job_df)))
stopifnot(all(job_df$found))

null_df <- job_df |>
    group_by(tag, covariate) |>
    group_split() |>
    map(\(jd) {
        map(jd$name, \(nm) {
            untar_one(nm, sprintf("%s_summary_betai.out.gz",
                                  str_replace(nm, "^baypass-pod", "tany-pod"))) |>
                read_table(col_types = "iidddd", skip = 1,
                           col_names = c("cov", "mrk", "m_beta", "sd_beta",
                                         "pip", "bf_db"),
                           progress = FALSE) |>
                select(mrk, bf_db)
        }) |>
            list_rbind() |>
            mutate(tag = jd$tag[[1]], covariate = jd$covariate[[1]])
    }) |>
    list_rbind()

med_df <- null_df |>
    group_by(tag, covariate, mrk) |>
    summarize(n_runs = n(), bf_db = median(bf_db), .groups = "drop")
stopifnot(all(med_df$n_runs == n_runs))

probs <- c(0.99, 0.999, 0.9999)
q_tbl <- function(x, statistic) {
    x |>
        group_by(tag, covariate) |>
        summarize(statistic = statistic,
                  snps = n(),
                  q99 = quantile(bf_db, probs[1]),
                  q999 = quantile(bf_db, probs[2]),
                  q9999 = quantile(bf_db, probs[3]),
                  max = max(bf_db),
                  at_ceiling = sum(bf_db >= bf_ceiling - 1e-3),
                  .groups = "drop")
}
cat("\n== 2. Null BF (dB) quantiles\n")
bind_rows(q_tbl(null_df, "single run"), q_tbl(med_df, "median of 10")) |>
    arrange(tag, covariate, statistic) |>
    mutate(across(q99:max, \(x) round(x, 2))) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- 3. null exceedance and expected false positives ----
# =============================================================================*


#' p_null = proportion of POD SNPs with median BF >= cutoff (0 if none: an
#' upper bound is then ~3 / n_pod, the 95% bound for zero events);
#' expected_fp = p_null x SNPs tested in the real setting.
exceed_df <- crossing(med_df |> distinct(tag, covariate), cutoff = cutoffs) |>
    mutate(n_pod = map2_int(tag, covariate, \(t, cv) {
        sum(med_df$tag == t & med_df$covariate == cv)
    }),
    n_exceed = pmap_int(list(tag, covariate, cutoff), \(t, cv, ct) {
        sum(med_df$bf_db[med_df$tag == t & med_df$covariate == cv] >= ct - 1e-3)
    }),
    p_null = n_exceed / n_pod,
    expected_fp = p_null * n_real[tag])

write_csv(exceed_df, "_data/baypass-pod-null.csv")

cat("\n== 3. Median BF >= cutoff: null proportion and expected false positives\n")
cat("(SNPs tested: ", paste(names(n_real), n_real, sep = " = ", collapse = ", "),
    ")\n", sep = "")
exceed_df |>
    mutate(cutoff = round(cutoff, 2), p_null = signif(p_null, 3),
           expected_fp = round(expected_fp, 1)) |>
    select(tag, covariate, cutoff, n_exceed, p_null, expected_fp) |>
    as.data.frame() |>
    print()



# =============================================================================*
# ---- 4. tail fit (sensitivity) ----
# =============================================================================*


#' The empirical p can't go below 1 / (n_pod + 1). As a sensitivity
#' analysis (owner, 2026-11-07; not used for candidates), the upper tail of
#' each setting x covariate's POD median BFs is fitted with a generalized
#' Pareto distribution (GPD; peaks over threshold, as for permutation
#' p-values in Knijnenburg et al. 2009, doi:10.1093/bioinformatics/btp211):
#' for exceedances y = x - u over threshold u,
#'   P(Y > y) = (1 + xi y / sigma)^(-1 / xi)   (exp(-y / sigma) if xi = 0),
#' with shape xi and scale sigma by maximum likelihood. For xi < 0 the tail
#' ends at u - sigma / xi (BFs have a ceiling, so xi < 0 is expected).
#' Tail p of a value x > u: (share of POD above u) x P(Y > x - u).
#'
#' Thresholds: POD quantiles 0.99 (primary: 4,000 exceedances at 400,000
#' SNPs), 0.995, and 0.999. The fit "works" (criteria fixed 2026-11-07,
#' before any POD results) if, for each setting x covariate:
#'   (a) the 95% confidence intervals of xi at the 3 thresholds overlap;
#'   (b) at the primary threshold, the POD exceedance counts at the 0.999,
#'       0.9999 quantiles and the maximum's level are inside the 95%
#'       Poisson interval of the counts the fit predicts.
#' If not, the plan is a second 400,000-SNP POD (`TODO.md`).
gpd_nll <- function(par, y) {
    s <- exp(par[1])
    xi <- par[2]
    if (abs(xi) < 1e-8) return(length(y) * log(s) + sum(y) / s)
    z <- 1 + xi * y / s
    if (any(z <= 0)) return(1e10)
    length(y) * log(s) + (1 / xi + 1) * sum(log(z))
}
gpd_fit <- function(y) {
    o <- optim(c(log(mean(y)), 0.1), gpd_nll, y = y, hessian = TRUE)
    se <- tryCatch(sqrt(diag(solve(o$hessian))), error = \(e) c(NA, NA))
    list(sigma = exp(o$par[1]), xi = o$par[2], se_xi = se[2],
         converged = o$convergence == 0)
}
gpd_surv <- function(y, sigma, xi) {
    if (abs(xi) < 1e-8) return(exp(-y / sigma))
    z <- pmax(1 + xi * y / sigma, 0)
    z^(-1 / xi)
}

u_probs <- c(0.99, 0.995, 0.999)
gpd_df <- crossing(med_df |> distinct(tag, covariate), u_prob = u_probs) |>
    mutate(fit = pmap(list(tag, covariate, u_prob), \(t, cv, up) {
        x <- med_df$bf_db[med_df$tag == t & med_df$covariate == cv]
        u <- unname(quantile(x, up))
        f <- gpd_fit(x[x > u] - u)
        tibble(u = u, n_u = sum(x > u), n_pod = length(x), sigma = f$sigma,
               xi = f$xi, se_xi = f$se_xi, converged = f$converged,
               endpoint = ifelse(f$xi < 0, u - f$sigma / f$xi, Inf))
    })) |>
    unnest(fit)

#' (b): observed vs. predicted POD exceedances at the primary threshold.
gof_df <- gpd_df |>
    filter(u_prob == u_probs[1]) |>
    mutate(chk = pmap(list(tag, covariate, u, n_u, sigma, xi), \(t, cv, u, n_u, s, xi) {
        x <- med_df$bf_db[med_df$tag == t & med_df$covariate == cv]
        lv <- c(quantile(x, c(0.999, 0.9999)), max = max(x))
        obs <- map_int(lv, \(l) sum(x >= l))
        pred <- n_u * gpd_surv(lv - u, s, xi)
        tibble(level = names(lv), bf_db = unname(lv), observed = obs,
               predicted = pred,
               lo = qpois(0.025, pred), hi = qpois(0.975, pred),
               inside = obs >= lo & obs <= hi)
    })) |>
    select(tag, covariate, chk) |>
    unnest(chk)

works_df <- gpd_df |>
    group_by(tag, covariate) |>
    summarize(xi_overlap = max(xi - 1.96 * se_xi) <= min(xi + 1.96 * se_xi),
              converged = all(converged), .groups = "drop") |>
    left_join(gof_df |>
                  group_by(tag, covariate) |>
                  summarize(gof_inside = all(inside), .groups = "drop"),
              by = c("tag", "covariate")) |>
    mutate(works = xi_overlap & converged & gof_inside)

write_csv(gpd_df, "_data/baypass-pod-gpd.csv")

cat("\n== 4. GPD tail fit of POD median BFs (sensitivity)\n")
gpd_df |>
    mutate(across(c(u, sigma, xi, se_xi, endpoint), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()
cat("\nObserved vs. predicted POD exceedances (threshold at quantile 0.99):\n")
gof_df |>
    mutate(across(c(bf_db, predicted), \(x) round(x, 2))) |>
    as.data.frame() |>
    print()
cat("\nDoes the fit work? (criteria in the script)\n")
print(as.data.frame(works_df))

#' QQ plots: exceedances over the primary threshold vs. fitted GPD.
gpd_qq <- gpd_df |>
    filter(u_prob == u_probs[1]) |>
    mutate(qq = pmap(list(tag, covariate, u, sigma, xi), \(t, cv, u, s, xi) {
        y <- sort(med_df$bf_db[med_df$tag == t & med_df$covariate == cv])
        y <- y[y > u] - u
        pp <- (seq_along(y) - 0.5) / length(y)
        fitted <- if (abs(xi) < 1e-8) -s * log(1 - pp) else
            s / xi * ((1 - pp)^(-xi) - 1)
        tibble(empirical = y + u, fitted = fitted + u)
    })) |>
    select(tag, covariate, qq) |>
    unnest(qq)
dir.create("_plots", showWarnings = FALSE)
p_qq <- gpd_qq |>
    ggplot(aes(fitted, empirical)) +
    geom_abline(linetype = 2, color = "gray50") +
    geom_point(size = 0.5) +
    facet_grid(covariate ~ tag) +
    labs(x = "GPD quantile (dB)", y = "POD median BF (dB)") +
    theme_bw()
ggsave("_plots/baypass-pod-gpd-qq.pdf", p_qq, width = 7, height = 6)



# =============================================================================*
# ---- 5. q-values for the real SNPs ----
# =============================================================================*


#' Empirical p = (1 + POD SNPs with median BF >= the real SNP's) /
#' (POD SNPs + 1), so the smallest possible p is 1 / (n_pod + 1). Storey's
#' q-values with `qvalue` defaults (null proportion pi0 estimated by
#' smoother over lambda). Also the same with the GPD tail p (section 4,
#' primary threshold) for values above the threshold: `p_gpd`, `q_gpd`
#' (sensitivity only). A value beyond the fitted endpoint gets p_gpd = 0;
#' how many do is reported.
if (compare_real) {
    bf_df <- read_rds("_data/baypass-bf.rds")
    stopifnot(all(bf_df$n_runs == n_runs))
    q_df <- bf_df |>
        group_by(tag, covariate) |>
        group_split() |>
        map(\(x) {
            null <- sort(med_df$bf_db[med_df$tag == x$tag[[1]] &
                                          med_df$covariate == x$covariate[[1]]])
            n_ge <- length(null) - findInterval(x$bf_db - 1e-9, null)
            p <- (1 + n_ge) / (length(null) + 1)
            qv <- qvalue::qvalue(p)
            g <- gpd_df |>
                filter(tag == x$tag[[1]], covariate == x$covariate[[1]],
                       u_prob == u_probs[1])
            p_gpd <- ifelse(x$bf_db > g$u,
                            g$n_u / g$n_pod * gpd_surv(x$bf_db - g$u, g$sigma, g$xi),
                            p)
            qg <- qvalue::qvalue(p_gpd)
            x |>
                select(tag, covariate, sub, mrk, contig, pos, bf_db) |>
                mutate(p_emp = p, q = qv$qvalues, pi0 = qv$pi0,
                       p_gpd = p_gpd, q_gpd = qg$qvalues,
                       beyond_endpoint = bf_db >= g$endpoint)
        }) |>
        list_rbind()
    write_rds(q_df, "_data/baypass-qvalues.rds", compress = "gz")

    cat("\n== 5. Candidates by Storey's q-value (q: empirical p; q_gpd: GPD tail p)\n")
    q_df |>
        group_by(tag, covariate) |>
        summarize(snps = n(),
                  pi0 = pi0[[1]],
                  min_p = min(p_emp),
                  q05 = sum(q <= 0.05),
                  q10 = sum(q <= 0.10),
                  q05_gpd = sum(q_gpd <= 0.05),
                  q10_gpd = sum(q_gpd <= 0.10),
                  beyond_endpoint = sum(beyond_endpoint),
                  min_bf_db_q05 = suppressWarnings(min(bf_db[q <= 0.05])),
                  min_bf_db_q10 = suppressWarnings(min(bf_db[q <= 0.10])),
                  .groups = "drop") |>
        mutate(across(c(pi0, min_bf_db_q05, min_bf_db_q10), \(x) round(x, 3)),
               min_p = signif(min_p, 3)) |>
        as.data.frame() |>
        print()
}

sink()
#' Printed tables pad lines with spaces; strip them.
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
