#'
#' Sensitivity analysis added after the main results (owner, 2026-11-19):
#' a generalized Pareto (GPD) tail fit to the POD null of z (median over 10
#' runs of |M_Beta| / SD_Beta), to give p-values below the empirical floor
#' 1 / (n_pod + 1) for single strong loci.
#'
#' Why z: the GPD fit to the POD median BFs failed (step 4 of
#' `_R/baypass-pod-threshold.R`), mainly from the pile-up of POD SNPs at
#' the BF ceiling; z has no ceiling. z had the same power as BF on the
#' spiked POD (`_data/baypass-pod-spike.txt`).
#'
#' Fit and criteria are those of step 4 of `_R/baypass-pod-threshold.R`
#' (fixed 2026-11-07), applied to z: GPD by maximum likelihood to the
#' exceedances of z over thresholds at POD quantiles 0.99 (primary), 0.995,
#' and 0.999. The fit works for a setting x covariate if (a) the 95%
#' intervals of the shape xi at the 3 thresholds overlap, (b) all 3 fits
#' converge, and (c) at the primary threshold, the POD counts at or above
#' the 0.999 and 0.9999 quantiles and the maximum are inside the 95%
#' Poisson interval of the counts the fit predicts. Only where it works are
#' real SNPs given GPD p-values. These criteria are fixed in this script
#' before the fits are run and before any real SNP is scored against them.
#'
#' Real p-values (z as the statistic throughout): z at or below the primary
#' threshold u gets the empirical p, (1 + POD SNPs with z >= its z) /
#' (n_pod + 1); z above u gets (n_u / n_pod) x P(Y > z - u) from the fit.
#' A z beyond the fitted endpoint (xi < 0) gets p = 0 and is counted. Then
#' Storey's q-values per setting x covariate, and calls at q <= 0.05 and
#' q <= 0.10.
#'
#' Reads `tany_dir("baypass-pod", "pod-null-stats.rds")`
#' (`_R/baypass-pod-stat-compare.R`) and `_data/baypass-bf.rds`
#' (`_R/baypass-read-bf.R`). Writes `_data/baypass-pod-gpd-z.txt`,
#' `_data/baypass-pod-gpd-z.csv`, `_plots/baypass-pod-gpd-z-qq.pdf`, and
#' `_data/baypass-qvalues-gpd-z.rds`.
#'

library(tidyverse)

#' Analyses the first 2026-11 rerun (18 pools; archived):
ARCHIVE_18POOLS <- TRUE
source("_R/data-dir.R")

#' GPD functions from `_R/baypass-pod-threshold.R` (only their definitions
#' are evaluated):
local({
    for (x in parse("_R/baypass-pod-threshold.R")) {
        if (is.call(x) && identical(x[[1]], as.name("<-")) &&
            as.character(x[[2]]) %in% c("gpd_nll", "gpd_fit", "gpd_surv")) {
            eval(x, envir = globalenv())
        }
    }
})
stopifnot(exists("gpd_fit"), exists("gpd_surv"))

out_file <- "_data/baypass-pod-gpd-z.txt"
sink(out_file, split = TRUE)

pod_null <- read_rds(tany_dir("baypass-pod", "pod-null-stats.rds"))
u_probs <- c(0.99, 0.995, 0.999)



# =============================================================================*
# ---- null fits and criteria ----
# =============================================================================*


gpd_df <- crossing(pod_null |> distinct(tag, covariate), u_prob = u_probs) |>
    mutate(fit = pmap(list(tag, covariate, u_prob), \(t, cv, up) {
        x <- pod_null$z[pod_null$tag == t & pod_null$covariate == cv]
        u <- unname(quantile(x, up))
        f <- gpd_fit(x[x > u] - u)
        tibble(u = u, n_u = sum(x > u), n_pod = length(x), sigma = f$sigma,
               xi = f$xi, se_xi = f$se_xi, converged = f$converged,
               endpoint = ifelse(f$xi < 0, u - f$sigma / f$xi, Inf))
    })) |>
    unnest(fit)

gof_df <- gpd_df |>
    filter(u_prob == u_probs[1]) |>
    mutate(chk = pmap(list(tag, covariate, u, n_u, sigma, xi), \(t, cv, u, n_u, s, xi) {
        x <- pod_null$z[pod_null$tag == t & pod_null$covariate == cv]
        lv <- c(quantile(x, c(0.999, 0.9999)), max = max(x))
        obs <- map_int(lv, \(l) sum(x >= l))
        pred <- n_u * gpd_surv(lv - u, s, xi)
        tibble(level = names(lv), z = unname(lv), observed = obs,
               predicted = pred, lo = qpois(0.025, pred),
               hi = qpois(0.975, pred), inside = obs >= lo & obs <= hi)
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

write_csv(gpd_df, "_data/baypass-pod-gpd-z.csv")

cat("== 1. GPD fits to the POD null of median z\n")
gpd_df |>
    mutate(across(c(u, sigma, xi, se_xi, endpoint), \(x) round(x, 3))) |>
    as.data.frame() |>
    print()
cat("\nObserved vs. predicted POD counts (threshold at quantile 0.99):\n")
gof_df |>
    mutate(across(c(z, predicted), \(x) round(x, 2))) |>
    as.data.frame() |>
    print()
cat("\nDoes the fit work? (criteria in the script header)\n")
print(as.data.frame(works_df))

qq_df <- gpd_df |>
    filter(u_prob == u_probs[1]) |>
    mutate(qq = pmap(list(tag, covariate, u, sigma, xi), \(t, cv, u, s, xi) {
        y <- sort(pod_null$z[pod_null$tag == t & pod_null$covariate == cv])
        y <- y[y > u] - u
        pp <- (seq_along(y) - 0.5) / length(y)
        fitted <- if (abs(xi) < 1e-8) -s * log(1 - pp) else
            s / xi * ((1 - pp)^(-xi) - 1)
        tibble(empirical = y + u, fitted = fitted + u)
    })) |>
    select(tag, covariate, qq) |>
    unnest(qq)
p_qq <- qq_df |>
    ggplot(aes(fitted, empirical)) +
    geom_abline(linetype = 2, color = "gray50") +
    geom_point(size = 0.5) +
    facet_grid(covariate ~ tag) +
    labs(x = "GPD quantile (z)", y = "POD median z") +
    theme_bw()
ggsave("_plots/baypass-pod-gpd-z-qq.pdf", p_qq, width = 7, height = 6)



# =============================================================================*
# ---- real SNPs, where the fit works ----
# =============================================================================*


bf_df <- read_rds("_data/archive-18pools/baypass-bf.rds")
ok <- works_df |> filter(works)
cat("\n== 2. q-values with GPD tail p (z), for the", nrow(ok),
    "setting x covariate pairs where the fit works\n")

q_df <- ok |>
    select(tag, covariate) |>
    pmap(\(tag, covariate) {
        x <- bf_df[bf_df$tag == tag & bf_df$covariate == covariate, ]
        nz <- sort(pod_null$z[pod_null$tag == tag &
                                  pod_null$covariate == covariate])
        g <- gpd_df |>
            filter(tag == !!tag, covariate == !!covariate,
                   u_prob == u_probs[1])
        p_emp <- (1 + length(nz) - findInterval(x$z, nz, left.open = TRUE)) /
            (length(nz) + 1)
        p <- ifelse(x$z > g$u,
                    g$n_u / g$n_pod * gpd_surv(x$z - g$u, g$sigma, g$xi),
                    p_emp)
        qv <- qvalue::qvalue(p)
        x |>
            select(tag, covariate, sub, mrk, contig, pos, bf_db, z) |>
            mutate(p_emp_z = p_emp, p_gpd_z = p, q_gpd_z = qv$qvalues,
                   pi0 = qv$pi0, beyond_endpoint = x$z >= g$endpoint)
    }) |>
    list_rbind()

if (nrow(q_df) > 0) {
    write_rds(q_df, "_data/baypass-qvalues-gpd-z.rds", compress = "gz")
    q_df |>
        group_by(tag, covariate) |>
        summarize(snps = n(), pi0 = round(pi0[[1]], 3),
                  min_p = signif(min(p_gpd_z), 3),
                  min_q = signif(min(q_gpd_z), 3),
                  q05 = sum(q_gpd_z <= 0.05), q10 = sum(q_gpd_z <= 0.10),
                  beyond_endpoint = sum(beyond_endpoint), .groups = "drop") |>
        as.data.frame() |>
        print()
    cat("\nSNPs at q <= 0.10:\n")
    q_df |>
        filter(q_gpd_z <= 0.10) |>
        arrange(tag, covariate, q_gpd_z) |>
        select(tag, covariate, contig, pos, bf_db, z, p_gpd_z, q_gpd_z) |>
        mutate(across(c(z, p_gpd_z, q_gpd_z), \(x) signif(x, 3))) |>
        as.data.frame() |>
        print()
} else {
    cat("No setting x covariate pair passes the fit criteria; no real SNP",
        "is scored.\n")
}

sink()
read_lines(out_file) |>
    str_trim("right") |>
    write_lines(out_file)
