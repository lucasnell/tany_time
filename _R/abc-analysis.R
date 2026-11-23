#'
#' ABC by random forests (`abcrf`; Pudlo et al. 2016; Raynal et al. 2019) for
#' crashes vs. crashes plus immigration, as fixed in `_R/abc-design.md`.
#'
#' Reads (copy from `/staging/l/lnell/dna/abc/` first):
#'   tany_dir("abc", "sims")                     abc-<variant>-<model>-b<batch>.tsv.gz
#'   tany_dir("abc", "abc-observed-stats.tsv")   observed statistics
#'
#' Steps, in the design's order:
#'   1. Reference table (variant `main`) and held-out sets (`holdout`);
#'      simulations with any missing statistic are dropped and counted.
#'   2. Prior-predictive check: is each observed statistic inside the
#'      range simulated under each model? If none of the models covers
#'      it, the design says to revise the model before any inference
#'      (`stop_if_outside`).
#'   3. Model choice among M1-M5 and between groups (immigration: M2, M3,
#'      M5; none: M1, M4), with out-of-bag prior error rates. `abcrf`
#'      gives a posterior probability only for the chosen (most voted)
#'      model or group, so each model's posterior comes from a forest for
#'      that model vs. all others: P(model) = the posterior if the model is
#'      chosen, else 1 - the posterior of "all others".
#'      Decision rule: immigration "supported" only if the immigration
#'      group is chosen with posterior >= 0.9.
#'   4. Parameters (log10) under each model with posterior >= 0.1:
#'      `regAbcrf`, median and 90% interval.
#'   5. Parameter recovery on the held-out simulations: coverage of the 90%
#'      intervals and width relative to the prior. m or I is "identifiable"
#'      only if coverage is 0.85-0.95 and the mean width is below the
#'      prior's.
#'   6. Posterior-predictive draws: 1,000 parameter sets from the most
#'      probable model, resampled from its reference table with the random
#'      forest weights for the immigration parameter (or scale, under
#'      models without immigration), to simulate with `abc-sim.py --fix`.
#'   7. Sensitivity runs available so far: leaving out SN-1-13 and SN-2-97
#'      (their statistics dropped; SN-1-95 is left out of everything); diversity statistics alone
#'      vs. temporal and among-pool statistics alone; and each simulated
#'      variant (`power`, `widen`, `older`, `mu`) whose files exist.
#'
#' Writes `_data/abc-results.txt` (report), `_data/abc-model-choice.csv`,
#' `_data/abc-params.csv`, `_data/abc-recovery.csv`,
#' `_data/abc-ppc-params.csv`, and `_plots/abc-prior-check.pdf`.
#'

library(tidyverse)
library(abcrf)

source("_R/data-dir.R")

#' For a test on stand-in data or a pilot (`ABC_TEST_DIR`, with the same
#' layout as `tany_dir("abc")`), outputs go to that folder instead.
abc_dir <- if (Sys.getenv("ABC_TEST_DIR") == "") tany_dir("abc") else
    Sys.getenv("ABC_TEST_DIR")
out_dir <- if (Sys.getenv("ABC_TEST_DIR") == "") "_data" else abc_dir
plot_dir <- if (Sys.getenv("ABC_TEST_DIR") == "") "_plots" else abc_dir
N_TREE <- 1000L
N_CORES <- max(parallel::detectCores() - 1L, 1L)
#' Design: called "supported" at a grouped posterior >= 0.9; parameters
#' estimated under models with posterior >= 0.1.
P_SUPPORT <- 0.9
P_PARAMS <- 0.1
stop_if_outside <- TRUE

set.seed(1393870127)

out_txt <- file.path(out_dir, "abc-results.txt")
sink(out_txt, split = TRUE)

say <- function(...) cat(..., "\n", sep = "")
show_df <- function(x) print(as.data.frame(x), row.names = FALSE)



# =============================================================================*
# ---- 1. reference table ----
# =============================================================================*


par_names <- c("scale", "floor", "Ne_c", "Ne_anc", "Ne_s", "T_s", "m", "I",
               "b", "Ne_old", "T_old")
models <- paste0("M", 1:5)
imm_models <- c("M2", "M3", "M5")

read_sims <- function(variant) {
    fns <- list.files(file.path(abc_dir, "sims"),
                      sprintf("^abc-%s-M[1-5]-b\\d+\\.tsv\\.gz$", variant),
                      full.names = TRUE)
    if (length(fns) == 0) return(NULL)
    fns |>
        map(\(f) read_tsv(f, show_col_types = FALSE, progress = FALSE,
                          col_types = cols(model = "c", .default = "d"))) |>
        list_rbind()
}

obs_df <- read_tsv(file.path(abc_dir, "abc-observed-stats.tsv"),
                   show_col_types = FALSE)
obs <- setNames(obs_df$value, obs_df$stat)
stat_names <- names(obs)
stopifnot(length(stat_names) == 258, !anyNA(obs))

ref <- read_sims("main")
hold <- read_sims("holdout")
stopifnot(!is.null(ref), all(stat_names %in% names(ref)))

#' Model formulas need syntactic names (pool names contain "-"):
syn <- function(x) {
    i <- names(x) %in% stat_names
    names(x)[i] <- make.names(names(x)[i])
    x
}
ref <- syn(ref)
if (!is.null(hold)) hold <- syn(hold)
names(obs) <- make.names(names(obs))
stat_names <- names(obs)

#' SN-1-95 is dropped from all analyses (owner, 2026-11-19; likely
#' contamination, `_data/sn-1-95-check.txt`). It's still a sampled pool in
#' the simulations, which doesn't change any other pool's statistics, so
#' only the statistics involving it are removed (235 of 258 remain).
excluded_pools <- make.names("SN-1-95_S59")
stat_names <- stat_names[!map_lgl(stat_names, \(s) any(str_detect(
    s, fixed(excluded_pools))))]
obs <- obs[stat_names]

clean <- function(x, label) {
    bad <- !complete.cases(x[, stat_names]) |
        !apply(is.finite(as.matrix(x[, stat_names])), 1, all)
    say(label, ": ", nrow(x), " simulations; dropped (missing or infinite ",
        "statistics): ", sum(bad))
    print(table(x$model[!bad]))
    x[!bad, ]
}
say("== 1. Simulations")
ref <- clean(ref, "Reference table")
if (!is.null(hold)) hold <- clean(hold, "Held out")
say("Time per simulation (min), reference table:")
ref |>
    group_by(model) |>
    summarize(mean = mean(t_all) / 60, max = max(t_all) / 60) |>
    mutate(across(c(mean, max), \(x) round(x, 1))) |>
    show_df()



# =============================================================================*
# ---- 2. prior-predictive check ----
# =============================================================================*


say("\n== 2. Prior-predictive check")
cover <- map(models, \(md) {
    x <- ref[ref$model == md, stat_names]
    tibble(stat = stat_names, model = md,
           inside = map2_lgl(x, obs, \(s, o) o >= min(s) & o <= max(s)),
           q = map2_dbl(x, obs, \(s, o) mean(s <= o)))
}) |>
    list_rbind()
cov_any <- cover |>
    group_by(stat) |>
    summarize(any_inside = any(inside), .groups = "drop")
say("Observed statistics inside the simulated range, by model:")
cover |>
    group_by(model) |>
    summarize(inside = sum(inside), of = n(),
              in_central_98 = sum(q >= 0.01 & q <= 0.99)) |>
    show_df()
outside <- cov_any$stat[!cov_any$any_inside]
say("Outside every model's range: ", length(outside), " of ",
    length(stat_names))
if (length(outside) > 0) {
    #' Direction and size: observed vs. the smallest and largest value over
    #' all simulations, and the time between the two pools (generations;
    #' 1977a = 0) for temporal F and F_ST.
    gen_of <- function(p) {
        yy <- as.integer(str_match(p, "^SN\\.[12]\\.(\\d{2})")[, 2])
        y <- ifelse(yy < 70, 2000L, 1900L) + yy
        2L * (y - 1977L) + as.integer(str_match(p, "^SN\\.([12])\\.")[, 2]) - 1L
    }
    out_df <- tibble(stat = outside,
                     observed = obs[outside],
                     sim_min = map_dbl(outside, \(s) min(ref[[s]])),
                     sim_max = map_dbl(outside, \(s) max(ref[[s]]))) |>
        mutate(side = ifelse(observed > sim_max, "above", "below"),
               pools = str_match(stat, "^(fc|fst)_(SN\\.\\d\\.\\d{2}_S\\d+)_(SN\\.\\d\\.\\d{2}_S\\d+)$"),
               gap = abs(gen_of(pools[, 3]) - gen_of(pools[, 4]))) |>
        select(-pools)
    say("Statistics outside every model's range:")
    out_df |>
        mutate(across(c(observed, sim_min, sim_max), \(x) signif(x, 3))) |>
        arrange(desc(observed / sim_max)) |>
        show_df()
    say("By statistic type and side:")
    out_df |>
        mutate(type = str_extract(stat, "^[a-z]+")) |>
        count(type, side) |>
        show_df()
    say("Statistics of the same types inside some model's range, by gap ",
        "(generations):")
    tibble(stat = stat_names) |>
        mutate(type = str_extract(stat, "^[a-z]+"),
               pools = str_match(stat, "^(fc|fst)_(SN\\.\\d\\.\\d{2}_S\\d+)_(SN\\.\\d\\.\\d{2}_S\\d+)$"),
               gap = abs(gen_of(pools[, 3]) - gen_of(pools[, 4])),
               outside = stat %in% outside) |>
        filter(type %in% c("fc", "fst")) |>
        mutate(gap_bin = cut(gap, c(0, 2, 5, 10, 20, 50))) |>
        count(type, gap_bin, outside) |>
        pivot_wider(names_from = outside, values_from = n, values_fill = 0,
                    names_prefix = "outside_") |>
        show_df()
}

#' First two principal components of the statistics (scaled), with the
#' observed data projected:
pca <- prcomp(ref[, stat_names], scale. = TRUE)
pc_obs <- predict(pca, as.data.frame(t(obs)))
p <- as_tibble(pca$x[, 1:2]) |>
    mutate(model = ref$model) |>
    slice_sample(n = min(nrow(ref), 20000)) |>
    ggplot(aes(PC1, PC2, color = model)) +
    geom_point(size = 0.3, alpha = 0.3) +
    geom_point(data = tibble(PC1 = pc_obs[1], PC2 = pc_obs[2]),
               color = "black", shape = 4, size = 4, stroke = 1.5) +
    scale_color_brewer(palette = "Dark2") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_classic()
ggsave(file.path(plot_dir, "abc-prior-check.pdf"), p, width = 6, height = 5)

#' `ABC_CHECK_ONLY=1`: stop here, after the prior-predictive check (for
#' pilots; no inference).
check_only <- Sys.getenv("ABC_CHECK_ONLY") == "1"
if (check_only && length(outside) == 0) {
    say("\nPrior-predictive check passed (ABC_CHECK_ONLY: stopping before ",
        "inference).")
    sink()
    read_lines(out_txt) |>
        str_subset("^Growing trees", negate = TRUE) |>
        str_trim("right") |>
        write_lines(out_txt)
    quit(save = "no")
}
if (length(outside) > 0 && stop_if_outside) {
    sink()
    read_lines(out_txt) |>
        str_subset("^Growing trees", negate = TRUE) |>
        str_trim("right") |>
        write_lines(out_txt)
    stop("Observed statistics outside every model's simulated range; ",
         "the design says to revise the model before inference.")
}



# =============================================================================*
# ---- 3. model choice ----
# =============================================================================*


#' One forest; returns the chosen model or group, its posterior, votes,
#' and the out-of-bag prior error rate.
fit_choice <- function(x, stats, group = list(), label = "") {
    d <- data.frame(modindex = factor(x$model), x[, stats])
    f <- abcrf(modindex ~ ., data = d, group = group, lda = TRUE,
               ntree = N_TREE, paral = TRUE, ncores = N_CORES)
    pr <- predict(f, obs = as.data.frame(t(obs[stats])), training = d,
                  ntree = N_TREE, paral = TRUE, ncores = N_CORES)
    list(forest = f, allocation = as.character(pr$allocation),
         post = as.numeric(pr$post.prob), votes = pr$vote[1, ], prior_err = f$prior.err,
         confusion = f$model.rf$confusion.matrix, label = label)
}

model_choice <- function(x, stats, label) {
    say("\n-- ", label, " (", length(stats), " statistics, ", nrow(x),
        " simulations)")
    all5 <- fit_choice(x, stats, label = "M1-M5")
    say("Five models: chosen ", all5$allocation, ", posterior ",
        signif(all5$post, 3), ", out-of-bag prior error ",
        signif(all5$prior_err, 3))
    say("Votes: ", paste(names(all5$votes), all5$votes, sep = " = ",
                         collapse = ", "))
    say("Out-of-bag confusion (rows: true model):")
    print(all5$confusion)
    grp <- fit_choice(x, stats,
                      group = list(c("M1", "M4"), imm_models),
                      label = "groups")
    #' With groups, `abcrf` labels them g1, g2, ... in the order given:
    imm_chosen <- grp$allocation == "g2"
    p_imm <- if (imm_chosen) grp$post else 1 - grp$post
    say("Immigration (M2+M3+M5) vs. none (M1+M4): chosen ",
        if (imm_chosen) "immigration" else "none", ", posterior ",
        signif(grp$post, 3), "; P(immigration) = ", signif(p_imm, 3),
        ", out-of-bag prior error ", signif(grp$prior_err, 3))
    say("Immigration 'supported' (immigration chosen, posterior >= ",
        P_SUPPORT, "): ", imm_chosen && grp$post >= P_SUPPORT)
    #' Each model vs. all others:
    p_model <- map_dbl(models, \(md) {
        b <- fit_choice(x, stats, group = list(md, setdiff(models, md)))
        chosen <- b$allocation == "g1"
        if (chosen) b$post else 1 - b$post
    })
    names(p_model) <- models
    say("Posterior of each model (vs. all others): ",
        paste(models, signif(p_model, 3), sep = " = ", collapse = ", "))
    tibble(analysis = label, n_stats = length(stats), n_sims = nrow(x),
           chosen_5 = all5$allocation, post_5 = all5$post,
           prior_err_5 = all5$prior_err, p_immigration = p_imm,
           prior_err_group = grp$prior_err,
           supported = imm_chosen && grp$post >= P_SUPPORT,
           !!!setNames(as.list(p_model), paste0("p_", models)))
}

say("\n== 3. Model choice")
mc <- model_choice(ref, stat_names, "main")
p_model <- unlist(mc[1, paste0("p_", models)])
names(p_model) <- models



# =============================================================================*
# ---- 4. parameters ----
# =============================================================================*


say("\n== 4. Parameters (log10; median and 90% interval)")
model_pars <- function(x) {
    p <- intersect(par_names, names(x))
    p[colSums(!is.na(x[, p])) > 0]
}
prior_bounds <- read_csv("_bash/04-diversity/abc-priors.csv",
                         show_col_types = FALSE)

fit_param <- function(x, par, stats) {
    d <- data.frame(y = log10(x[[par]]), x[, stats])
    regAbcrf(y ~ ., data = d, ntree = N_TREE, paral = TRUE, ncores = N_CORES)
}

est_models <- models[p_model >= P_PARAMS]
par_fits <- list()
par_est <- map(est_models, \(md) {
    x <- ref[ref$model == md, ]
    map(model_pars(x), \(par) {
        f <- fit_param(x, par, stat_names)
        par_fits[[paste(md, par)]] <<- f
        d <- data.frame(y = log10(x[[par]]), x[, stat_names])
        pr <- predict(f, obs = as.data.frame(t(obs)), training = d,
                      quantiles = c(0.05, 0.95), paral = TRUE,
                      ncores = N_CORES)
        tibble(model = md, param = par, median = pr$med[1],
               lo90 = pr$quantiles[1, 1], hi90 = pr$quantiles[1, 2])
    }) |>
        list_rbind()
}) |>
    list_rbind() |>
    left_join(prior_bounds |> transmute(param, prior_lo = log10(lo),
                                        prior_hi = log10(hi)),
              by = "param")
par_est |> mutate(across(where(is.double), \(x) round(x, 2))) |> show_df()



# =============================================================================*
# ---- 5. parameter recovery ----
# =============================================================================*


say("\n== 5. Parameter recovery on held-out simulations")
recov <- tibble()
if (!is.null(hold)) {
    key_pars <- c("m", "I", "scale", "floor")
    recov <- map(models, \(md) {
        x <- ref[ref$model == md, ]
        h <- hold[hold$model == md, ]
        map(intersect(key_pars, model_pars(x)), \(par) {
            f <- par_fits[[paste(md, par)]]
            if (is.null(f)) f <- fit_param(x, par, stat_names)
            d <- data.frame(y = log10(x[[par]]), x[, stat_names])
            pr <- predict(f, obs = h[, stat_names], training = d,
                          quantiles = c(0.05, 0.95), paral = TRUE,
                          ncores = N_CORES)
            truth <- log10(h[[par]])
            b <- prior_bounds[prior_bounds$param == par, ]
            tibble(model = md, param = par, n = nrow(h),
                   coverage90 = mean(truth >= pr$quantiles[, 1] &
                                         truth <= pr$quantiles[, 2]),
                   mean_width = mean(pr$quantiles[, 2] - pr$quantiles[, 1]),
                   prior_width = log10(b$hi) - log10(b$lo),
                   rmse = sqrt(mean((pr$med - truth)^2)))
        }) |>
            list_rbind()
    }) |>
        list_rbind() |>
        mutate(identifiable = coverage90 >= 0.85 & coverage90 <= 0.95 &
                   mean_width < prior_width)
    recov |> mutate(across(where(is.double), \(x) round(x, 3))) |> show_df()
    say("Immigration parameter identifiable (m or I; coverage 0.85-0.95 and ",
        "narrower than the prior): see rows m and I above")
} else {
    say("No held-out simulations found")
}



# =============================================================================*
# ---- 6. posterior-predictive draws ----
# =============================================================================*


say("\n== 6. Posterior-predictive parameter draws")
best <- models[which.max(p_model)]
x <- ref[ref$model == best, ]
wpar <- if ("m" %in% names(x) && any(!is.na(x$m))) "m" else
    if ("I" %in% names(x) && any(!is.na(x$I))) "I" else "scale"
if (!wpar %in% model_pars(x)) wpar <- model_pars(x)[1]
f <- par_fits[[paste(best, wpar)]]
if (is.null(f)) f <- fit_param(x, wpar, stat_names)
d <- data.frame(y = log10(x[[wpar]]), x[, stat_names])
w <- predict(f, obs = as.data.frame(t(obs)), training = d,
             rf.weights = TRUE, paral = TRUE, ncores = N_CORES)$weights
w <- as.numeric(w)
draw <- sample.int(nrow(x), 1000, replace = TRUE, prob = w)
ppc <- x[draw, c("model", model_pars(x))]
say("Model ", best, "; weights from the forest for ", wpar,
    "; effective sample size ", round(1 / sum((w / sum(w))^2)))
write_csv(ppc, file.path(out_dir, "abc-ppc-params.csv"))



# =============================================================================*
# ---- 7. sensitivity ----
# =============================================================================*


say("\n== 7. Sensitivity")
drop_pools <- make.names(c("SN-1-13_S56", "SN-2-97_S37"))
st_keep <- stat_names[!map_lgl(stat_names, \(s) any(str_detect(
    s, fixed(drop_pools))))]
div_stats <- stat_names[str_detect(stat_names, "^(pi|thw|tdr|rare)_")]
tf_stats <- stat_names[str_detect(stat_names, "^(fc|fst)_")]
mc <- bind_rows(
    mc,
    model_choice(ref, st_keep, "without SN-1-13, SN-2-97"),
    model_choice(ref, div_stats, "diversity statistics only"),
    model_choice(ref, tf_stats, "temporal and among-pool statistics only"))
for (v in c("power", "widen", "older", "mu")) {
    xv <- read_sims(v)
    if (is.null(xv)) {
        say("\nVariant ", v, ": not simulated yet")
        next
    }
    mc <- bind_rows(mc, model_choice(clean(syn(xv), v), stat_names,
                                     paste("variant", v)))
}

write_csv(mc, file.path(out_dir, "abc-model-choice.csv"))
write_csv(par_est, file.path(out_dir, "abc-params.csv"))
write_csv(recov, file.path(out_dir, "abc-recovery.csv"))

sink()
read_lines(out_txt) |>
    str_subset("^Growing trees", negate = TRUE) |>
    str_trim("right") |>
    write_lines(out_txt)
