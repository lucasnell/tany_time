#!/bin/bash

#' Check consistency of Ω (population covariance) across the thinned
#' subsamples of the main missing-data setting, using the FMD distance
#' (Förstner and Moonen 2003; BayPass `fmd.dist`). Gautier (2015) suggests
#' pairwise FMD < 0.5 indicates consistent estimates. Covariate runs use
#' Ω from all SNPs (`sub0`) regardless, so this is a diagnostic: it reports
#' PASS or WARN but doesn't fail.
#'
#' FMD uses the generalized eigenvalues of (A, B), computed in base R as
#' the eigenvalues of solve(B, A) (same values as `geigen::geigen` on the
#' 2025 Ω matrices).
#'
#' Output: `omega-check.tsv` in `PARENT_DIR_OUT` (and printed).
#'


. /app/.bashrc
. bin/helpers.sh || exit 1
conda activate baypass-env

# Edit these for your system:
export INIT_DIR="/staging/l/lnell/dna/baypass/init"
export PARENT_DIR_OUT="/staging/l/lnell/dna/baypass/init"

# Main missing-data setting (see `_R/make-job-lists.R`):
export MP_TAG="mp075"
export MAX_FMD=0.5

for k in 0 1 2 3; do
    require_files ${INIT_DIR}/baypass-init-${MP_TAG}-sub${k}.tar.gz
done

for k in 0 1 2 3; do
    name=baypass-init-${MP_TAG}-sub${k}
    tar -xzf ${INIT_DIR}/${name}.tar.gz \
        ${name}/tany-init-${MP_TAG}-sub${k}_mat_omega.out
    check_exit_status "extract Ω (sub${k})" $?
    mv ${name}/tany-init-${MP_TAG}-sub${k}_mat_omega.out omega_sub${k}.out
    rm -r ${name}
done

Rscript - << REOF
fmd <- function(a, b) {
    sqrt(sum(log(Re(eigen(solve(b, a), only.values = TRUE)\$values))^2))
}
om <- lapply(0:3, function(k) as.matrix(read.table(sprintf("omega_sub%i.out", k))))
pr <- combn(4, 2)
out <- data.frame(a = sprintf("sub%i", pr[1, ] - 1),
                  b = sprintf("sub%i", pr[2, ] - 1),
                  fmd = apply(pr, 2, function(j) fmd(om[[j[1]]], om[[j[2]]])))
write.table(out, "omega-check.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
subs <- out\$a != "sub0"
if (all(out\$fmd[subs] < ${MAX_FMD})) {
    cat("OMEGA CHECK PASSED: all subsample pairs have FMD <", ${MAX_FMD}, "\n")
} else {
    cat("OMEGA CHECK WARN: some subsample pairs have FMD >=", ${MAX_FMD},
        "(covariate runs use Ω from all SNPs regardless)\n")
}
REOF
check_exit_status "FMD distances" $?

cp omega-check.tsv ${PARENT_DIR_OUT}/
rm omega_sub*.out omega-check.tsv
safe_exit 0
