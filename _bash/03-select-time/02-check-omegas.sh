#!/bin/bash

#' Check consistency of Ω (population covariance) across the 3 thinned
#' subsamples of each missing-data setting, using the FMD distance
#' (Förstner and Moonen 2003; BayPass `fmd.dist`). Gautier (2015) suggests
#' pairwise FMD < 0.5 indicates consistent estimates. Each covariate run
#' uses its own subsample's Ω, so this is a diagnostic: it reports PASS or
#' WARN but doesn't fail. If subsample Ωs disagree, covariate results may
#' differ among subsamples for that reason alone.
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

export MP_TAGS=(mp050 mp075 mp100)
export MAX_FMD=0.5

for t in ${MP_TAGS[@]}; do
    for k in 1 2 3; do
        require_files ${INIT_DIR}/baypass-init-${t}-sub${k}.tar.gz
    done
done

for t in ${MP_TAGS[@]}; do
    for k in 1 2 3; do
        name=baypass-init-${t}-sub${k}
        tar -xzf ${INIT_DIR}/${name}.tar.gz \
            ${name}/tany-init-${t}-sub${k}_mat_omega.out
        check_exit_status "extract Ω (${t} sub${k})" $?
        mv ${name}/tany-init-${t}-sub${k}_mat_omega.out omega_${t}_sub${k}.out
        rm -r ${name}
    done
done

Rscript - << REOF
fmd <- function(a, b) {
    sqrt(sum(log(Re(eigen(solve(b, a), only.values = TRUE)\$values))^2))
}
out <- do.call(rbind, lapply(c("mp050", "mp075", "mp100"), function(t) {
    om <- lapply(1:3, function(k) {
        as.matrix(read.table(sprintf("omega_%s_sub%i.out", t, k)))
    })
    pr <- combn(3, 2)
    data.frame(tag = t, a = sprintf("sub%i", pr[1, ]),
               b = sprintf("sub%i", pr[2, ]),
               fmd = apply(pr, 2, function(j) fmd(om[[j[1]]], om[[j[2]]])))
}))
write.table(out, "omega-check.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
if (all(out\$fmd < ${MAX_FMD})) {
    cat("OMEGA CHECK PASSED: all subsample pairs have FMD <", ${MAX_FMD}, "\n")
} else {
    cat("OMEGA CHECK WARN: some subsample pairs have FMD >=", ${MAX_FMD},
        "(each covariate run uses its own subsample's Ω)\n")
}
REOF
check_exit_status "FMD distances" $?

cp omega-check.tsv ${PARENT_DIR_OUT}/
rm omega_*_sub*.out omega-check.tsv
safe_exit 0
