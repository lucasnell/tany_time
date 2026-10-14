#!/bin/bash

#' Run once on the submit server (not a job): put Dorado 0.9.6 and the
#' R9.4.1 sup basecalling model on staging for `02-dorado`.
#'
#' Dorado 0.9.6 is the last release that reads FAST5 and has R9.4.1 models
#' (both removed in Dorado 1.0.0; see Dorado's CHANGELOG). Needs internet.
#'
#' Output (`/staging/l/lnell/dorado/`): `dorado-0.9.6-linux-x64.tar.gz`,
#' `dna_r9.4.1_e8_sup@v3.6.tar.gz`.
#'

set -e

OUT=/staging/l/lnell/dorado
DORADO=dorado-0.9.6-linux-x64
MODEL=dna_r9.4.1_e8_sup@v3.6

mkdir -p ${OUT}
TMP=$(mktemp -d)
cd ${TMP}
curl -sSL -O https://cdn.oxfordnanoportal.com/software/analysis/${DORADO}.tar.gz
tar -xzf ${DORADO}.tar.gz
./${DORADO}/bin/dorado --version
./${DORADO}/bin/dorado download --model ${MODEL}
tar -czf ${MODEL}.tar.gz ${MODEL}
cp ${DORADO}.tar.gz ${MODEL}.tar.gz ${OUT}/
ls -l ${OUT}
cd /
rm -rf ${TMP}
