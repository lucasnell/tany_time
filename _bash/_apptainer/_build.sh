#!/bin/bash

#'
#' Script to build the container on the cluster.
#'


#'
#' Copy the files from the repo clone on the cluster (see
#' `_apptainer/README.md`):
#'
#' cd /home/lnell/tany_time && git pull
#' mkdir -p /staging/l/lnell/_apptainer/
#' cp -r _bash/_apptainer/{*.yml,*.def,_build.sh} /staging/l/lnell/_apptainer/
#'



export IN_DEF=tany_time.def
export OUT_SIF=tany_time.sif

mkdir working

mv *.sh *.yml *.def ./working/

cd working

apptainer build ${OUT_SIF} ${IN_DEF}

## If you want to test it:
# apptainer shell -e ${OUT_SIF}

mv ${OUT_SIF} /staging/l/lnell/

cd ..
rm -r working
