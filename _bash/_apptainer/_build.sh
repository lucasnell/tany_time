#!/bin/bash

#'
#' Script to build the container on the cluster.
#'


export IN_DEF=tany_time.def
export OUT_SIF=tany_time.sif

mkdir working

mv *.sh *.yml *.py *.def ./working/

cd working

apptainer build ${OUT_SIF} ${IN_DEF}

## If you want to test it:
# apptainer shell -e ${OUT_SIF}

mv ${OUT_SIF} /staging/lnell/

cd ..
rm -r working
