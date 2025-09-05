#!/bin/bash

#'
#' Script to build the container on the cluster.
#'


#'
#' If necessary, to send all these scripts to the cluster from your machine:
#'
#' cd ~/GitHub/Wisconsin/tany_time/_bash/_apptainer
#' scp *.yml helpers.sh *.def *.py lnell@ap2001.chtc.wisc.edu:/home/lnell/_apptainer/
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
