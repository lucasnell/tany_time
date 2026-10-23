#!/bin/bash

#' Submit a step with email notification, on the submit server.
#'
#' Usage, from the step's folder (e.g., `_bash/01-variant-calling/`):
#'   ../submit.sh <step>.sub [Complete|Error|Always|Never]
#'
#' Makes `logs/<step>/`, then runs `condor_submit` with `notification` and
#' `notify_user` appended. HTCondor sends one email per job, so the default
#' mode is Complete (an email when each job ends, with its exit code) for
#' steps that queue one job, and Error (an email only for jobs that fail)
#' for steps that queue more. The address is read from NOTIFY_EMAIL (set it
#' in `~/.bashrc` on the submit server; not stored in this public repo).
#' Without NOTIFY_EMAIL, submits as plain `condor_submit` and says so.
#'

SUB=$1
if [ -z "${SUB}" ] || [ ! -f "${SUB}" ]; then
    echo "Usage: $0 <step>.sub [Complete|Error|Always|Never]" 1>&2
    exit 2
fi
STEP=$(basename ${SUB} .sub)
mkdir -p logs/${STEP}

if [ -z "${NOTIFY_EMAIL}" ]; then
    echo "NOTIFY_EMAIL not set; submitting without email" 1>&2
    condor_submit ${SUB}
    exit $?
fi

MODE=$2
if [ -z "${MODE}" ]; then
    # Jobs this file would queue (dry run writes one ad per job)
    N=$(condor_submit -dry-run - ${SUB} 2> /dev/null | grep -c "^ProcId *=")
    if [ "${N}" == "1" ]; then MODE=Complete; else MODE=Error; fi
    echo "${N} job(s) queued: notification = ${MODE}"
fi

condor_submit -append "notification = ${MODE}" \
    -append "notify_user = ${NOTIFY_EMAIL}" ${SUB}
