#'
#' Helper functions used numerous times throughout:
#'


#' Find the number of threads on CHTC cluster and print to stdout.
#' Apparently sometimes `$_CONDOR_MACHINE_AD` isn't present, so based on some
#' testing I found a number of other objects to use.
#' If none are found (and produce integers greater than 0), use `nproc`.
#'
#' Usage:
#' count_threads
#' or to define an object:
#' THREADS=$(count_threads)
#'
count_threads () {
    local nt=""
    if [ ! -z "${_CONDOR_MACHINE_AD}" ] && [ -f "${_CONDOR_MACHINE_AD}" ]; then
        nt=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')
    fi
    if [ ! -z ${nt} ] && [[ ${nt} =~ ^[0-9]+$ ]] && (( nt >= 1 )); then
        echo ${nt}
    elif [ ! -z ${OMP_NUM_THREADS} ] && [[ ${OMP_NUM_THREADS} =~ ^[0-9]+$ ]] &&
        (( OMP_NUM_THREADS >= 1 )); then
        echo ${OMP_NUM_THREADS}
    elif [ ! -z ${JULIA_NUM_THREADS} ] && [[ ${JULIA_NUM_THREADS} =~ ^[0-9]+$ ]] &&
        (( JULIA_NUM_THREADS >= 1 )); then
        echo ${JULIA_NUM_THREADS}
    elif [ ! -z ${MKL_NUM_THREADS} ] && [[ ${MKL_NUM_THREADS} =~ ^[0-9]+$ ]] &&
        (( MKL_NUM_THREADS >= 1 )); then
        echo ${MKL_NUM_THREADS}
    elif [ ! -z ${NUMEXPR_NUM_THREADS} ] && [[ ${NUMEXPR_NUM_THREADS} =~ ^[0-9]+$ ]] &&
        (( NUMEXPR_NUM_THREADS >= 1 )); then
        echo ${NUMEXPR_NUM_THREADS}
    elif [ ! -z ${OPENBLAS_NUM_THREADS} ] && [[ ${OPENBLAS_NUM_THREADS} =~ ^[0-9]+$ ]] &&
        (( OPENBLAS_NUM_THREADS >= 1 )); then
        echo ${OPENBLAS_NUM_THREADS}
    else
        nproc
    fi

    return 0
}




#' Check previous command's exit status.
#' If != 0, then...
#' ... non-interactive session: send error message, print any *.log or
#'                              *.err files, and exit.
#' ... interactive session: just break and ignore errors if the context
#'                          isn't sensible.
#'
#' Usage:
#' OperationX ... options ...
#' check_exit_status "OperationX" $?
#'
check_exit_status () {
    if [ "$2" != "0" ]; then
        echo "Step $1 failed with exit status $2" 1>&2
        local LOGS=$(ls *.err *.out 2> /dev/null)
        for f in $LOGS; do
            echo -e "\n\nFILE: " ${f} "\n\n" 1>&2
            cat ${f} 1>&2
        done
        out_loc=$(readlink -f /proc/self/fd/0)
        if [ "$out_loc" != "/dev/null" ]; then
            # (interactive shell)
            break 2> /dev/null
        else
            # (non-interactive shell)
            if [ -z "$_CONDOR_SCRATCH_DIR" ]; then
                echo -n "No condor scratch directory detected. " 1>&2
                echo "No files removed." 1>&2
            else
                cd $_CONDOR_SCRATCH_DIR
                # Remove everything except for the files / folders you start with:
                ## shopt -s extglob
                ## TO_RM=$(ls -d !(*condor*|tmp|var|*docker*) 2> /dev/null)
                TO_RM=$(find . \! \( -name \*condor\* -o -name tmp -o -name var -o -name \*docker\* \))
                rm -rf ${TO_RM} 2> /dev/null
                ## shopt -u extglob
            fi
            exit $2
        fi
        return 0
    fi
    if [[ "$1" != "null" ]]; then
        echo "Checked step $1"
    fi
    return 0

}





#' Get a read name from a tar file.
#'
#' Usage for 1st of pair:
#' read_tar_name ${TAR_FILE} 1
#' Usage for 2nd of pair:
#' read_tar_name ${TAR_FILE} 2
read_tar_name () {
    if (( $# != 2 )); then
        echo -e "read_tar_names requires 2 args" 1>&2
        return 1
    fi
    if (( $2 != 1 )) && (( $2 != 2 )); then
        echo -e "read_tar_names 2nd arg must be 1 or 2" 1>&2
        return 1
    fi
    local READS_ARR=($(tar -tf $1))
    local SORTED_READS=""
    IFS=$'\n' SORTED_READS=($(sort <<<"${READS_ARR[*]}")); unset IFS
    echo ${SORTED_READS[$(( $2 - 1 ))]}
    return 0
}




#' Check on BAM file with `bamtools stats`, check status of this call,
#' then output the file name.
#' THIS MUST BE RUN INSIDE `main-env`!
#'
#' Usage:
#' call_bam_stats ${BAM_FILE} "Label for file"
#'
call_bam_stats () {
    local B=$1
    local S=${B/.bam/.stats}
    bamtools stats -in $B | tee $S
    check_exit_status "bamtools stats $2" $?
    echo -e "FILE:" $B "\n**********************************************"
}
