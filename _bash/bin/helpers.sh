#'
#' Helper functions used numerous times throughout.
#'
#' Jobs transfer `_bash/bin/` from the repo and source this file after
#' `. /app/.bashrc` (`. bin/helpers.sh`). Sourcing it puts this folder first
#' on `PATH`, so the job runs the repo's scripts. No project scripts are in
#' the container.
#'

_BIN_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
chmod +x "${_BIN_DIR}"/*
export PATH="${_BIN_DIR}:${PATH}"
unset _BIN_DIR


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



#' Exit with given status but only in non-interactive job.
#'
#' Usage:
#' # for error:
#' safe_exit 1
#' # for normal exit:
#' safe_exit 0
#'
safe_exit () {

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
        if [ -f tany_time.sif ]; then rm tany_time.sif; fi
        exit $1
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
        safe_exit $2
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




#' Check that input files exist.
#' Prints every missing file to stderr, then exits (via `safe_exit 111`)
#' if any are missing.
#'
#' Usage:
#' require_files ${FILE1} ${FILE2} ...
#'
require_files () {
    local missing=0
    local f
    for f in "$@"; do
        if [ ! -f "${f}" ]; then
            echo "${f} does not exist! " 1>&2
            missing=1
        fi
    done
    if (( missing != 0 )); then
        safe_exit 111
        return 1
    fi
    return 0
}




#' Copy a gzipped FASTA file to the working directory and decompress it.
#' The output file name is the input's base name without `.gz`.
#' If the second argument is "upper", sequences are converted to uppercase
#' (i.e., softmasking is removed) and header lines are trimmed to their
#' first word.
#'
#' Usage:
#' stage_genome ${GENOME_FULL_PATH}
#' stage_genome ${GENOME_FULL_PATH} upper
#'
stage_genome () {
    local out=$(basename ${1%.gz})
    if [ "$2" == "upper" ]; then
        gunzip -c $1 \
            | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' \
            > ${out}
    else
        gunzip -c $1 > ${out}
    fi
    check_exit_status "stage genome ${out}" $?
    return 0
}




#' Print sequence names from a FASTA file (optionally gzipped) in order,
#' one per line, keeping only the first word of each header.
#' Contig names must not contain spaces.
#'
#' Usage:
#' CONTIG_NAMES=($(contig_names ${GENOME}))
#'
contig_names () {
    if [[ $1 == *.gz ]]; then
        gunzip -c $1
    else
        cat $1
    fi | grep "^>" | sed 's/>//g' | sed 's/[[:space:]].*$//'
}




#' Split a file whose first column is a contig name into one file per
#' contig in the working directory, named `<contig><suffix>`.
#' Exits if any contig isn't in the array `CONTIG_NAMES` (which must be
#' defined beforehand, e.g., using `contig_names`).
#' The optional third argument is an awk expression for what to print
#' (default is the whole line, `$0`).
#'
#' Usage:
#' split_by_contig ${PILEUP} _mp.txt
#' split_by_contig ${SNP} .snp '$2'
#'
split_by_contig () {
    local input=$(readlink -f $1)
    local suffix=$2
    local expr=${3:-'$0'}
    local tmp=$(mktemp -d -p ./)
    local f contig
    (cd ${tmp} && awk "{if (last != \$1) close(last); print ${expr} >> \$1; last = \$1}" ${input})
    check_exit_status "split ${1} by contig" $?
    for f in ${tmp}/*; do
        [ -e "${f}" ] || continue
        contig=$(basename ${f})
        if [[ ! " ${CONTIG_NAMES[*]} " =~ " ${contig} " ]]; then
            echo "Contig ${contig} does not exist in reference! " 1>&2
            safe_exit 1
            return 1
        fi
        mv ${f} ./${contig}${suffix}
    done
    rm -r ${tmp}
    return 0
}
