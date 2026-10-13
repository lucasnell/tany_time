#!/bin/bash

#' Build a Kraken2 database from SILVA SSU rRNA (SSURef NR99; Bacteria,
#' Archaea, and Eukarya), for detecting eukaryotic parasites (fungi,
#' oomycetes, microsporidia, nematodes, mites, ...) in unmapped reads
#' (`06b-kraken`). Follows Kraken2's own `16S_silva_installation.sh`, but
#' with SILVA 138.2 (Kraken2 2.17.2 uses 138.1).
#'
#' Compute nodes may lack internet access, so first download these files
#' into `SILVA_DIR` on staging (e.g., with `wget` on the submit server):
#'   https://ftp.arb-silva.de/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
#'   https://ftp.arb-silva.de/release_138_2/Exports/taxonomy/tax_slv_ssu_138.2.acc_taxid.gz
#'   https://ftp.arb-silva.de/release_138_2/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz
#'
#' Output: `${DB_NAME}.tar.gz` in `PARENT_DIR_OUT`.
#'


. /app/.bashrc
conda activate kraken-env

export THREADS=$(count_threads)

# Edit these for your system:
export SILVA_DIR="/staging/lnell/silva"
export PARENT_DIR_OUT="/staging/lnell/kraken2-db"

export SILVA_VERSION="138.2"
export DB_NAME="silva-${SILVA_VERSION}"
export FASTA="SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva.fasta"
export TAXO_PREFIX="tax_slv_ssu_${SILVA_VERSION}"

require_files ${SILVA_DIR}/${FASTA}.gz \
    ${SILVA_DIR}/${TAXO_PREFIX}.acc_taxid.gz \
    ${SILVA_DIR}/${TAXO_PREFIX}.txt.gz

# Script from the Kraken2 install that converts SILVA taxonomy to
# names.dmp/nodes.dmp:
export TAX_SCRIPT=$(find ${CONDA_PREFIX} -name build_silva_taxonomy.pl | head -n 1)
require_files ${TAX_SCRIPT}

mkdir -p ${DB_NAME}/data ${DB_NAME}/taxonomy ${DB_NAME}/library
cd ${DB_NAME}/data
for f in ${FASTA} ${TAXO_PREFIX}.acc_taxid ${TAXO_PREFIX}.txt; do
    gunzip -c ${SILVA_DIR}/${f}.gz > ${f}
    check_exit_status "gunzip ${f}" $?
done
perl ${TAX_SCRIPT} ${TAXO_PREFIX}.txt
check_exit_status "build_silva_taxonomy.pl" $?
cd ..
mv data/names.dmp data/nodes.dmp taxonomy/
mv data/${TAXO_PREFIX}.acc_taxid seqid2taxid.map
# RNA to DNA alphabet in sequences:
sed -e '/^>/!y/U/T/' data/${FASTA} > library/silva.fna
check_exit_status "convert sequences" $?
rm -r data
cd ..

kraken2-build --db ${DB_NAME} --build --threads ${THREADS}
check_exit_status "kraken2-build" $?
kraken2-build --db ${DB_NAME} --clean
check_exit_status "kraken2-build --clean" $?

tar -czf ${DB_NAME}.tar.gz ${DB_NAME}
mkdir -p ${PARENT_DIR_OUT}
mv ${DB_NAME}.tar.gz ${PARENT_DIR_OUT}/
rm -r ${DB_NAME}

safe_exit 0
