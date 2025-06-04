#!/usr/bin/bash

# Script to run ialign

# Parse in arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -d|--structure-dir)
      STRUCTDIR="$2"
      shift
      shift
      ;;
    --workdir)
      WORKDIR="$2"
      shift
      shift
      ;;
  esac
done

if [ -d "${WORKDIR}/tmpfolder" ]; then rm -rf "${WORKDIR}/tmpfolder"; fi
if [ -f "${WORKDIR}/result" ]; then rm "${WORKDIR}/result"; fi
if [ -f "${WORKDIR}/result_report" ]; then rm "${WORKDIR}/result_report"; fi

#foldseek createdb ${STRUCTDIR} ${WORKDIR}/af2_models_db --db-extraction-mode 1

foldseek easy-complexsearch ${WORKDIR}/af2_models_db ${WORKDIR}/af2_models_db ${WORKDIR}/result ${WORKDIR}/tmpfolder -e inf --db-extraction-mode 1 --exhaustive-search 1 --prefilter-mode 2 --tmscore-threshold 0.0 --remove-tmp-files 1 --compressed 1 --format-mode 4 --format-output "query,target,fident,evalue,bits,complexqtmscore,complexttmscore,lddt,qtmscore,ttmscore,alntmscore,rmsd,prob"
