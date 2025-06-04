#!/usr/bin/bash

# Script to run ialign

# Parse in arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -s1|--structure1)
      STRUCT1PATH="$2"
      shift
      shift
      ;;
    -s2|--structure2)
      STRUCT2PATH="$2"
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

if [ -d "${WORKDIR}/tmpfolder" ]; then rm -Rf "${WORKDIR}/tmpfolder"; fi
if [ -d "${WORKDIR}/result" ]; then rm -Rf "${WORKDIR}/result"; fi
if [ -d "${WORKDIR}/result_report" ]; then rm -Rf "${WORKDIR}/result_report"; fi

#cp $STRUCT2PATH "${WORKDIR}/"
#STRUCT2PATH_NEW="${WORKDIR}/$(basename $STRUCT2PATH)"
#IFS='.' read -r NAME EXTENSION <<< "$(basename $STRUCT2PATH)"
#STRUCT2PATH_FIXED="${WORKDIR}/${NAME}_fixed.${EXTENSION}"

#source /cluster/home/dbaptista/.bashrc
#conda activate vh_struct_pred

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/rewrite_pdb.py -f $STRUCT2PATH_NEW

#foldseek easy-complexsearch $STRUCT1PATH $STRUCT2PATH_FIXED ${WORKDIR}/result ${WORKDIR}/tmpfolder --db-extraction-mode 1 --exhaustive-search 1 --format-mode 4 --format-output "query,target,fident,evalue,bits,complexqtmscore,complexttmscore"

foldseek easy-complexsearch $STRUCT1PATH $STRUCT2PATH ${WORKDIR}/result ${WORKDIR}/tmpfolder --db-extraction-mode 1 --exhaustive-search 1 -e inf --format-mode 4 --format-output "query,target,fident,evalue,bits,complexqtmscore,complexttmscore"

#rm -r "tmpfolder"
