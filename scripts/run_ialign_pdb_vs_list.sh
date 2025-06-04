#!/usr/bin/bash

# Script to run ialign

# Parse in arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -p|--pdb_file)
      PDBFILEPATH="$2"
      shift
      shift
      ;;
    -l|--pdb_list)
      PDBLISTPATH="$2"
      shift
      shift
      ;;
    -w|--workdir)
      WORKDIR="$2"
      shift
      shift
      ;;
    -n|--normalization)
      NORMALIZATION="$2"
      shift
      shift
      ;;
    -m|--scoring_metric)
      SCORINGMETRIC="$2"
      shift
      shift
      ;;
    --dc)
      DISTANCECUTOFF="$2"
      shift
      shift
      ;;
    --minp)
      MINP="$2"
      shift
      shift
      ;;
    --mini)
      MINI="$2"
      shift
      shift
      ;;
  esac
done

TMPDIR="/cluster/scratch/dbaptista/ialign_tmpdir"
cd $TMPDIR
if [ -d "${WORKDIR}" ]; then rm -Rf ${WORKDIR}; fi
mkdir ${WORKDIR}

perl /cluster/project/beltrao/dbaptista/tools/ialign/bin/ialign.pl -s -a 0 -w ${WORKDIR} -o ${WORKDIR}/results.txt -e ${SCORINGMETRIC} -n ${NORMALIZATION} -dc ${DISTANCECUTOFF} -minp ${MINP} -mini ${MINI} -l ${PDBLISTPATH} > ${TMPDIR}/stdout.txt

cat ${WORKDIR}/results.txt
#rm *pdb *txt *lst
