#!/usr/bin/bash

# Script to run ialign

# Parse in arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -s1|--structure1)
      STRUCT1PATH="$2"
      shift;
      shift;
      ;;
    -s2|--structure2)
      STRUCT2PATH="$2"
      shift;
      shift;
      ;;
    -n|--normalization)
      NORMALIZATION="$2"
      shift;
      shift;
      ;;
    -m|--scoring_metric)
      SCORINGMETRIC="$2"
      shift;
      shift;
      ;;
    --dc)
      DISTANCECUTOFF="$2"
      shift;
      shift;
      ;;
    --minp)
      MINP="$2"
      shift;
      shift;
      ;;
    --mini)
      MINI="$2"
      shift;
      shift;
      ;;
  esac
done

echo "$PWD"


TMPDIR="/cluster/scratch/dbaptista/ialign_tmpdir"
echo "$TMPDIR"
cd $TMPDIR
if [ -d "outputs" ]; then rm -Rf outputs; fi
mkdir outputs

echo "$PWD"

cp $STRUCT1PATH ${TMPDIR}/
STRUCT1PATH_NEW="${TMPDIR}/$(basename $STRUCT1PATH)"
IFS='.' read -r NAME1 EXTENSION1 <<< "$(basename $STRUCT1PATH)"
STRUCT1PATH_FIXED="${TMPDIR}/${NAME1}_fixed.${EXTENSION1}"

cp $STRUCT2PATH ${TMPDIR}/
STRUCT2PATH_NEW="${TMPDIR}/$(basename $STRUCT2PATH)"
IFS='.' read -r NAME2 EXTENSION2 <<< "$(basename $STRUCT2PATH)"
STRUCT2PATH_FIXED="${TMPDIR}/${NAME2}_fixed.${EXTENSION2}"

source /cluster/home/dbaptista/.bashrc
conda activate vh_struct_pred

python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/rewrite_pdb.py -f ${STRUCT1PATH_NEW} -o ${STRUCT1PATH_FIXED}
python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/rewrite_pdb.py -f ${STRUCT2PATH_NEW} -o ${STRUCT2PATH_FIXED}

perl /cluster/project/beltrao/dbaptista/tools/ialign/bin/ialign.pl -s -a 0 -w outputs -o outputs/results.txt -e ${SCORINGMETRIC} -n ${NORMALIZATION} -dc ${DISTANCECUTOFF} -minp ${MINP} -mini ${MINI} -p1 ${STRUCT1PATH_FIXED} -p2 ${STRUCT2PATH_FIXED} > ${TMPDIR}/stdout.txt

cat outputs/results.txt
rm *pdb *txt
