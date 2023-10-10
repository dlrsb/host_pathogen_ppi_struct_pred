#!/usr/bin/bash

# Script to run DockQ. The predicted structure PDB file is renumbered beforehand.

SORTNATIVE=False

# Parse in arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--model)
      MODELPATH="$2"
      shift
      shift
      ;;
    -n|--native)
      NATIVEPATH="$2"
      shift
      shift
      ;;
    --native_chain1)
      NATIVECH1="$2"
      shift
      shift
      ;;
    --model_chain1)
      MODELCH1="$2"
      shift
      shift
      ;;
    --native_chain2)
      NATIVECH2="$2"
      shift
      shift
      ;;
    --model_chain2)
      MODELCH2="$2"
      shift
      shift
      ;;
	--sort_native)
      SORTNATIVE=True
      shift
      shift
      ;;
  esac
done

source /cluster/home/dbaptista/.bashrc
conda activate vh_struct_pred
echo "$PWD"
wdir="$PWD"
if [ "$SORTNATIVE" = True ]; then
	DIR=$(dirname "$NATIVEPATH")
	FILENAME=$(basename "$NATIVEPATH")
	cd ${DIR}
	echo "$PWD"
	pdb_splitmodel ${NATIVEPATH}
	IFS='.' read -r NAME EXTENSION <<< "$FILENAME"
	SPLITFILE="${DIR}/${NAME}_1.${EXTENSION}"
	SORTEDFILE="${DIR}/${NAME}_sorted.${EXTENSION}"
	echo "$SPLITFILE"
	echo "$SORTEDFILE"
	pdb_sort -C ${SPLITFILE} > ${SORTEDFILE}
	rm ${SPLITFILE}
	NATIVEPATH=${SORTEDFILE}
	echo "$NATIVEPATH"
	cd ${wdir}
	echo "$PWD"
fi

pdb_delinsertion ${NATIVEPATH} | pdb_delhetatm > ${NATIVEPATH}.modified
conda deactivate

source /cluster/apps/local/env2lmod.sh
module load gcc/6.3.0 python/3.8.5 emboss/6.6.0

# Renumber chains in both model and native and then fix numbering in model by aligning its sequences with the native sequences
# DockQ does this too, but only when len(model_chains) > 2 or len(native_chains)> 2
perl /cluster/project/beltrao/dbaptista/tools/DockQ/scripts/renumber_pdb.pl ${MODELPATH}

perl /cluster/project/beltrao/dbaptista/tools/DockQ/scripts/renumber_pdb.pl ${NATIVEPATH}.modified

perl /cluster/project/beltrao/dbaptista/tools/DockQ/scripts/fix_numbering.pl ${MODELPATH}.renum ${NATIVEPATH}.modified.renum

python /cluster/project/beltrao/dbaptista/tools/DockQ/DockQ.py ${MODELPATH}.renum.fixed ${NATIVEPATH}.modified.renum -short -model_chain1 ${MODELCH1} -model_chain2 ${MODELCH2} -native_chain1 ${NATIVECH1} -native_chain2 ${NATIVECH2}

rm ${NATIVEPATH}.modified ${MODELPATH}.renum ${NATIVEPATH}.modified.renum ${MODELPATH}.renum.fixed