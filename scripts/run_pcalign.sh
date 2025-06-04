#!/usr/bin/bash

# Script to run PCalign

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
  esac
done

WORKDIR="/cluster/project/beltrao/dbaptista/tools/pcalign"
cd ${WORKDIR}
cp ${STRUCT1PATH} ${WORKDIR}
cp ${STRUCT2PATH} ${WORKDIR}
STRUCT1=$(basename "$STRUCT1PATH")
STRUCT2=$(basename "$STRUCT2PATH")
${WORKDIR}/PCprepare ${STRUCT1}
${WORKDIR}/PCprepare ${STRUCT2}
${WORKDIR}/PCalign ${STRUCT1} ${STRUCT2}
rm *.pdb* *.hash* *.points* *.map*
