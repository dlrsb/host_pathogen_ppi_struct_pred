#!/bin/bash

# Script to download files from RCSB http file download services. Modified.
# Use the -h switch to get help on usage.

if ! command -v curl &> /dev/null
then
    echo "'curl' could not be found. You need to install 'curl' for this script to work."
    exit 1
fi

PROGNAME=$0
BASE_URL="https://files.rcsb.org/download"

usage() {
  cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]

 -f <file>: the input file containing a comma-separated list of PDB ids
 -o  <dir>: the output dir, default: current dir
EOF
  exit 1
}

download() {
  url="$BASE_URL/$1"
  out=$2/$1
  IFS=. read -r entry pdb extension <<< "$1"
  assembly="${pdb: -1}"
  echo "Downloading $url to $out"
  curl -s -f $url -o $out
  if [ $? -eq 0 ]; then 
	gzip -d $out
	mv $2/${entry}.pdb${assembly} $2/${entry}-${assembly}.pdb
  else
	echo "Failed to download $url"
	ciffile=${entry}-assembly${assembly}.cif.gz
	newurl="$BASE_URL/${ciffile}"
	newout="$2/${ciffile}"
	curl -s -f $newurl -o $newout || echo "Failed to download $newurl"
	gzip -d $newout
	/cluster/project/beltrao/dbaptista/tools/maxit-v11.100-prod-src/bin/maxit -i $2/${entry}-assembly${assembly}.cif -o 2
	mv $PWD/${entry}-assembly${assembly}.cif.pdb $2/${entry}-${assembly}.pdb
  fi
}

listfile=""
outdir="."
while getopts f:o:cpaxsmr o
do
  case $o in
    (f) listfile=$OPTARG;;
    (o) outdir=$OPTARG;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"

if [ "$listfile" == "" ]
then
  echo "Parameter -f must be provided"
  exit 1
fi
contents=$(cat $listfile)

mkdir -p $outdir/native_pdb_files
mkdir -p $outdir/fasta_files

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"

for token in "${tokens[@]}"
do
	IFS=- read -r name bioassembly <<< "$token"
	if [ ! -f $outdir/native_pdb_files/${name}-${bioassembly}.pdb ]; then
		download ${name}.pdb${bioassembly}.gz $outdir/native_pdb_files
		
		python src/extract_sequences.py -i $outdir/native_pdb_files/${name}-${bioassembly}.pdb -o $outdir/fasta_files
		
		# preprocess file using pdb_tools
		# pdb_delhetatm ${name}-${bioassembly}.pdb | pdb_delelem -H | pdb_reres > ${name}-${bioassembly}_modified.pdb
	fi

done








