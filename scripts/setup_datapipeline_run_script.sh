#!/bin/bash

# Initialize variables
FASTAFILE="undefined.fasta"
WORKDIR=$PWD
MSADIR="${PWD}/msas"
MAX_TEMPLATE_DATE=$(date +'%Y-%m-%d')
MULTIMER=False
BATCH_SYS="SLURM"

print_help()
{
   # Display Help
   echo "Script to create sbatch script to run the AlphaFold2 DataPipeline on Euler."
   echo
   echo "Syntax: setup_datapipeline_run_script.sh [-f fastafile] [-w working directory] [-s shareholder group] [-o output directory] [--multimer][--max_template_date Y-M-D] [--reduced_dbs]"
   echo "options:"
   echo "-h                     print help and exit"
   echo "-f                     FASTA filename"
   echo "-s                     shareholder group for the use of GPUs. Mandatory for the submissions of scripts with SLURM"
   echo "-w                     working directory"
   echo "-o                     output directory (where MSAs will be stored)"
   echo "--multimer             use the multimer pipeline to create files that will be need for prediction with AlphaFold-multimer"
   echo "--reduced_dbs          use settings with reduced hardware requirements"
   echo "--max_template_date    format: "
   echo
   exit 1
}

# Print help if not options are provided
if [[ $# -eq 0 ]];then
    print_help
    exit 1
fi

REDUCED_DBS=False

# Parse in arguments
while [[ $# -gt 0 ]]; do
    case $1 in
	 -h|--help)
          # Print help and exit
          print_help
          exit
	  ;;
        -f|--fastafile)
          # Get absolute path
          FASTAFILE=$(readlink -f $2)
          # Get the protein name
          fastaname=$(basename -- "$FASTAFILE")
          PROTEIN="${fastaname%.*}"
          echo "  Reading $FASTAFILE"
          echo "  Protein name:              $PROTEIN"
          shift;
          shift;
          ;;
        -w|--workdir)
          # Users can specify a work directory, e.g., $SCRATCH/alphafold_tests
          # Otherwise it will use the current directy as a work directory
          WORKDIR="$2"
          shift;
          shift;
          ;;
        -s|--shareholder)
          # For the submission of SLURM jobs, the shareholder group is mandatory
          SHAREHOLDER_GROUP="$2"
          shift;
          shift;
          ;;
        -o|--output)
          # Users can specify a work directory, e.g., $SCRATCH/alphafold_tests
          # Otherwise it will use the current directy as a work directory
          MSADIR="$2"
          shift;
          shift;
          ;;
        --max_template_date)
          # The max template date of the databases to use for pair representation
          # This could affect the accuracy of the outcome
          MAX_TEMPLATE_DATE="$2"
          shift;
          shift;
          ;;
        --multimer)
          MULTIMER=True
          shift;
          ;;
        --reduced_dbs)
          REDUCED_DBS=True
          shift;
          ;;
        * )
          print_help
          exit 1
    esac
done


if [[ $BATCH_SYS = "SLURM" &&  $SHAREHOLDER_GROUP = "" ]]; then
        echo
        echo -e "Please provide your shareholder group with the -s option"
        echo -e "This parameter is mandatory when requesting GPUs with SLURM"
        echo -e "You can display all the groups you are a part of on Euler using the my_share_info command"
        echo
        print_help
fi

if [ "$MULTIMER" = True ]; then
    echo "  Protein type:              multimer"
    OPTIONS="--model_preset=multimer --pdb_seqres_database_path=\$DATA_DIR/pdb_seqres/pdb_seqres.txt --uniprot_database_path=\$DATA_DIR/uniprot/uniprot.fasta \\"$'\n'
else
    echo "  Protein type:              monomer"
    OPTIONS="--model_preset=monomer --pdb70_database_path=\$DATA_DIR/pdb70/pdb70 \\"$'\n'
fi

if [ "$REDUCED_DBS" = True ]; then
    OPTIONS+="--db_preset=reduced_dbs \\"$'\n'
    OPTIONS+="--small_bfd_database_path=\$DATA_DIR/small_bfd/bfd-first_non_consensus_sequences.fasta \\"$'\n'
else
    OPTIONS+="--db_preset=full_dbs \\"$'\n'
    OPTIONS+="--bfd_database_path=\$DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\"$'\n'
    OPTIONS+="--uniref30_database_path=\$DATA_DIR/uniref30/UniRef30_2021_03 \\"$'\n'
fi

# Determine the sequence length
# The required total GPU mem depends on the sum of the number of amino acids
# The required total CPU mem depends on the max of the number of amino acids
sum_aa=$(cat $FASTAFILE | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' | awk ' { getline aa; sum+=length(aa); } END { print sum } ')
max_aa=$(cat $FASTAFILE | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' | awk ' BEGIN {max=0} { getline aa; if (length(aa) > max) {max=length(aa)}} END { print max } ')
echo "  Number of amino acids:"
echo "                    sum:     $sum_aa"
echo "                    max:     $max_aa"

# Estimate the required computing resources
NCPUS=8
if (( "$sum_aa" < 200 )); then
    RUNTIME="24:00" #"04:00"
    TOTAL_CPU_MEM_MB=120000
    TOTAL_SCRATCH_MB=120000
elif (( "$sum_aa" >= 200 )) && (( "$sum_aa" < 1500 )); then
    RUNTIME="24:00"
    TOTAL_CPU_MEM_MB=120000
    TOTAL_SCRATCH_MB=120000
elif (( "$sum_aa" >= 1500 )) && (( "$sum_aa" < 2500 )); then
    RUNTIME="120:00" #"24:00"
    TOTAL_CPU_MEM_MB=240000
    TOTAL_SCRATCH_MB=240000
elif (( "$sum_aa" >= 2500 )) && (( "$sum_aa" < 3500 )); then
    RUNTIME="120:00" #"48:00"
    TOTAL_CPU_MEM_MB=480000
    TOTAL_SCRATCH_MB=240000
elif (( "$sum_aa" >= 3500 )); then
    RUNTIME="120:00"
    TOTAL_CPU_MEM_MB=480000 #640000
    TOTAL_SCRATCH_MB=320000
fi

echo -e "    Estimate required resources, please do not hesitate to adjust if required: "
echo -e "    Run time:            " $RUNTIME
echo -e "    Number of CPUs:      " $NCPUS
echo -e "    Total CPU memory:    " $TOTAL_CPU_MEM_MB
echo -e "    Total scratch space: " $TOTAL_SCRATCH_MB

########################################
# Output a SLURM run script
########################################

mkdir -p $WORKDIR
RUNSCRIPT=$WORKDIR/"af2_msas_$PROTEIN.sbatch"
echo -e "  Output a SLURM run script to run the AlphaFold2 DataPipeline: $RUNSCRIPT"

RUNTIME="${RUNTIME}":00" "

cat <<EOF > $RUNSCRIPT
#!/usr/bin/bash
#SBATCH -n $NCPUS
#SBATCH --time=$RUNTIME
#SBATCH --mem-per-cpu=$((TOTAL_CPU_MEM_MB/NCPUS))
#SBATCH --ntasks-per-node=$NCPUS
#SBATCH --nodes=1
#SBATCH --tmp=$TOTAL_SCRATCH_MB
#SBATCH -A $SHAREHOLDER_GROUP
#SBATCH --partition $SHAREHOLDER_GROUP
#SBATCH -J af2_msas_$PROTEIN
#SBATCH -e $WORKDIR/af2_msas_$PROTEIN.err.txt
#SBATCH -o $WORKDIR/af2_msas_$PROTEIN.out.txt

set -e

source /cluster/apps/local/env2lmod.sh
module load gcc/6.3.0 openmpi/4.0.2 alphafold/2.3.1
source /cluster/apps/nss/alphafold/venv_alphafold_2.3.1/bin/activate

# Define paths to databases and output directory
DATA_DIR=/cluster/project/alphafold
OUTPUT_DIR=\${TMPDIR}/output

python /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/run_alphafold_datapipeline.py \\
--data_dir=\$DATA_DIR \\
--output_dir=\$OUTPUT_DIR \\
--max_template_date="$MAX_TEMPLATE_DATE" \\
--uniref90_database_path=\$DATA_DIR/uniref90/uniref90.fasta \\
--mgnify_database_path=\$DATA_DIR/mgnify/mgy_clusters_2022_05.fa \\
--template_mmcif_dir=\$DATA_DIR/pdb_mmcif/mmcif_files \\
--obsolete_pdbs_path=\$DATA_DIR/pdb_mmcif/obsolete.1.dat \\
$OPTIONS --fasta_paths=$FASTAFILE

rsync -av \$TMPDIR/output/$PROTEIN $MSADIR

touch $WORKDIR/af2_msas_$PROTEIN.done

EOF
