#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8000 
#SBATCH --nodes=1
#SBATCH --tmp=20000
#SBATCH -A es_biol
#SBATCH -J score_all
#SBATCH -e ./score_all.err
#SBATCH -o ./score_all.out

source /cluster/home/dbaptista/.bashrc
conda activate vh_struct_pred

python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_af3.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq_af3_notemplates_dockqselectedresidues.csv --interface-cutoff 8 --select-residues-before-dockq /cluster/work/beltrao/jjaenes/25.01.28_af3_virus_no_templates/alphafold3_predictions

python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_af3.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq_af3_notemplates_dockqselectedresidues.csv --interface-cutoff 8 --select-residues-before-dockq /cluster/work/beltrao/jjaenes/25.01.28_af3_bacteria_no_templates/alphafold3_predictions

python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_af3.py -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/hpidb_af3_results.csv --interface-cutoff 8 /cluster/work/beltrao/jjaenes/25.02.06_hpidb_af3_input/alphafold3_predictions

python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_af3.py -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/krogan_lab_datasets_af3_results.csv --interface-cutoff 8 /cluster/work/beltrao/jjaenes/25.02.06_krogan_lab_af3_input/alphafold3_predictions

touch ./score_all.done
