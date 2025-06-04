#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8000 
#SBATCH --nodes=1
#SBATCH --tmp=20000
#SBATCH -A es_biol
#SBATCH -J score_all
#SBATCH -e ./score_all.err
#SBATCH -o ./score_all.out

source /cluster/apps/local/env2lmod.sh
module load gcc/6.3.0

source /cluster/home/dbaptista/.bashrc
conda activate vh_struct_pred

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq.csv --interface-cutoff 8 --rename-chains  /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq_onlybestmodel.csv --interface-cutoff 8 --rename-chains --only-best-model  /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq_onlybestmodel.csv --interface-cutoff 8 --rename-chains  --only-best-model /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq_5models.csv --interface-cutoff 8 --rename-chains  /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq_5models.csv --interface-cutoff 8 --rename-chains  /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues.csv --interface-cutoff 8 --rename-chains --only-best-model  /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues.csv --interface-cutoff 8 --rename-chains  --only-best-model /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv --interface-cutoff 8 --rename-chains --only-best-model --select-residues-before-dockq /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq

#python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv --interface-cutoff 8 --rename-chains --only-best-model --select-residues-before-dockq /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq

python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/virus_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv --interface-cutoff 8 --rename-chains --only-best-model --select-residues-before-dockq /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/virus_mammalia_dimers_fullseq

python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -n /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/native_pdb_files -f /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/bacteria_mammalia_dimers/fullseq_fasta_files/protein_pairs -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq_onlybestmodel_dockqselectedresidues_recalc.csv --interface-cutoff 8 --rename-chains --only-best-model --select-residues-before-dockq /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/bacteria_mammalia_dimers_fullseq

touch ./score_all.done
