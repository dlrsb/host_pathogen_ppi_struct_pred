#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=12000 
#SBATCH --nodes=1
#SBATCH --tmp=20000
#SBATCH -A es_biol
#SBATCH -J mmseqs_cluster
#SBATCH -e ./mmseqs_cluster.err
#SBATCH -o ./mmseqs_cluster.out

module load stack/2024-06 gcc/12.2.0 mmseqs2/14-7e284

rm -r $SCRATCH/tmp
mkdir $SCRATCH/tmp

mmseqs easy-linclust /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/all_seqs_min10aa_filenames_as_ids.fasta $SCRATCH/clusterRes70 $SCRATCH/tmp --cluster-mode 1 -c 0.80 --min-seq-id 0.7

mmseqs easy-linclust /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/data/all_seqs_min10aa_filenames_as_ids.fasta $SCRATCH/clusterRes95cov $SCRATCH/tmp --cluster-mode 1 -c 0.90 --min-seq-id 0.95

touch ./mmseqs_cluster.done
