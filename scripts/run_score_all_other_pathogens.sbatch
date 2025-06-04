#!/usr/bin/bash
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=32000 
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

arr=(Chlamydia-HEK293T CVB3 Dengue-HEK293T Ebola-HEK293T EV71-HEK293T HBV-HUH7 HCV-HEK293T HCV-HUH7 HIV-HEK293T HPV-HEK293T IAV-A549 IAV-HTBE IAV-THP1 KSHV-HEK293T Mtb-HEK293T sars_mers_sars2_human_ppi WNV-HEK293T Zika-HEK293T)

for i in ${!arr[@]}; do 
    echo ${arr[$i]}
    python3 /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/src/score_all.py -o /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/${arr[$i]}_results_onlybestmodel.csv --interface-cutoff 8 --only-best-model  /cluster/project/beltrao/dbaptista/host_pathogen_ppi_struct_pred/results/${arr[$i]}
done


touch ./score_all.done
