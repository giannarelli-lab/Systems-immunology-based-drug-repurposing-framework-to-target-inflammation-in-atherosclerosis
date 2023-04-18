#!/bin/bash -l

#SBATCH -J merge
#SBATCH --mem=300000
#SBATCH -n 1
#SBATCH --array=50
#SBATCH -p fn_medium
#SBATCH --export=ALL
#SBATCH --time=1200
#SBATCH --output="merge_data-%A_%a.out"

#date
d1=$(date +%s)

##### merge files as a job
path='/gpfs/Data/Saracatinib_allsamples'

lane1=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_merge.txt | cut -d' ' -f 1`
lane2=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_merge.txt | cut -d' ' -f 2`
merged=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_merge.txt | cut -d' ' -f 3`


cat $path/fastq/$lane1.fastq.gz $path/fastq/$lane2.fastq.gz > $path/merged_data/$merged.fastq.gz

echo $merged

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)