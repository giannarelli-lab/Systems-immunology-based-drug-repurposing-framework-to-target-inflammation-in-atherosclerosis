#!/bin/bash -l

#SBATCH -J fastqc
#SBATCH --mem=250000
#SBATCH -n 8
#SBATCH --array=1-36
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=600
#SBATCH --output="fastqc1_raw_data-%A_%a.out"

#date
d1=$(date +%s)

module load fastqc/0.11.7


##### run fastqc in raw reads ####
path='/gpfs/Data/2022-06-17/fastq'
path1='/gpfs/Data/RNAseq/06102022'

sample=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path1/raw_data_sample_names.txt`
input_file=$path/$sample.fastq.gz
output_folder=$path1/fastqc_raw_data/

echo $HOSTNAME
echo $sample
 
module load fastqc/0.11.7
 
fastqc $input_file --outdir $output_folder
 
#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
