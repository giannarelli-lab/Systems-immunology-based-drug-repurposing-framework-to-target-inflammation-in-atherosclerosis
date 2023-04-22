#!/bin/bash -l

#SBATCH -J fastqc
#SBATCH --mem=350000
#SBATCH -n 8
#SBATCH --array=369-392
#SBATCH -p fn_medium
#SBATCH --export=ALL
#SBATCH --time=2400
#SBATCH --output="fastqc1_raw_data-%A_%a.out"

#date
d1=$(date +%s)

module load fastqc/0.11.7


##### run fastqc in raw reads ####
path='/gpfs/Data/Saracatinib_allsamples'
path1='/gpfs/Data/Saracatinib_allsamples/fastq'

sample=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/raw_data_sample_names.txt`
input_file=$path1/$sample.fastq.gz
output_folder=$path/fastqc_qc_data/

echo $HOSTNAME
echo $sample
 
module load fastqc/0.11.7
 
fastqc $input_file --outdir $output_folder
 
#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
