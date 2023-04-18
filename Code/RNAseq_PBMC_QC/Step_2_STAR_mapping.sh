#!/bin/bash -l
 
#SBATCH -J STAR_map
#SBATCH --mem=350000
#SBATCH -n 8
#SBATCH --array=1-26
#SBATCH -p fn_short
#SBATCH --export=ALL
#SBATCH --time=2400
#SBATCH --output="STAR_map-%A_%a.out"
 
#date
d1=$(date +%s)

path='/gpfs/Data/Saracatinib_allsamples'

read1=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 1`
read2=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 2`
output_name=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 3`


echo $HOSTNAME
echo $read1
echo $read2
echo $output_name
 
module load star/2.6.1d
 
STAR --runMode alignReads --readFilesCommand zcat --genomeDir /gpfs/Data/genomes/hg38/STAR --readFilesIn $path/merged_data/$read1.fastq.gz $path/merged_data/$read2.fastq.gz --runThreadN 20 --outFileNamePrefix $path/STAR_mapped_data/$output_name. --outSAMprimaryFlag AllBestScore --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate  --sjdbOverhang 50 --quantMode GeneCounts

#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
