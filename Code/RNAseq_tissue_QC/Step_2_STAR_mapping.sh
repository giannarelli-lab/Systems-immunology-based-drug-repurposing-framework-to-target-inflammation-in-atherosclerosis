#!/bin/bash -l
 
#SBATCH -J STAR_map
#SBATCH --mem=250000
#SBATCH -n 8
#SBATCH --array=1-9
#SBATCH -p cpu_short
#SBATCH --export=ALL
#SBATCH --time=600
#SBATCH --output="STAR_map-%A_%a.out"
 
#date
d1=$(date +%s)

path='/gpfs/Data/RNAseq/06102022'

read1=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 1`
read2=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 2`
output_name=`sed -n "${SLURM_ARRAY_TASK_ID}p" $path/key_to_STAR.txt | cut -d' ' -f 3`

echo $HOSTNAME
echo $read1
echo $read2
echo $output_name
 
module load star/2.7.7a
 
STAR --genomeDir /gpfs/Data/genomes/hg38/STAR --readFilesCommand zcat --readFilesIn $path/merged_data/$read1.fastq.gz $path/merged_data/$read2.fastq.gz --outSAMprimaryFlag AllBestScore --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $path/STAR_mapped_data/$output_name. --quantMode GeneCounts
STAR   --runMode genomeGenerate   --runThreadN 20   --genomeDir /gpfs/Data/genomes/hg38/STAR_51_2.7.7a   --genomeFastaFiles /gpfs/Data/genomes/hg38/GRCh38.primary_assembly.genome.fa      --sjdbGTFfile /gpfs/Data/genomes/hg38/gencode.v36.annotation.gtf   --sjdbOverhang 50


#date
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
