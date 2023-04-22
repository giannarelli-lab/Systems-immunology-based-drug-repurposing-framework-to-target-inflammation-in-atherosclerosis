#!/bin/bash -l
 
#SBATCH -J FeatCount
#SBATCH --mem=300000
#SBATCH -n 12
#SBATCH -p fn_medium
#SBATCH --export=ALL
#SBATCH --time=24:00:00
#SBATCH --output="FeatCounts-%A.out"

#date
d1=$(date +%s)
 
echo $HOSTNAME
 
myoutfile='/gpfs/Data/Saracatinib_allsamples/FeatureCounts'
mybams=`cat /gpfs/Data/Saracatinib_allsamples/key_to_FeatureCounts.txt`

module load subread/1.6.3
 
featureCounts -T 12 -s 2 -p -g gene_id --extraAttributes gene_type,gene_name -a /gpfs/Data/genomes/Homo_sapiens.GRCh38.107.chr.gtf -o $myoutfile/featurecounts.txt $mybams
 
d2=$(date +%s)
sec=$(( ( $d2 - $d1)))
hour=$(echo - | awk '{print '$sec'/3600}')
echo Runtime: $hour hours \($sec\s\)
