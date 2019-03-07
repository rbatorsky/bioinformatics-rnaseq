#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --nodes=1			
#SBATCH --cpus-per-task=1	
#SBATCH --mem=20Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load subread/1.6.3

REF=/cluster/tufts/bio/data/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/sacCer3.gtf
DIR=$1

featureCounts \
-a $REF \
-o WT_1_featurecounts_gene_results.txt \
$DIR/*bam

