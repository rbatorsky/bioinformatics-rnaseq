#!/bin/bash
module load anaconda
source activate /cluster/tufts/bio/tools/conda_envs/rseqc/3.0.0

mkdir -p bamqc

BAM=$1
BAM_BASE=$( basename $BAM)
REF_DIR=/cluster/tufts/bio/data/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/

bam_stat.py -i $BAM > bamqc/${BAM_BASE%.bam}.bam_stat.txt
read_distribution.py -r ${REF_DIR}/sacCer3.sgdGene.wholegene.bed  -i $BAM > bamqc/${BAM_BASE%.bam}.read_dist.txt

## need to unload and load again because of samtools version conflict
source deactivate
module load samtools/1.2

samtools flagstat $BAM > bamqc/${BAM_BASE}_flagstat_WT_1.txt

source activate /cluster/tufts/bio/tools/conda_envs/multiqc/1.7

multiqc bamqc -o multiqc
