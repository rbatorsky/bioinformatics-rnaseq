#!/bin/bash

module load fastqc/0.11.8
module load anaconda
source activate /cluster/tufts/bio/tools/conda_envs/rnaseq_course

mkdir -p fastqc
for fastq in WT_1/*fastq.gz;
do 
	fastqc $fastq --extract -o fastqc
done

multiqc fastqc --dirs -o multiqc

