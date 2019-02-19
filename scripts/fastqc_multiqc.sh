#!/bin/bash

module load fastqc/0.11.8
module load anaconda
source activate /cluster/tufts/bio/tools/conda_envs/multiqc/1.7

mkdir -p fastqc
for fastq in WT_1/*fastq.gz;
do 
	fastqc $fastq --extract -o fastqc
done

multiqc fastqc --dirs -o multiqc

