#!/bin/bash

module unload openmpi
module load fastqc/0.11.8
module load subread/1.6.3
module load anaconda
module load samtools/1.2
module load STAR/2.7.0a

mkdir -p course_data
cp /cluster/tufts/bio/tools/tutorials/bioinformatics-rnaseq-spring-20109/* course_data/
