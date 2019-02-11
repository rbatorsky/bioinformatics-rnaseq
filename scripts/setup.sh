#!/bin/bash
module unload openmpi
module load fastqc/0.11.8
module load subread/1.6.3
module load anaconda
module load samtools/1.2
module load STAR/2.7.0a
source activate /cluster/kappa/90-days-archive/bio/tools/conda_envs/rnaseq_course
mkdir course_data
cp /cluster/kappa/90-days-archive/bio/tools/tutorials/bioinformatics-rnaseq-spring-20109/* course_data/
