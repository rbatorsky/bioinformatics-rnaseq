#!/bin/bash
#SBATCH --job-name=job		    
#SBATCH --nodes=1			    
#SBATCH --cpus-per-task=12		
#SBATCH --mem=20Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# load the module 
module load STAR/2.5.2b

# create a directory to store the index in
REF_DIR=/cluster/tufts/bio/data/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3
STAR_DIR=${REF_DIR}/Sequence/STAR

# Run STAR in "genomeGenerate" mode
STAR --runMode genomeGenerate \
--genomeDir $STAR_DIR \
--genomeFastaFiles ${REF_DIR}/Sequence/WholeGenomeFasta/genome.fa \
--sjdbGTFfile ${REF_DIR}/Annotation/Genes/sacCer3.gtf \
--sjdbOverhang 49 \
--runThreadN 12 
