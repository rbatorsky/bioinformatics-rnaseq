#!/bin/bash
#SBATCH --job-name=star_job           
#SBATCH --nodes=1                       
#SBATCH --cpus-per-task=12               
#SBATCH --mem=64Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

## Use STAR aligner to align all fastq files in a directory
## This step must be done for each sample
module load STAR/2.7.0a
mkdir -p STAR

## File is run as sbatch sbatch_align_star.sh <folder>
## Where <folder> contains fastq.gz files for 1 sample
SAMPLE=$1

## Obtain a comma separated list of files
#FILES=`ls -m WT_1/*fastq.gz | sed 's/ //g'`
#FILES=`echo $FILES `
FILES=`ls -m ${SAMPLE}/*fastq.gz | tr -d ' ' | tr -d '\n'`

## Name the output file, for example if the folder is /cluster/tufts/sample_1
## The output will have prefix sample_1_
OUT=$(basename $SAMPLE)

echo "Starting to align: $FILES"
echo "Output file will have prefix: $OUT"

REF_DIR=/cluster/tufts/bio/data/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/STAR

# execute STAR in the runMode "alignReads"
STAR --genomeDir ${REF_DIR} \
--readFilesIn $FILES \
--readFilesCommand zcat \
--outFileNamePrefix STAR/${OUT}_ \
--outFilterMultimapNmax 1 \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--runThreadN 12 \
--alignIntronMin 1 \
--alignIntronMax 2500

# generate the bam index
module load samtools/1.2
samtools index STAR/${OUT}_Aligned.sortedByCoord.out.bam
