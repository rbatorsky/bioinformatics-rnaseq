#!/bin/bash

mkdir -p WT_1

while read $acc; 
do 
	echo $acc
	egrep $acc samples_at_ENA.txt| cut -f11 | xargs wget
done< <(cat ERP004763_sample_mapping.tsv | awk '$4 == 1' | awk '$3 == WT'  | awk '{print $1}')

mv *fastq.gz WT_1