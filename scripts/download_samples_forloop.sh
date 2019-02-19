#!/bin/bash

mkdir -p WT_1

for acc in ERR458493 ERR458494 ERR458495 ERR458496 ERR458497 ERR458498;
do
	echo $acc
	egrep ${acc} samples_at_ENA.txt | cut -f11 | xargs wget
done

mv *fastq.gz WT_1