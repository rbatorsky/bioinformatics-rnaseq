#!/bin/bash

## example bash script to download 5 replicates of WT and SNF2 in a foor loop

for rep in 1 2 3 4 5;
do
	for type in WT SNF2;
	do
		mkdir -p ${type}_${rep}
		cd ${type}_${rep}

		while read p; 
			do 
   		 		egrep $p ../samples_at_ENA.txt| cut -f11 | xargs sbatch ../sbatch_wget.sh 
				done< <(cat ../ERP004763_sample_mapping.tsv | awk -v rep="$rep" '$4 == rep' | awk -v type="$type" '$3 == type'  | awk '{print $1}')
    			
		cd ../
	done
done
