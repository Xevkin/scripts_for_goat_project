#!/bin/sh
#Quick script to run DoC on a bam
#reference_path input_bam

out=$(ls $2 |  cut -f1 -d'_')

root=`echo $i | cut -f1 -d'.'`


if [ -f ${root}.bam.bai ]; then

	echo "Indexing " $2;

	samtools index -@ 24 $2

fi


/raid_md0/Software/gatk \
 DepthOfCoverage \
-R $1 \
--omit-depth-output-at-each-base \
--omit-interval-statistics \
-I $2 \
-O "DoC-autosomes_"$out \
-L /home/kdaly/autosomes_sheep.list

for i in $(ls DoC* | grep -v "sample_summary"); do rm $i; done
