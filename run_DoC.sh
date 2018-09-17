#!/bin/sh
#Quick script to run DoC on a bam
#reference_path input_bam

out=$(ls $2 |  cut -f1 -d'_')

samtools index -@ 24 $2

java -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $1 \
--omitDepthOutputAtEachBase \
--omitIntervalStatistics \
-I $2 \
-o "DoC_"$out \
-nt 24

for i in $(ls DoC* | grep -v "sample_summary"); do rm $i; done
