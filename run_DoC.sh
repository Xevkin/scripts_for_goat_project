#!/bin/sh
#Quick script to run DoC on an indexed bam
#reference_path input_index_bam

out=$(ls $2 |  cut -f1 -d'_')

java -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $1 \
--omitDepthOutputAtEachBase \
--omitIntervalStatistics \
-I $2 \
-o "DoC_"$out \
-nt 24
