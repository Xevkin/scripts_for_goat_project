#!/bin/sh
#Quick script to run DoC on an indexed bam
#reference_path input_index_bam

out=$(ls $2 |  cut -f1 -d'_')

java -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $1 \
--omitDepthOutputAtEachBase \
-I $2 \
-o "DoC_"$out

