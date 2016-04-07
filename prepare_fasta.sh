#!/bin/sh


STEM=`echo $1 | cut -f1 -d'.'`

java -jar /research/picard-tools-1.119/CreateSequenceDictionary.jar "R="$1 "O="$STEM".dict"

samtools faidx $1

bwa index -a bwtsw $1
