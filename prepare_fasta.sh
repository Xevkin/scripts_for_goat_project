#!/bin/sh


STEM=`echo $1 | cut -f1 -d'.'`

java -jar /home/kdaly/programs/picard/picard.jar CreateSequenceDictionary  "R="$1 "O="$STEM".dict"

samtools faidx $1

bwa index -a bwtsw $1
