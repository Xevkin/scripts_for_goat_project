#!/bin/bash

# input 1 is the fastq file to filter

# input 2 is P7 to match

INPUT=$1

OUT=`echo $INPUT | rev | cut -f1 -d'/' | rev | cut -f1 -d'.'`

P7=$2

INDEX=`echo ${2}+`

FINAL_OUT=`echo ${OUT}-P7-exact.fastq.gz | sed -e "s/_R1-P7-exact/-P7-exact_R1/g" | sed -e "s/_R2-P7-exact/-P7-exact_R2/g"`

zcat ${INPUT} | grep -A 3 $INDEX | grep -v "^--$" | gzip -c - > ${FINAL_OUT}
