#!/bin/bash

# input 1 is the fastq file to filter

# input 2 is P7 to match

# input 3 is P5

INPUT=$1

OUT=`echo $INPUT | rev | cut -f1 -d'/' | rev | cut -f1 -d'.'`

P7=$2

P5=$3

INDEX=`echo ${2}+${3}`

FINAL_OUT=`echo ${OUT}_exact.fastq.gz | sed -e "s/_R1_exact/-exact_R1/g" | sed -e "s/_R2_exact/-exact_R2/g"`

zcat ${INPUT} | grep -A 3 $INDEX | grep -v "^--$" | gzip -c - > ${FINAL_OUT}
