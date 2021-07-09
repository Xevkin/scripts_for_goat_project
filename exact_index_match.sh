#!/bin/bash

# input 1 is the fastq file to filter

# input 2 is P7+P5 to match

INPUT=$1

OUT=`echo $INPUT | rev | cut -f1 -d'/' | rev | cut -f1 -d'.'`

INDEX=$2

zcat ${INPUT} | grep -A 3 $INDEX | grep -v "^--$" | gzip -c - > ${OUT}_exact.fastq.gz
