#!/bin/bash

## Owner: Marta Verdugo
## Created: 10/03/2015
## Last updated: 10/03/2015
# script for checking indexes for multiple samples
 
touch indexes_check.txt

for i in $(ls *.fastq.gz) 
do 
	echo $i  >> indexes_check.txt
	zcat $i | grep -P "^@MG" | cut -f 10 -d ":"| sort | uniq -c  >> indexes_check.txt 2> indexes_check.log; done 
 



