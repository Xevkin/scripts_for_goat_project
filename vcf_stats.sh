#!/bin/sh

#filters heterozygous positions from the given vcf file]
#Very inefficient right now - creates additional files which seems unnecessary
#first section of code returns percentage of heterozgous positions

#file to process
echo "Input file to get stats for:"
read file
#count the number of header lines
awk '{if ($1 ~ /^#/ ) print $0}' $file > tmp

echo "Number of header lines is:"
HEADER=$(wc -l tmp| cut -f1 -d' ')
echo $HEADER

#filter for heterozygous positions
awk '{if ($1 !~ /^#/ && $10 ~ /0\/1/) print $0}' $file > tmp

#get a count of the number of heterozygous positions
echo "Number of heterozygous positions:"
HETS=$(wc -l tmp | cut -f1 -d' ')
echo $HETS

#count the number of lines in the file
echo "Number of lines in the vcf file is:"
VCF=$(wc -l $file | cut -f1 -d' ')

echo $VCF

echo "Total number of variants is:"
VARIANTS=$(($VCF-$HEADER))
echo $VARIANTS

echo "Proportion of heterozygous sites:"
echo "$HETS $VARIANTS" | awk '{printf "%.3f \n", 100*$1/$2}'

rm tmp
