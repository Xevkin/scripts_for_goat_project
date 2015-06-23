#!/bin/sh

#input must be in the format <script> <vcf to convert to plink> <output_name>

vcftools --vcf $1 --plink --out $2

#now produce bed, fam and bim files

plink --file $2 --make-bed  --out $2
