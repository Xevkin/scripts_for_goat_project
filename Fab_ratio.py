#!/usr/bin/env python

#Script will calculated F(A|B) ratio
#ratio is the proportion of times that for a randomly chosen chromosome in ancient A, it
#shares a het/derived allele that is present in B
#Supposedly, as we are randomly selecting a chromosome in A we are not affected by population demographies of the B ancestral pop
#We are for pop A, however


import sys
import random

A_vcf = sys.argv[1]

B_vcf = sys.argv[2]



with open(B_vcf) as line:

	if not line.startswith("#"):

		vcf_entries = line.split("\t")
		
		print vcf_entries[len(vcf_entries)-1].split(":")[0]
