#!/usr/bin/python
#script takes the corrected (FP) ROH file and applies it

#we want to take the snp positions from the former, by chromosome, remove everything but those positions in the vcf

#must give chromosome number

import sys

def get_variants(cor_roh, chr_num):

	last_var = 0 
	
	
	chr_variants = []

	for chrs in range(0, int(chr_num.rstrip("\n"))):

		chr_variants.append([str(chrs + 1), []]) 

	with open(cor_roh) as file:
		
		n = 0

		for line in file:
	
			n += 1
							
			#go through each line in the roh file

			if not line.startswith("Chr"):

				split_line = line.strip("\n").split(" ")
					
				#as chromosome positions are in "CELA" format, need to cut off "CELA_"
				
				split_line[0] = split_line[0].split("_")[1]
					
				if (split_line[0] == "X"):
	
					split_line[0] = "30"
				
				#check if we've reached a new chromosome to reset the counter
									
				if (int(split_line[1]) <  int(last_var)):

					n = 1
					
				#for each chromosome
				
				for chr in chr_variants:
						
					if (str(split_line[0]) == chr[0]):
							
						if (n == 1):

							chr[1].append(split_line[1])
								
						chr[1].append(split_line[2])  
						
				last_var = split_line[2]			
	
	return chr_variants

def filter_vcf(vcf_file, variants):

	entries_to_return = []
		
	with open(vcf_file) as file:

		for line in file:

			if not line.startswith("#"):

				split_line = line.strip("\n").split("\t")

				chr = split_line[0].split("_")[1]
		
				if (chr == "X"):

					chr = "30"
				
				for variant_list in variants:
					
					if (variant_list[0] == chr):

						if (split_line[1] in variant_list[1]):

							entries_to_return.append(line.rstrip("\n"))			
							print line.rstrip("\n")
					

def main(cor_roh, vcf_file, chr_num):

	filter_vcf(vcf_file, get_variants(cor_roh, chr_num))



main(sys.argv[1], sys.argv[2], sys.argv[3])
