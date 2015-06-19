#!/usr/bin/python

#script takes input vcf, parses for transversions (A<->T/C;G<->T/C), and prints those lines 

import sys

A_G = ["A", "G"]

T_C = ["T", "C"]

input = sys.argv[1]

with open(input) as file:

	for line in file:

	#print the line if it begins with #

		if line.startswith(r"#"):

			print line

		else:

			split_line = line.split("\t") 

			if (((split_line[3] in A_G) and (split_line[4] in T_C)) or ((split_line[4] in A_G) and (split_line[3] in T_C))):

				print line.rstrip("\n") 
	
