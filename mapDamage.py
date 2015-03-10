#this script is similar to the miseq_handling script, and should be merged at a later date
#similar set-up to that script: input file list, then pipe to mapDamage
from subprocess import call
import sys

#we will use this later to check if the input files actually exist
import os.path

#reference genomes
goat_ref = "~/goat/data/reference_genomes/goat_CHIR1_0/goat_CHIR1_0.fasta"

#sys.argv[1] refers to the list input
#create list that will carry any lines in the file which fail at any stage
failures = []
f = open(sys.argv[1])
#cycle through each line in the input file
for lines in f:
	
	current_file = lines.strip()
	
	stem = current_file.split(".")[0]
	print stem
	call("gunzip " + current_file,shell=True)

	call("mapDamage --merge-reference-sequences -i " + stem + ".bam -d ~/goat/results/2015_02_06/mapDamage/" + stem + " -r " + goat_ref + " -t " + stem,shell=True)

	call("gzip " + stem + ".fa",shell=True)	
