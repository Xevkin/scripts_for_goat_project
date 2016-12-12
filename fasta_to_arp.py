#!/usr/bin/env python

import sys

sequence = ""
sample = ""

print "[Profile]\n\nTitle=\"\"\nNbSamples=1\nGenotypicData=0\nDataType=DNA\nLocusSeparator=NONE\nMissingData='N'\n"

print "[Data]\n\n[[Samples]]\n\nSampleName=\"\"\nSampleSize=\nSampleData= {\n\n"

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith(">"):

			if (len(sequence) > 0):

				print sample.split(" ")[0].replace(">","") + "\t1\t" + sequence

			sample = line.strip()

			sequence = ""

		else:

			sequence = sequence + line.strip()


#print the final sequences

print sample.split(" ")[0].replace(">","") + "\t1\t" + sequence


print "\n}\n\n[[Structure]]\n\nStructureName=\"\"\nNbGroups=1\nGroup={\n\"\"\n}"
