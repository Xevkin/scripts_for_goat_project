#!/usr/bin/env python

import sys

#expects single input file
#input file: for each line: fastq.gz name i.e. Tube1, sample name,number of indexes, indexes, hiseq number (hiseqX), lane number (laneX), PCR number (PCRX), and platform
#comma seperated

for line in open(sys.argv[1]):

	entry = line.rstrip("\n").split(",")
	tube = entry[0]
	sample = entry[1]
	index_number = entry[2]
	last_index = 3 + int(index_number)
	indexes = entry[3:last_index]
	hiseq_number = entry[last_index]
	lane_number = entry[last_index + 1]
	pcr_number = entry[last_index + 2]
	platform = entry[last_index + 3]

	count = 1
	for i in indexes:

		print tube + "_" + str(count) + "_1.fastq.gz" + "\t" + r"@RG\tID:" + sample + "-index" + i + "-hiseq" + hiseq_number + "-lane" + lane_number + r"\tSM:" + sample + r"\tPL:" + platform + r"\tLB:" + sample + "-index" + i + "-PCR" + pcr_number + "\tlane" + lane_number + "\t" + sample 

		count += 1
