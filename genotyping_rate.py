from __future__ import division

#takes ped input file and <>, gives both per sample genotyping rate and per breed

import sys

filter = sys.argv[2].rstrip("\n")

ped_file = sys.argv[1]

with open(ped_file) as f:

	for line in f:

		split_line = line.split(" ")

		#first six columns are not genotypes

		genotypes = split_line[7:]

		site_num = len(genotypes)

		sites_called = genotypes.count("A") + genotypes.count("G") + genotypes.count("C") + genotypes.count("T")

		genotyping_rate = sites_called / site_num
		
		#print genotyping_rate

		if (genotyping_rate < float(filter)):

			print split_line[1]
