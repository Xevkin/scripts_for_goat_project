#!/usr/bin/python
from __future__ import division

import sys

total_bases = 0

total_N = 0

with open(sys.argv[1]) as file:

	for line in file:

		if not line.startswith(">"):

			line = line.rstrip("\n").upper()

			total_bases = total_bases + len(line)

			total_N = total_N + line.count("N")

proportion = total_N / total_bases


print "Total N: " + str(total_N) 
print "Total bases: " + str(total_bases)
print "Proportion of N in genome: " + str(proportion) 
