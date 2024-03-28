#!/usr/bin/env python2

#python2 decircularize.py input.fa N sample_name

#remove N bps from both ends

import sys

remove = int(sys.argv[2])

dump = ""

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith(">"):

			print ">" + sys.argv[3].rstrip("\n")

		else:

			dump = dump + line.rstrip("\n")

	dump = dump[remove:-remove]

	while (len(dump) > 70):

		print dump[0:70]

		dump = dump[70:]

	print dump
