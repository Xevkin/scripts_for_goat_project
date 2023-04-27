#!/usr/bin/env python

#remove N bps from both ends

import sys

remove = int(sys.argv[2])

dump = ""

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith(">"):

			print ">" + sys.argv[1].rstrip("\n").split("_")[0]

		else:

			dump = dump + line.rstrip("\n")

	dump = dump[remove:-remove]

	while (len(dump) > 70):

		print dump[0:70]

		dump = dump[70:]

	print dump
