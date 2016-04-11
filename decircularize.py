#!/bin/sh

#remove 15 bps from both ends

import sys

dump = ''

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith(">"):

			print line.strip("\n")

		else:
			
			dump = dump + line.rstrip("\n")

	dump = dump[15:-15]
	
	
	while (len(dump) > 70):

		print dump[0:70]

		dump = dump[70:]		

	print dump
