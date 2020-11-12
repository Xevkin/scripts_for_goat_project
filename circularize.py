#!/usr/bin/env python

#adds 15 bps from the start to the end and vice versa

import sys

dump = ''

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith(">"):

			print line.strip("\n")

		else:
			
			dump = dump + line.rstrip()

	head = dump[0:15]

	tail = dump[-15:]

	dump = tail + dump + head
	
	
	while (len(dump) > 70):

		print dump[0:70]

		dump = dump[70:]		

	print dump
