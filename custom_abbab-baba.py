import gzip

import sys

input_file = sys.argv[1]

with gzip.open(input_file) as file:

	for line in file:

		#remove Ns
		if any("N" in s for s in line.split()[3:7]):

			continue

		else:

			#remove transitions
			if (line.split()[2:3] == "C" and line.split()[3:4] == "T") or (line.split()[2:3] == "C" and line.split()[4:5] == "T") or (line.split()[2:3] == "G" and line.split()[3:4] == "A") or (line.split()[2:3] == "G" and line.split()[4:5] == "A"):

				continue

		if ((line.split()[2:3] ==  line.split()[3:4]) and  (line.split()[4:5] ==  line.split()[5:6]) and (line.split()[2:3] ==  line.split()[4:5])) or ((line.split()[2:3] ==  line.split()[4:5]) and  (line.split()[3:4] ==  line.split()[5:6]) and (line.split()[2:3] ==  line.split()[3:4])):

				print line
