import gzip

import sys

input_file = sys.argv[1]

with gzip.open(input_file) as file:

	for line in file:

		if any("N" in s for s in line.split()[3:7]):

			skip

		else:

			print line

