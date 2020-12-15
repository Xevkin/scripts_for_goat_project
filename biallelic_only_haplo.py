import sys

with open(sys.argv[1]) as file:

	for line in file:

		line = line.rstrip("\n")

		if line.startswith("chr"):

			print line

		else:

			spline = line.split("\t")

			count = 0

			merged_spline = "".join(spline[3:len(spline)])

			for base in ["A", "G", "T", "C"]:

				if base in merged_spline:

					count += 1

			if count == 2:

				print line
