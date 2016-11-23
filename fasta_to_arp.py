import sys

sequence = ""
sample = ""

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith(">"):

			if (len(sequence) > 0):

				print sample.replace(">","") + "\t1\t" + sequence

			sample = line.strip()

			sequence = ""

		else:

			sequence = sequence + line.strip()
