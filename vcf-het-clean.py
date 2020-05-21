import sys
import gzip

with gzip.open(sys.argv[1], "rt") as file:

	for line in file:

		line = line.rstrip("\n")

		if line.startswith("#"):

			print line

			continue

		spline = line.split("\t")

		root = spline[0:9]

		for individual in spline[9:]:

			if not individual.startswith("0/1"):

				root.append(individual)

			else:

				allele_depth = individual.split(":")[6].split(",")

				if "0" in allele_depth:

					root.append("./." + individual.split(":")[1:])

				else:


					root.append(individual)


		print "\t".join(root)
