import sys

min=int(sys.argv[2])

max=int(sys.argv[3])

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith("#"):

			continue

		split_vcf_line = line.split("\t")

		new_line = split_vcf_line[0:9]

		for individual in split_vcf_line[9:]:

			DP = int(individual.split(":")[2].rstrip("\n"))

			SP = int(individual.split(":")[3].rstrip("\n"))

			if (DP < min) or (DP > max) or (SP > 13):

				blank_genotype="./."

				new_entry = blank_genotype + ":" + ":".join(individual.split(":")[1:])

				new_line.append(new_entry)

			else:

				new_line.append(individual)


		print "\t".join(new_line).rstrip("\n")


