import sys

#give the input file and the column of the ancestral individual

anc_col = int(sys.argv[2])

with open(sys.argv[1]) as file:

	for line in file:

		line = line.rstrip("\n")

		split_line = line.split("\t")

		if line.startswith("chr"):

			spl1 = [split_line[1]]

			new = spl1 + split_line[3:]

			print "\t".join(new)

		else:

			split_line = split_line[0:len(split_line)-1]

			ancestral = split_line[anc_col]

			sampled = split_line[3:]

			sampled_alleles = set(sampled)

			sampled_alleles.remove(ancestral)

			if "N" in sampled_alleles:

				sampled_alleles.remove("N")

			new_line = split_line[1] + "\t"

			for ind in sampled:

				sampled_alleles_list=list(sampled_alleles)

				if (ind == ancestral):

					new_line = new_line + "0\t"

				elif (ind == "N"):

					new_line = new_line + "NA\t"

				elif (ind == sampled_alleles_list[0]):

					new_line = new_line + "1\t"

				else:

					new_line = new_line + "2\t"


			print new_line
