import sys

count_list = []

#input is a vcf file
vcf_file = sys.argv[1]

with open(vcf_file) as file:

	for line in file:

		if line.startswith("#C"):

			continue

		spline = line.rstrip("\n").split("\t")

		if count_list == []:

			for individual in spline[9:]:

				count_list.append(0)

		#counter to keep track of individuals
		i = 0

		for individual in spline[9:]:

			if individual.startswith("0/1"):

				count_list[i] = count_list[i] + 1

			i = i + 1



for ind in count_list:

	print ind
