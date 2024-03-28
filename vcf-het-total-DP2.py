import sys

count_list = []

total_list = []

#input is a vcf file
vcf_file = sys.argv[1]

with open(vcf_file) as file:

	for line in file:

		if line.startswith("#"):

			continue

		spline = line.rstrip("\n").split("\t")

		if count_list == []:

			for individual in spline[9:]:

				count_list.append(0)

				total_list.append(0)

		#counter to keep track of individuals
		i = 0

		for individual in spline[9:]:

			#for each individual, if the individual is covered, add to the total site count
			#note that I am skipping over HET sites with 1 read
			if not individual.startswith("./."):

				total_list[i] = total_list[i] + 1

			if individual.startswith("0/1"):

				AD = individual.split(":")[6].split(",")

				#remove any het sites in which there are not at least one high quality base observed for each allele
				if AD.count("0") == 0:

					count_list[i] = count_list[i] + 1

				else:

					total_list[i] = total_list[i] - 1


			i = i + 1



print count_list

print total_list
