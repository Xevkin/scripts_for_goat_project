from __future__ import division

import sys

#modified to set het sites to missing if < 0.3 reads support het

#vcf and depth file
# depth file should have one samples per row, tab separated
# sample, mean depth, minimum depth

depth_file = sys.argv[2]

sample=[]

mean_depth=[]

min_depth=[]

#store minimum and mean coverages
with open(depth_file) as file:

	for line in file:

		splitline=line.rstrip("\n").split("\t")

		sample.append(splitline[0])

		mean_depth.append(float(splitline[1]))

		min_depth.append(int(splitline[2]))

#print the assigned coverages
#for i in range(0,len(sample)):

#	sys.stderr.write(sample[i] + " " + str(mean_depth[i]) + " " + str(min_depth[i]) + "\n")

with open(sys.argv[1]) as file:

	for line in file:

		if line.startswith("#"):

			print(line.rstrip("\n"))

			continue

		split_vcf_line = line.split("\t")

		new_line = split_vcf_line[0:9]

		n = 0

		for individual in split_vcf_line[9:]:

			DP = int(individual.split(":")[2].rstrip("\n"))

			SP = int(individual.split(":")[3].rstrip("\n"))

			if (DP < min_depth[n]) or (DP > int(round(mean_depth[n]*2))) or (SP > 13):

				#come back to this
				#print str(DP) + " " + str(min_depth[n]) + " " + str(int(round(mean_depth[n]*2)))

				blank_genotype="./."

				new_entry = blank_genotype + ":" + ":".join(individual.split(":")[1:])

				new_line.append(new_entry)

			elif (individual.split(":")[0].rstrip("\n") == "0/1"):

				ref = float(individual.split(":")[6].rstrip("\n").split(",")[0])

				alt = float(individual.split(":")[6].rstrip("\n").split(",")[1])

				total = ref + alt

				prop = float(float(alt)/float(total))

				if (prop < 0.3):

					blank_genotype="./."

	                                new_entry = blank_genotype + ":" + ":".join(individual.split(":")[1:])

        	                        new_line.append(new_entry)

				else:

					new_line.append(individual)

			else:

                        	new_line.append(individual)


			n += 1

		print "\t".join(new_line).rstrip("\n")
