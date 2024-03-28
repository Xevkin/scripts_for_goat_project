import random
import sys

with open(sys.argv[1]) as file:

	for line in file:

		line = line.rstrip("\n")

		split_line = line.split("\t")

		new_line = split_line[0:6]

		site_count =  len(split_line[6:]) / 2

		for site in range(6, 6 + (site_count*2), 2):

			selected = split_line[site+random.randint(0,1)]

			new_line.append(selected)

			new_line.append(selected)

		print "\t".join(new_line)
