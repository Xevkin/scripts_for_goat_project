import sys

#regions: chr region_start region_end
regions = sys.argv[1]

#chr: chr length
chrom = sys.argv[2]

with open(chrom) as chromosomes:

		for line in chromosomes:

			start = 0

			split_chromosomes = line.rstrip("\n").split(" ")

			with open(regions) as region_file:

				for region in region_file:

					split_region = region.rstrip("\n").split(" ")

					if (split_region[0] == split_chromosomes[0]):

						print split_region[0] + " " + str(start) + " " + str(int(split_region[1]) - 1)

						start = int(split_region[2]) + 1

				print split_chromosomes[0] + " " + str(start) + " " + str(int(split_chromosomes[1]) - 1 )
