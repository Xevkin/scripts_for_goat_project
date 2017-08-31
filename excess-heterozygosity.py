import sys

chr = sys.argv[2]

region_start = int(sys.argv[3])

region_end = int(sys.argv[4].rstrip("\n"))

#use base 1 columns and python will put in base 0!!

to_extract = []

for i in sys.argv[5:]:

	to_extract.append(int(i)-1)

bin_start = region_start

bin_end = region_start + 99

#danger count tracks how many problematic het-rich sites are in a bin
danger_count = 0

danger="no"

with open(sys.argv[1]) as file:

	for line in file:

		het_count = 0

		split_line = line.split("\t")

		if line.startswith("#"):

			continue

		if line.startswith(chr):

			#skip if line is smaller than bin start
			if (split_line[1] < bin_start):

				continue

			#end the loop if we are past the region of interest

			if (int(split_line[1]) > region_end):

				#reset as we have to check at the end of the script

				het_count = 0

				break

			#if the bin is finished, we need to update

			while (int(split_line[1]) > bin_end & bin_end <= region_end):

				print "updating"

				print str(danger_count)

				#we need to check how many danger sites were found

				if (danger_count >= int(round((bin_end - bin_start) / 10)) ):

					danger="yes"

				#reset danger counter for next bin
				danger_count = 0

				#move the bin forward
				bin_start += 100

				bin_end += 100

				if (bin_end > region_end):

					bin_end = region_end

			#check if bin 
			if ((int(split_line[1]) >= bin_start) & (int(split_line[1]) <= bin_end)):

				for ind in [split_line[i] for i in to_extract]:

					if ind.startswith("0/1"):

						het_count += 1

						print "adding a het:"  + str(het_count)

			#check hets at the end of line
			if (het_count >= (round(len(to_extract) / 2)) ):

				danger_count += 1

if (danger_count >= int(round((bin_end - bin_start) / 10)) ):

	danger="yes"

print str(danger_count)

if (danger=="yes"):

	print "Too SNPpy region/s detected"
