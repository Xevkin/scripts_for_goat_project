import sys

#converts to bed co-ordinates!! 0 based

with open(sys.argv[1]) as file:

	for line in file:

		split_line = line.split("\t")

		new_start = int(split_line[1])-50001

		if (new_start <= 0):

			new_start = "0"

		else:

			new_start = str(new_start)

		new_end = str(int(split_line[2])+49999)

		print "\t".join([split_line[0],new_start,new_end])
