import sys

#print out the bins that do not have sufficient sites
BLOCK_SIZE = 500000

MINIMUM = 50

chromosomes = []

with open("/home/kdaly/chromosomes.length") as file:

	for line in file:

		split_line = line.rstrip("\n").split("\t")

		chromosomes.append(split_line)

for CHROMOSOME in chromosomes:

	CURRENT_COUNT = 0

	SWITCH = 0

	with open(sys.argv[1]) as file:

		for line in file:

			if SWITCH == 1:

				break

			line = line.rstrip("\n")

			spline = line.split("\t")

			CHR = int(spline[0])

			POS = int(spline[1])

			HET = int(spline[2])

			#if we moved past the chromosome,print out the remaining blocks

			if (CHR == int(CHROMOSOME[0]) + 1) and SWITCH == 0:

				SWITCH = 1

				while BLOCK[1] < int(CHROMOSOME[1]):

					print str(CHR - 1) + " " + str(BLOCK[0]) + " " + str(BLOCK[1]) + " " + str(CURRENT_COUNT)

					CURRENT_COUNT = 0

					BLOCK[0] = BLOCK[1]

					BLOCK[1] = BLOCK[1] + BLOCK_SIZE

				BLOCK[1] = int(CHROMOSOME[1])

				print str(CHR - 1) + " " + str(BLOCK[0]) + " " + str(BLOCK[1]) + " " + str(CURRENT_COUNT)

				BLOCK = [0,BLOCK_SIZE]

			# if we're not on the right csome, just move on
			if CHR != int(CHROMOSOME[0]):

				continue

			#if we've moved past the block, start moving the block forward
			while POS > BLOCK[1]:

				print str(CHR) + " " + str(BLOCK[0]) + " " + str(BLOCK[1]) + " " + str(CURRENT_COUNT)

				CURRENT_COUNT = 0

				BLOCK[0] = BLOCK[1]

				BLOCK[1] = BLOCK[1] + BLOCK_SIZE

			#if we're in the correct block, print the variant count
			if POS >= BLOCK[0] and POS < BLOCK[1]:

				CURRENT_COUNT = CURRENT_COUNT + 1
