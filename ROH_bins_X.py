import sys

#input 1 is a het site file, tab separated: chrom, position, 0 for HOM, 1 for HET
#input 2 is length of window

#chromosomes lengths for ARS1 are defined
chromosomes = []

with open("/home/kdaly/xchromosomes.length") as file:

	for line in file:

		split_line = line.rstrip("\n").split("\t")

		chromosomes.append(split_line)

BIN_LENGTH = int(sys.argv[2]) - 1

BIN = 0

for CHR in chromosomes:

	CURRENT_CHR = int(CHR[0])

	CHR_END = int(CHR[1])

	BIN_START = 1

	BIN_END = BIN_START + BIN_LENGTH

	HET_COUNT = 0

	VAR_COUNT = 0

	SWITCH = 0

	with open(sys.argv[1]) as file:

		for line in file:

			line = line.rstrip("\n")

			split_line = line.split("\t")

			SAMPLE_CHR = int(split_line[0].replace("NW_017189516.1","1").replace("NW_017189517.1","2"))

			SAMPLE_POS = int(split_line[1])

			SAMPLE_HET = int(split_line[2])

			#if we have just come to the end of the chromosome update the chromosome and print the remaining windows
			if (SAMPLE_CHR == CURRENT_CHR + 1) and (SWITCH != 1):

				while BIN_END < CHR_END:

					if (HET_COUNT != 0):

						print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t" + str(HET_COUNT) + "\t1\t" + str(VAR_COUNT)

						HET_COUNT = 0

						VAR_COUNT = 0

					else:

						print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t0\t0\t" + str(VAR_COUNT)

						VAR_COUNT = 0

					BIN_START = BIN_START + BIN_LENGTH + 1

	                                BIN_END = BIN_START + BIN_LENGTH

				BIN_END = CHR_END

				if HET_COUNT != 0:

					BIN = "1"

				else:

					BIN = "0"

				print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t" + str(HET_COUNT) + "\t" + BIN + "\t" + str(VAR_COUNT)

				VAR_COUNT = 0

				SWITCH = 1

				continue

			#if the variant chromosome does not match the current chromosome, skip the variant
			if not SAMPLE_CHR == CURRENT_CHR:

				continue

			#keep printing empty windows if are not up to the variant
			while SAMPLE_POS > BIN_END:

				#if the last window had a het, print the window and reset the het
				if (HET_COUNT != 0):

					print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t" + str(HET_COUNT) + "\t1\t" + str(VAR_COUNT)

					HET_COUNT = 0

					VAR_COUNT = 0

				#otherwise, print empty windows
				else:

					print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t0\t0\t" + str(VAR_COUNT)

					VAR_COUNT = 0

				#update the bin
				BIN_START = BIN_START + BIN_LENGTH + 1

				BIN_END = BIN_START + BIN_LENGTH

			#if we are up to the variant, update the window as having an additional het and move onto the next variant
			if  (SAMPLE_POS <= BIN_END ) and (SAMPLE_POS >= BIN_START):

				VAR_COUNT = VAR_COUNT + 1

				if (SAMPLE_HET == 1):

					HET_COUNT = HET_COUNT + 1

					continue

#for the very last window
while BIN_END < CHR_END:

	if (HET_COUNT != 0):

		print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t" + str(HET_COUNT) + "\t1\t" + str(VAR_COUNT)

		HET_COUNT = 0

		VAR_COUNT =0

	else:

		print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t0\t0\t" + str(VAR_COUNT)

		VAR_COUNT = 0

	BIN_START = BIN_START + BIN_LENGTH + 1

	BIN_END = BIN_START + BIN_LENGTH

BIN_END = CHR_END

if HET_COUNT != 0:

	BIN = "1"

else:

	BIN = "0"

print str(CURRENT_CHR) + "\t" +  str(BIN_START) + "\t" + str(BIN_END) + "\t" + str(HET_COUNT) + "\t" + BIN + "\t" + str(VAR_COUNT)
