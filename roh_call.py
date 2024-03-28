import sys


BLOCK = []

LAST_CHR = 1

LAST_CHR2 = 0

LAST_END = 1

COUNT = 0

NUM_BIN = 1

ROH_LIST = []

HET_LIMIT = int(sys.argv[2])

with open(sys.argv[1]) as file:

	for line in file:

		line = line.rstrip("\n")

		split_line = line.split("\t")

		CURRENT_CHR = int(split_line[0])

		HET_COUNT = int(split_line[3])

		if CURRENT_CHR > LAST_CHR:

			BLOCK = []

		LAST_CHR = CURRENT_CHR

		if (LAST_CHR2 == 0):

			LAST_CHR2 = CURRENT_CHR

			LAST_START = int(split_line[1])

		if len(BLOCK) < 10:

			BLOCK.append(HET_COUNT)

			continue

		else:

			BLOCK = BLOCK[1:11]

			BLOCK.append(HET_COUNT)

		HET_SUM = sum(BLOCK)

		if (HET_SUM < HET_LIMIT) and (LAST_END > 0):

			if ((int(split_line[1]) - 1) == LAST_END) :

				NUM_BIN = NUM_BIN + 1

				CURRENT_START  = LAST_START

			else:

				ROH_LIST.append([LAST_CHR2, LAST_START, LAST_END , NUM_BIN])

				NUM_BIN = 1

				LAST_CHR2 = CURRENT_CHR

				LAST_START = int(split_line[1])

			#print " ".join(split_line[0:2]) + " " + str(LAST_END) + " " + str(HET_SUM)

			LAST_END = int(split_line[2])

for ROH in range(1, len(ROH_LIST)):

	print " ".join(str(i) for i in ROH_LIST[ROH])
