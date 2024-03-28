import sys

#SMUDGE factor if missing windows
SMUDGE = int(sys.argv[3])

#minimum no. varying sites in window - discard otherwise
MINSITES = int(sys.argv[4])

BLOCK = []

LAST_CHR = 1

LAST_CHR2 = 0

LAST_END = 1

COUNT = 0

NUM_BIN = 1

ROH_LIST = []

LAST_SAMPLE = ""

HET_LIMIT = float(sys.argv[2])

with open(sys.argv[1]) as file:

	for line in file:

		line = line.rstrip("\n")

		split_line = line.split(" ")

		if int(split_line[5]) < MINSITES:

			continue

		split_line.append("A_not_roh")

		CURRENT_SAMPLE = str(split_line[6])

		if CURRENT_SAMPLE != LAST_SAMPLE:

			for x in BLOCK:

				print " ".join(x)

			BLOCK = []

			SITE_BLOCK = []

		LAST_SAMPLE = CURRENT_SAMPLE

		CURRENT_CHR = int(split_line[0])

		HET_PROP = float(split_line[8])

		HET_SUM = 0

		if CURRENT_CHR > LAST_CHR:

			for x in BLOCK:

				print " ".join(x)

			BLOCK = []

			SITE_BLOCK = []

		LAST_CHR = CURRENT_CHR

		if (LAST_CHR2 == 0):

			LAST_CHR2 = CURRENT_CHR

			LAST_START = int(split_line[1])

		#within 10 windows i.e. 5000000
		if ((LAST_END + SMUDGE) <= (int(split_line[1])-1) ):

			for x in BLOCK:

				print " ".join(x)

			BLOCK = []

		LAST_END2 = LAST_END

		LAST_END = int(split_line[2])

		if len(BLOCK) < 10:

			BLOCK.append(split_line)

			continue

		else:

			print " ".join(BLOCK[0])

			BLOCK = BLOCK[1:11]

			BLOCK.append(split_line)

		for x in BLOCK:

			if float(x[8]) > HET_LIMIT:

				HET_SUM = 1

		if len(BLOCK) == 10:

			#if there are a 10 500K windows with lower than the het prop limit
			if (HET_SUM != 1) and (LAST_END > 0):

				if ((int(split_line[1]) - 1) <= (LAST_END2 + SMUDGE)) :

					for x in range(0,10):

						if BLOCK[x][-1] == "A_not_roh":

							BLOCK[x][-1] = "B_roh"

					CURRENT_START  = LAST_START

#					LAST_CHR2 = CURRENT_CHR

#					LAST_START = int(split_line[1])
#
#				print " ".join(split_line[0:2]) + " " + str(LAST_END) + " " + str(HET_SUM)

for x in BLOCK:

	print " ".join(x)
