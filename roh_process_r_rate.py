import sys

#SMUDGE factor if missing windows
SMUDGE = int(sys.argv[3])

BLOCK = []

LAST_CHR = 1

LAST_CHR2 = 0

LAST_END = 1

COUNT = 0

NUM_BIN = 1

ROH_LIST = []

LAST_SAMPLE = ""

RATE_LIMIT = float(sys.argv[2])

RUN_OF_WINDOWS = int(sys.argv[4])

TOLERANCE = int(sys.argv[5])

with open(sys.argv[1]) as file:

	for line in file:

		line = line.rstrip("\n")

		split_line = line.split(" ")

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

		RATE_BREAK = 0

		if CURRENT_CHR > LAST_CHR:

			for x in BLOCK:

				print " ".join(x)

			BLOCK = []

			SITE_BLOCK = []

		LAST_CHR = CURRENT_CHR

		if (LAST_CHR2 == 0):

			LAST_CHR2 = CURRENT_CHR

			LAST_START = int(split_line[1])

		#within the run of windows which we define in RUN_OF_WINDOWS
		if ((LAST_END + SMUDGE) <= (int(split_line[1])-1) ):

			for x in BLOCK:

				print " ".join(x)

			BLOCK = []

		LAST_END2 = LAST_END

		LAST_END = int(split_line[2])

		#define the number of ROH windows to combine
		if len(BLOCK) < RUN_OF_WINDOWS:

			BLOCK.append(split_line)

			continue

		else:

			print " ".join(BLOCK[0])

			BLOCK = BLOCK[1:(RUN_OF_WINDOWS + 1)]

			BLOCK.append(split_line)

		for x in BLOCK:

			#in our current chunk of X windows, is there any window with too high a rate?
			#if so, set the RATE_BREAK to 1
			if float(x[8]) > RATE_LIMIT:

				RATE_BREAK = RATE_BREAK + 1

		if len(BLOCK) == RUN_OF_WINDOWS:

			#if there are a X 500K windows with lower than the het prop limit
			#^if we haven't broke the rate limit
			if (RATE_BREAK <= TOLERANCE) and (LAST_END > 0):

				if ((int(split_line[1]) - 1) <= (LAST_END2 + SMUDGE)) :

					for x in range(0,RUN_OF_WINDOWS):

						if BLOCK[x][-1] == "A_not_roh":

							BLOCK[x][-1] = "B_roh"

					CURRENT_START  = LAST_START

for x in BLOCK:

	print " ".join(x)
