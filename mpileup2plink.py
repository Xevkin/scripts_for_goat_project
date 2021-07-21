import sys

import re

from random import randint

#list of samples
SAMPLE_FILE = sys.argv[2]

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		if LINE.startswith("##"):

			continue

		SPLINE = LINE.rstrip("\n").split("\t")

		REF = SPLINE[2].upper()

		TO_PRINT = SPLINE[0] +  " " + SPLINE[0] + "_" + SPLINE[1] + " 0 " + SPLINE[1]

		for SAMPLE in range(1,(len(SPLINE[4:])/5)+1):

			BASE_COUNT_POS = 3 + (5*(SAMPLE-1))

			BASES_POS =  4 + (5*(SAMPLE-1))

			BASES = SPLINE[BASES_POS]

			BASES_COUNT = int(SPLINE[BASE_COUNT_POS])

			#remove beginning and end of read info
			BASES = re.sub(r'\^.','',BASES).replace("$","")

			#deal with Ns
			BASES_COUNT = BASES_COUNT - BASES.count("N") - BASES.count("*")

			BASES = BASES.replace("*","N").replace("N","")

			if BASES_COUNT in [0,-1]:

				TO_PRINT = TO_PRINT + " 0 0"

			elif BASES_COUNT == 1:

				if BASES.startswith(",") or BASES.startswith("."):

					TO_PRINT = TO_PRINT + " " + REF + " " + REF

				elif BASES.upper() in ["A","G","T", "C"]:

					TO_PRINT = TO_PRINT + " " + BASES[0].upper() + " " + BASES[0].upper()

				elif len(BASES) != 1 and ("+" in BASES or "-" in BASES):

					TO_PRINT = TO_PRINT + " " + BASES[0].upper() + " " + BASES[0].upper()

				else:

					TO_PRINT = TO_PRINT + " error"

					print BASES

					print BASES_COUNT

			else:

				if "+" in BASES or "-" in BASES:

					TO_REMOVE = []

					SWITCH = "OFF"

					INDEX = -1

					for CHAR in BASES:

						INDEX += 1

						if SWITCH == "ON":

							TO_REMOVE.append(INDEX)

							for INDEL_BASES in range(1,(int(CHAR)+2)):

								TO_REMOVE.append(INDEX - 1 + INDEL_BASES)

							SWITCH = "OFF"

						if CHAR == "+" or CHAR == "-":

							TO_REMOVE.append(INDEX)

					NEW_BASES = [BASES for BASE in BASES]

					for POS in sorted(TO_REMOVE, reverse=True):

						del(NEW_BASES[POS])

					BASES = " ".join(NEW_BASES)

				BASES = BASES.replace(",",REF).replace(".",REF)

				RANDOM_BASE = BASES[randint(0,len(BASES)-1)].upper()

				TO_PRINT = TO_PRINT + " " + RANDOM_BASE + " " + RANDOM_BASE

		print TO_PRINT

COUNT = 0

with open(SAMPLE_FILE) as FILE:

	for LINE in FILE:

		LINE = LINE.rstrip("\n")

		if COUNT == 0:

			f = open(SAMPLE_FILE.split(".")[0] + ".tfam", "w")

		else:

			f = open(SAMPLE_FILE.split(".")[0] + ".tfam", "a")

		f.write(LINE + " " + LINE + " 0 0 0 0\n")

		f.close()

		COUNT += 1

