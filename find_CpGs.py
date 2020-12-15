import sys

LINE_LENGTH = 0

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		LINE = LINE.rstrip("\n").upper()

		if LINE.startswith(">"):

			CHR = LINE.split(" ")[0].replace(">","")

			#count for position
			CHR_LINE = 1

			#update the last base
			LAST_BASE = "N"

			continue

		if LINE_LENGTH == 0:

			LINE_LENGTH = len(LINE)

		if (LINE.startswith("G") and LAST_BASE == "C") or (LINE.startswith("C") and LAST_BASE == "G"):

			print CHR + "\t" + str( (CHR_LINE-1)*LINE_LENGTH )

			print CHR + "\t" + str( ((CHR_LINE-1)*LINE_LENGTH) + 1)

		CURRENT_BASE = 1

		while len(LINE) > 1:

			if ("".join(LINE[0:2]) == "GC") or ("".join(LINE[0:2]) == "CG"):

				print CHR + "\t" + str( ((CHR_LINE-1)*LINE_LENGTH) + CURRENT_BASE )

				print CHR + "\t" + str( ((CHR_LINE-1)*LINE_LENGTH) + CURRENT_BASE + 1)

			LINE = LINE[1:]

			CURRENT_BASE += 1


		CHR_LINE += 1

		LAST_BASE = LINE[-1]
