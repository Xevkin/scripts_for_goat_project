import sys

import re

MATCH_COUNT = 0

MISMATCH_COUNT = 0

TRANS_MISMATCH_COUNT = 0

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		CALL="_"

		SPLINE = LINE.rstrip("\t").split("\t")

		if int(SPLINE[1]) >= 16616:

			continue

		BASES = SPLINE[4]

		BASES = re.sub("\^[ABCDEFGHIJKLMNO]","", BASES)

		BASES = re.sub("[+-]1[ACGTcgat]", "", BASES)

		BASES = re.sub("[+-]2[ACGTcgat][ACGTcgat]", "", BASES)

		BASES = re.sub("[+-]3[ACGTcgat][ACGTcgat][ACGTcgat]", "", BASES)

		BASES =	BASES.replace("$","").replace('*',"").replace(",",".").replace(".",SPLINE[2]).upper()

		if  len(BASES) < 3:

			continue

		A_count = BASES.count("A")

		T_count = BASES.count("T")

		C_count = BASES.count("C")

		G_count = BASES.count("G")

		if all(i < A_count for i in [T_count, C_count, G_count]):

			CALL = "A"

			MATCH_COUNT = MATCH_COUNT + A_count

			MISMATCH_COUNT = MISMATCH_COUNT + G_count + C_count + T_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + G_count

			continue

                if all(i < T_count for i in [A_count, C_count, G_count]):

                        CALL = "T"

                        MATCH_COUNT = MATCH_COUNT + T_count

                        MISMATCH_COUNT = MISMATCH_COUNT +  G_count + C_count + A_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + A_count

			continue

                if all(i < C_count for i in [T_count, A_count, G_count]):

                        CALl = "C"

                        MATCH_COUNT = MATCH_COUNT + C_count

                        MISMATCH_COUNT = MISMATCH_COUNT +  G_count + A_count + T_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + T_count

			continue

                if all(i < G_count for i in [T_count, C_count, A_count]):

                        CALL = "G"

                        MATCH_COUNT = MATCH_COUNT + G_count

                        MISMATCH_COUNT = MISMATCH_COUNT +  A_count + C_count + T_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + A_count

			continue

		BASES.replace(CALL, "")

		BASES_SET=set(list(BASES))

		if len(BASES_SET) > 1:

			if CALL != "_":

				NUCS = ["A",  "G", "T", "C"].remove(CALL)

			if len(BASES_SET) == 2:

				MATCH_COUNT = MATCH_COUNT + BASES.count(max(BASES))

				MISMATCH_COUNT = MISMATCH_COUNT + BASES.count(max(BASES))

				if BASES_SET not in [set(['A', 'G']), set(['C', 'T'])]:

					TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + BASES.count(max(BASES))

			elif len(BASES_SET) >= 3:

				MATCH_COUNT = MATCH_COUNT + BASES.count(max(BASES))

				BASES_REMOVED = BASES.replace(max(BASES),"")

				MISMATCH_COUNT = MISMATCH_COUNT + sum([BASES_REMOVED.count(i) for i in ["A",  "G", "T", "C"]])


	print sys.argv[1].split("/")[-1].split("_")[0] + "," + str(MATCH_COUNT) + "," + str(MISMATCH_COUNT) + "," +  str(MISMATCH_COUNT - TRANS_MISMATCH_COUNT)  + "," +  str(MATCH_COUNT + MISMATCH_COUNT)

