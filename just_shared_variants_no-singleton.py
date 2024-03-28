from __future__ import division

import sys

import gzip

#inputs are gzip_haplo_file sample_of_interest comma_seperated_list_of_comparisons
SOI = sys.argv[2]

SAMPLES = sys.argv[3].rstrip("\n").split(",")

SITES_TOTAL = []

DIFFER_TOTAL = []

DIFFER_TRANSV_TOTAL = []

with gzip.open(sys.argv[1]) as FILE:

	for LINE in FILE:

		LINE=LINE.rstrip("\n")

		SPLINE = LINE.split("\t")

		if LINE.startswith("chr"):

			SAMPLE_INDEXES = [SPLINE.index(X) for X in SAMPLES]

			if SOI not in SPLINE:

				print "Sample of interest " + SOI + " is not in header"

				raise SystemExit

			else:

				SOI_INDEX = SPLINE.index(SOI)

			for SAMPLE in SAMPLES:

				if SAMPLE not in SPLINE:

					print SAMPLE + " not in header"

					raise SystemExit

				else:

					SITES_TOTAL.append(0)

					DIFFER_TOTAL.append(0)

					DIFFER_TRANSV_TOTAL.append(0)

			continue

		else:

			if SPLINE[SOI_INDEX] == "N":

				continue

			else:

				COUNTER = -1

				for SAMPLE in SAMPLES:

					COUNTER += 1

					if not SPLINE[SAMPLE_INDEXES[COUNTER]] == "N":

						SITES_TOTAL[COUNTER] += 1

						#controlling for error in the sample/gat genome (not Caprid as they will hav eunique ancestry)
						if SPLINE[SAMPLE_INDEXES[COUNTER]] != SPLINE[SOI_INDEX] and (SPLINE[3:].count(SPLINE[SAMPLE_INDEXES[COUNTER]]) > 1):

							DIFFER_TOTAL[COUNTER] += 1

							if SPLINE[SAMPLE_INDEXES[COUNTER]]+SPLINE[SOI_INDEX] not in ["GA","AG", "CT", "TC"]:

								DIFFER_TRANSV_TOTAL[COUNTER] += 1

COUNTER = -1

for SAMPLE in SAMPLES:

	COUNTER += 1

	print SAMPLE + " "+   " ".join([str(X) for X in [SITES_TOTAL[COUNTER], DIFFER_TOTAL[COUNTER], DIFFER_TRANSV_TOTAL[COUNTER]]])

