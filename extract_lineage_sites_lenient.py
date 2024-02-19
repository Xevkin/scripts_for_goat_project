from __future__ import division

import sys

# python2 extract_lineage_sites.py core.vcf X,Y,Z

# X Y Z must be the base 1 column they appear in the vcf

COLUMNS = [int(x) - 1 for x in sys.argv[2].split(",")]

with open(sys.argv[1]) as VCF_FILE:

	for LINE in VCF_FILE:

		SPLINE = LINE.split()

		# get the number of samples in the vcf

		if LINE.startswith("#C"):

			SAMPLE_COUNT = len(SPLINE) - 9

			continue

		if LINE.startswith("#"):

			continue

		# snp sites only

		if SPLINE[7] != "TYPE=snp":

			continue

		# bialleleic sites only

		if SPLINE[4].count(",") > 0:

			continue

		#create lists to store our variables
		TARGET = []

		BACKGROUND = []

		# for each sample, check if the sample is in the columns of interest

		# remember the sample columns start in column 9 (base 0)

		for SAMPLE in range(9,SAMPLE_COUNT + 9):

			if SAMPLE in COLUMNS:

				# if so, add to the target list for this site

				TARGET.append(SPLINE[SAMPLE])

			else:

				# otherwise, add to the background list for this site
				BACKGROUND.append(SPLINE[SAMPLE])

		# allow either one TARGET or one BACKGROUND site to match

		most_common_target = max(set(TARGET),key=TARGET.count)

		most_common_back = max(set(BACKGROUND),key=BACKGROUND.count)

		if (  ((TARGET.count(most_common_target) == (len(TARGET) - 1 )) and  (BACKGROUND.count(BACKGROUND[0]) == len(BACKGROUND)) ) or ((TARGET.count(TARGET[0]) == len(TARGET)) and  (BACKGROUND.count(most_common_back) == (len(BACKGROUND)-1)) ) ) and (TARGET[0] != BACKGROUND[0]):

			print " ".join([SPLINE[x] for x in [0,1,3,4]]) + " " + str(SPLINE[int(TARGET[TARGET.index(most_common_target)])+3])
