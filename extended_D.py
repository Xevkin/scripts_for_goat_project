#!/usr/bin/python

from __future__ import division

import sys

import gzip

import random

import numpy

import time
#reading each line should take 0.3sec per 1k for 10 tests, 1.5s/kline for 50 tests, 3.2s/klines for 100 tests

#input is gzipped haplo file, a space-seperated file with each line giving {H1, H2 group, H3 ind, list of H4s groups, grouped with _, the H5 outgroup groups} if single mode

BOOTS_LIST = []

#create the bootstrap regions
BOOTS=[]

with open("/home/kdaly/boots.list") as FILE:

	for LINE in FILE:

		SPLINE = LINE.rstrip("\n").split(" ")

		BOOTS.append([int(X) for X in SPLINE] + [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ])

HAPLO_FILE = sys.argv[1]

TESTS = []

with open(sys.argv[2]) as FILE:

	for LINE in FILE:

		TESTS.append(LINE.rstrip("\n").split(" "))

		BOOTS_LIST.append(BOOTS)

H1_GROUP_LIST = []

H1_INDICES_LIST = []

H2_GROUP_LIST = []

H2_INDICES_LIST = []

H3_LIST = []

H4_GROUPS_LIST = []

H4_LIST = []

H5_LIST = []

H4_INDICES_LIST = []

H5_INDICES_LIST = []

H4_GROUP_INDICES_LIST = []

ALL_INDS_LIST = []

SINGLE = sys.argv[3]

SINGLE_TOTAL_COUNT_LIST = []

SINGLE_B_COUNT_LIST = []

SINGLE_B_COUNT_TRANSV_LIST = []

TEST_COUNTER = -1

for TEST in TESTS:

	TEST_COUNTER += 1

	H1_GROUP_LIST.append(TEST[0].split(","))

	H1_INDICES_LIST.append([])

	H2_INDICES_LIST.append([])

	H2_GROUP_LIST.append(TEST[1].split(","))

	H3_LIST.append(TEST[2])

	H4_INDICES_LIST.append([])

	H4_GROUP_INDICES_LIST.append([])

	H5_INDICES_LIST.append([])

	H4_GROUPS_LIST.append([X.split(",") for X in TEST[3].split("_")])

	H4_LIST.append(TEST[3].replace("_",",").split(","))

	H5_LIST.append(TEST[4].split(","))

	ALL_INDS_LIST.append([])

	ALL_INDS_LIST[TEST_COUNTER].extend(H1_GROUP_LIST[TEST_COUNTER])

	ALL_INDS_LIST[TEST_COUNTER].extend(H2_GROUP_LIST[TEST_COUNTER])

	ALL_INDS_LIST[TEST_COUNTER].extend([H3_LIST[TEST_COUNTER]])

	ALL_INDS_LIST[TEST_COUNTER].extend(H5_LIST[TEST_COUNTER])

	ALL_INDS_LIST[TEST_COUNTER].extend(H4_LIST[TEST_COUNTER])

	SINGLE_TOTAL_COUNT_LIST.append(0)

	SINGLE_B_COUNT_LIST.append(0)

	SINGLE_B_COUNT_TRANSV_LIST.append(0)

with gzip.open(HAPLO_FILE) as FILE:

	for LINE in FILE:

		LINE = LINE.rstrip("\n")

		SPLINE = LINE.split("\t")

		TEST_COUNTER = -1

		for TEST in TESTS:

			TEST_COUNTER += 1

			if LINE.startswith("chr"):

				for IND in ALL_INDS_LIST[TEST_COUNTER]:

					if IND not in SPLINE:

						print IND + " is not in the header. Exiting..."

						raise SystemExit

				for H1_INDIV in H1_GROUP_LIST[TEST_COUNTER]:

					H1_INDICES_LIST[TEST_COUNTER].append(SPLINE.index(H1_INDIV))

				for H4_INDIV in H4_LIST[TEST_COUNTER]:

					H4_INDICES_LIST[TEST_COUNTER].append(SPLINE.index(H4_INDIV))

				for H2_INDIV in H2_GROUP_LIST[TEST_COUNTER]:

					H2_INDICES_LIST[TEST_COUNTER].append(SPLINE.index(H2_INDIV))

				for GROUP in H4_GROUPS_LIST[TEST_COUNTER]:

					H4_GROUP_INDICES_LIST[TEST_COUNTER].append([SPLINE.index(X) for X in GROUP])

				H3_INDEX = SPLINE.index(H3_LIST[TEST_COUNTER])

				for H5_INDIV in H5_LIST[TEST_COUNTER]:

					H5_INDICES_LIST[TEST_COUNTER].append(SPLINE.index(H5_INDIV))

			#now start the actual ABAAAA BAAAA etc calculations
			else:
				if  ([SPLINE[INDEX] for INDEX in [H3_INDEX]].count("N") == 0)  and ([SPLINE[INDEX] for INDEX in H5_INDICES_LIST[TEST_COUNTER] ].count("N") != len(H5_INDICES_LIST[TEST_COUNTER]))  and ([SPLINE[INDEX] for INDEX in H2_INDICES_LIST[TEST_COUNTER] ].count("N") != len(H2_INDICES_LIST[TEST_COUNTER])):

					if SINGLE not in ["yes", "Yes", "Single"]:

						if not ([SPLINE[INDEX] for INDEX in H1_INDICES_LIST[TEST_COUNTER] ].count("N") != len(H1_INDICES_LIST[TEST_COUNTER])):

							continue

					POSITIONS = [int(X) for X in SPLINE[0:2]]

					#filter on removing derived alleles present in our H4s
					H4_BASES = [SPLINE[INDEX] for INDEX in H4_INDICES_LIST[TEST_COUNTER]]

					H5_BASES = [SPLINE[INDEX] for INDEX in H5_INDICES_LIST[TEST_COUNTER]]

					H2_BASES = [SPLINE[INDEX] for INDEX in H2_INDICES_LIST[TEST_COUNTER]]

					if (len(''.join(set(H5_BASES)).replace("N","")) == 1):

						ANC_ALLELE = ''.join(set(H5_BASES)).replace("N","")

					else:

						continue

					DER_ALLELE = SPLINE[H3_INDEX]

					if ANC_ALLELE == DER_ALLELE:

						continue

					H1_BASES = [SPLINE[INDEX] for INDEX in H1_INDICES_LIST[TEST_COUNTER]]

					#randomly sample a H1 base
					H1_BASE = "N"

					if SINGLE not in ["yes", "Yes", "Single", "single"]:

						while H1_BASE == "N":

							H1_BASE = H1_BASES[random.randint(0,len(H1_BASES)-1)]

					H2_BASE = "N"

					while H2_BASE == "N":

						 H2_BASE = H2_BASES[random.randint(0,len(H2_BASES)-1)]

					#get the derived allele frequency of present bases
					DER_FREQ = H4_BASES.count(DER_ALLELE) /  len(list(filter(lambda a: a != "N", H4_LIST[TEST_COUNTER]) ) )

					if SINGLE in ["single", "Yes", "yes", "Single"]:

						MISSING = "NO"

						for GROUP in H4_GROUP_INDICES_LIST[TEST_COUNTER]:

							GROUP_BASES = [SPLINE[INDEX] for INDEX in GROUP]

							if GROUP_BASES.count("N") == len(GROUP):

								MISSING = "YES"

						if ( H4_BASES.count("N") != len(H4) )  and (MISSING == "NO"):

							SINGLE_TOTAL_COUNT_LIST[TEST_COUNTER] += 1

							if  (H2_BASE == DER_ALLELE) and (DER_FREQ == 0):

								SINGLE_B_COUNT_LIST[TEST_COUNTER] += 1

								if  ["GA", "AG", "CT", "TC"].count("".join([ANC_ALLELE,DER_ALLELE])) == 0:

									SINGLE_B_COUNT_TRANSV_LIST[TEST_COUNTER] += 1

						continue

					BOOT_COUNTER = -1

					#for each bootstrap region
					for BOOT in BOOTS:

						BOOT_COUNTER += 1

						BOOT_SPLIT = BOOT

						if not ( (BOOT[0] == POSITIONS[0]) and ( POSITIONS[1] >= BOOT[1]) and (POSITIONS[1] <= BOOT[2]) ):

							continue

						else:

						#if derived allele frequency is lower than cutoff

							COUNTER = -1

							for DER_CUTOFF in [0.00, 0.1, 0.2, 0.3, 0.4, 0.5,  0.6,  0.7,  0.8,  0.9, 1.0]:

								COUNTER += 1

								#this is so we can filter on no missing sites in our defined groups
								MISSING = "NO"

								for GROUP in H4_GROUP_INDICES_LIST[TEST_COUNTER]:

									GROUP_BASES = [SPLINE[INDEX] for INDEX in GROUP]

									#if any of the groups are missing data, skip
									if GROUP_BASES.count("N") == len(GROUP):

										MISSING = "YES"

								if ( H4_BASES.count("N") != len(H4_LIST[TEST_COUNTER]) ) and ( DER_FREQ <= DER_CUTOFF ) and (MISSING == "NO"):

									# AB sites
									if ( H1_BASE == ANC_ALLELE) and (H1_BASE != H2_BASE) and (H2_BASE == DER_ALLELE) :

										BOOTS_LIST[TEST_COUNTER][BOOT_COUNTER][3][COUNTER] += 1

										if ["GA", "AG", "CT", "TC"].count("".join([ANC_ALLELE,DER_ALLELE])) == 0:

											 BOOTS_LIST[TEST_COUNTER][BOOT_COUNTER][4][COUNTER]  += 1

									if ( H2_BASE == ANC_ALLELE) and (H1_BASE != H2_BASE) and (H1_BASE == DER_ALLELE):

										BOOTS_LIST[TEST_COUNTER][BOOT_COUNTER][5][COUNTER]  += 1

										if ["GA", "AG", "CT", "TC"].count("".join([ANC_ALLELE,DER_ALLELE])) == 0:

											 BOOTS_LIST[TEST_COUNTER][BOOT_COUNTER][6][COUNTER]  += 1

TEST_COUNTER = -1

for TEST in TESTS:

	TEST_COUNTER += 1

	if SINGLE in ["yes","single", "Single"]:

		print " _".join(H2_GROUP_LIST[TEST_COUNTER]) + " " +  " ".join([str(X) for X in [SINGLE_B_COUNT_LIST[TEST_COUNTER], SINGLE_B_COUNT_TRANSV_LIST[TEST_COUNTER], SINGLE_TOTAL_COUNT_LIST[TEST_COUNTER]]])

		exit()

	#list of the 11 different cutoffs and extended D estimates
	OUTPUT = [[] , [], [], [], [] , [], [], [], [] , [], []]

	TRANSV_OUTPUT = OUTPUT

	#calculate the final extended D too
	DE_OUT = [ [0, 0, 0, 0] , [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0] , [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0] , [0, 0, 0, 0],[0, 0, 0, 0]]

	DE_FINAL = []

	#for each bootstrap region
	for BLOCK in BOOTS_LIST[TEST_COUNTER]:

		#for each cutoff level
		for CUTOFF in range(0,11):

			DE_OUT[CUTOFF][0] += BLOCK[3][CUTOFF]

			DE_OUT[CUTOFF][1] += BLOCK[4][CUTOFF]

			DE_OUT[CUTOFF][2] += BLOCK[5][CUTOFF]

			DE_OUT[CUTOFF][3] += BLOCK[6][CUTOFF]

	for CUTOFF in range(0, 11):

		#skip blocks with no data
		if DE_OUT[CUTOFF][0] == 0 and DE_OUT[CUTOFF][2] == 0:

			continue

		else:

			EXTENDED_D = (DE_OUT[CUTOFF][0] - DE_OUT[CUTOFF][2]) / (DE_OUT[CUTOFF][0] + DE_OUT[CUTOFF][2])

			EXTENDED_D_TRANSV = (DE_OUT[CUTOFF][1] - DE_OUT[CUTOFF][3]) / (DE_OUT[CUTOFF][1] + DE_OUT[CUTOFF][3])

			DE_FINAL.append([EXTENDED_D, EXTENDED_D_TRANSV])

	#1000 replicates
	for REPLICATE in range(0,1000):

		#keep track of the subsampled blocks
		REPLICATE_BLOCK=[]

		AB_TOTAL=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		BA_TOTAL=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		AB_TRANSV_TOTAL=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		BA_TRANSV_TOTAL=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		#build a subsample of the same number of total bootstrap regions
		for SUBSAMPLE in range(0, len(BOOTS_LIST[TEST_COUNTER])):

			#randomly sample with replacement
			SUBSAMPLE_INDEX = random.randrange(0, len(BOOTS_LIST[TEST_COUNTER]))

		REPLICATE_BLOCK = BOOTS[SUBSAMPLE_INDEX]

		#continue to resample if the selection block has no data
		while (REPLICATE_BLOCK[3][CUTOFF] == 0) and  (REPLICATE_BLOCK[5][CUTOFF] == 0):

			REPLICATE_BLOCK = BOOTS[random.randrange(0, len(BOOTS))]

		#now generate total of AB BA sites for each of the 11 cutoffs
		for CUTOFF in range(0,11):

			AB_TOTAL[CUTOFF] += REPLICATE_BLOCK[3][CUTOFF]

			AB_TRANSV_TOTAL[CUTOFF] += REPLICATE_BLOCK[4][CUTOFF]

			BA_TOTAL[CUTOFF] += REPLICATE_BLOCK[5][CUTOFF]

			BA_TRANSV_TOTAL[CUTOFF] += REPLICATE_BLOCK[6][CUTOFF]

		for CUTOFF in range(0,11):

			if (AB_TOTAL[CUTOFF] + BA_TOTAL[CUTOFF]) == 0:

				EXTENDED_D = "NA"

			else:

				EXTENDED_D = (AB_TOTAL[CUTOFF] - BA_TOTAL[CUTOFF]) / (AB_TOTAL[CUTOFF] + BA_TOTAL[CUTOFF])

			if (AB_TRANSV_TOTAL[CUTOFF] + BA_TRANSV_TOTAL[CUTOFF]) == 0:

				EXTENDED_D_TRANSV = "NA"

			else:

				EXTENDED_D_TRANSV = (AB_TRANSV_TOTAL[CUTOFF] - BA_TRANSV_TOTAL[CUTOFF]) / (AB_TRANSV_TOTAL[CUTOFF] + BA_TRANSV_TOTAL[CUTOFF])

				OUTPUT[CUTOFF].append(EXTENDED_D)

				OUTPUT[CUTOFF].append(EXTENDED_D_TRANSV)

	TO_PRINT = []

	for CUTOFF in range(0, 11):

		DE_FINAL_VALUE = DE_FINAL[CUTOFF][0]

		DE_TRANSV_FINAL_VALUE = DE_FINAL[CUTOFF][1]

		STD = numpy.std(numpy.array(OUTPUT[CUTOFF]))

		STD_TRANSV = numpy.std(numpy.array(TRANSV_OUTPUT[CUTOFF]))

		if STD == 0:

			Z = "NA"

		else:

			 Z = DE_FINAL_VALUE / STD

		if STD_TRANSV == 0:

			Z_TRANSV = "NA"

		else:

			Z_TRANSV = DE_TRANSV_FINAL_VALUE / STD_TRANSV

		TO_PRINT.extend([str(X) for X in [DE_FINAL_VALUE, STD, Z, DE_OUT[CUTOFF][0], DE_OUT[CUTOFF][2], DE_TRANSV_FINAL_VALUE, STD_TRANSV, Z_TRANSV, DE_OUT[CUTOFF][1], DE_OUT[CUTOFF][3]]])

	print " ".join(TEST) + " " + " ".join(TO_PRINT)
