from __future__ import division

import sys

import numpy as np

import random

def main():

	#take samples to test
	#going by order in sample file, outgroup in position one
	#H1,H2,H3
	SAMPLES = sys.argv[2].rstrip("\n").split(",")

	#check if we are doing a 4/5 pop test
	TEST=len(SAMPLES)+1

	#print out the samples
	H1 = int(SAMPLES[0])+2
	H2 = int(SAMPLES[1])+2
	H3 = int(SAMPLES[2])+2
	if TEST == 5: H4 = int(SAMPLES[3])+2

	#lists to store the ABBA BABA sites
	ABBA = []
	BABA = []
	ABABA = []
	BAABA = []

	#bootstrapping
	#Going with 5Mb regions
	#assuming bootstrap regions are found at /home/kdaly/bootstrap.list
	#create a list storing the bootstrap region information

	BOOTSTRAPS = []
	with open("/home/kdaly/bootstrap.list") as FILE:

		for LINE in FILE:

			LINE = LINE.rstrip("\n")

			REGION = [0,0,0,0,0]

			REGION[0] = int(LINE.split(":")[0])

			REGION[1] = int(LINE.split(":")[1].split("-")[0])

			REGION[2] = int(LINE.split(":")[1].split("-")[1])

			BOOTSTRAPS.append(REGION)

			ABBA.append(0)

			BABA.append(0)

			ABABA.append(0)

			BAABA.append(0)

	#go through each position in the haplo file
	with open(sys.argv[1]) as FILE:

		for LINE in FILE:

			LINE = LINE.rstrip("\n")

			SPLINE = LINE.split("\t")

			#skip over header, printing the samples of interest
			if LINE.startswith("chr"):

				print "H1 is " + SPLINE[H1]

				print "H2 is " + SPLINE[H2]

				print "H3 is " + SPLINE[H3]

				if TEST == 5: print "H4 is " + SPLINE[H4]

				continue

			#skip if missing ancestral information
			#ancestral genome MUST be in position 1
			if SPLINE[3] == "N":

				continue

			#get position and chromosome for later
			CHR = int(SPLINE[0])

			POS = int(SPLINE[1])

			#keep bootstrap estimates
			#create a counter to loop through each block
			COUNTER = 0

			#for each line, interate over the bs regions
			for BOOTSTRAP in BOOTSTRAPS:

				#if the position matches the current bs region
				if CHR == BOOTSTRAP[0] and POS >= BOOTSTRAP[1] and POS <= BOOTSTRAP[2]:

					#if we are doing a 4 pop test:
					if TEST == 4:

						#add to the current ABBA BABA count for the current bs region
						BOOTSTRAPS[COUNTER][3],BOOTSTRAPS[COUNTER][4] = GET_ABBA_BABA(SPLINE,H1,H2,H3,BOOTSTRAPS[COUNTER][3],BOOTSTRAPS[COUNTER][4])

						#also add to the ABBA BABA lists
						ABBA[COUNTER] =  BOOTSTRAPS[COUNTER][3]

						BABA[COUNTER] =  BOOTSTRAPS[COUNTER][4]

					elif TEST == 5:

						#add to ABABA BAABA counter
						BOOTSTRAPS[COUNTER][3],BOOTSTRAPS[COUNTER][4] = GET_ABBA_BABA(SPLINE,H1,H2,H3,BOOTSTRAPS[COUNTER][3],BOOTSTRAPS[COUNTER][4],H4)

						ABABA[COUNTER] = BOOTSTRAPS[COUNTER][3]

						BAABA[COUNTER] = BOOTSTRAPS[COUNTER][4]

					continue

				#update the counter to keep track of the current bs region
				COUNTER += 1

	#new counter, for tracking bs regions as we filter out those with no data
	COUNTER = 0

	#we want to exclude bs regions with no information
	TO_INCLUDE = []

	for BLOCK in BOOTSTRAPS:

		#if missing both ABBA or BABA information, skip
		if BLOCK[3] != 0 or BLOCK[4] != 0:

			TO_INCLUDE.append(COUNTER)

		#update counter for new bs region
		COUNTER += 1

	#strip out ABBA / BABA information from bs regions we are keeping; also sum

	ABBA_TO_INCLUDE = list(map(ABBA.__getitem__,TO_INCLUDE))

	ABBA_TO_INCLUDE_SUM = sum(ABBA_TO_INCLUDE)

	BABA_TO_INCLUDE = list(map(BABA.__getitem__,TO_INCLUDE))

	BABA_TO_INCLUDE_SUM = sum(BABA_TO_INCLUDE)

	BOOT_NUM = len(TO_INCLUDE)

	#if we are doing a 5 pop test, do the same for ABABA BAABA
	if TEST == 5:

		ABABA_TO_INCLUDE = list(map(ABABA.__getitem__,TO_INCLUDE))

		ABABA_TO_INCLUDE_SUM = sum(ABABA_TO_INCLUDE)

		BAABA_TO_INCLUDE = list(map(BAABA.__getitem__,TO_INCLUDE))

		BAABA_TO_INCLUDE_SUM = sum(BAABA_TO_INCLUDE)

	#new counter for jk regions we are keeping
	COUNTER = 0

	#list for each bs estimate of D
	BOOT_D = []

	#if we are doing a four pop test
	if TEST == 4:

		#for each bs block
		for BLOCK in ABBA_TO_INCLUDE:

			ABBA_SUBSET = []

			BABA_SUBSET = []

			#generate a new set of regions by sampling with replacement
			#1000 bootstrap replicates
			for I in range(0,1000):

				RANDOM_INDEX = random.randrange(0,BOOT_NUM)

				ABBA_SUBSET.append(ABBA_TO_INCLUDE[RANDOM_INDEX])

				BABA_SUBSET.append(BABA_TO_INCLUDE[RANDOM_INDEX])

			#add the bs D estimate to the bs D list
			BOOT_D.append(D(sum(ABBA_SUBSET), sum(BABA_SUBSET) ))

			#update counter for each bs region
			COUNTER += 1

		#get the total D estimate

		D_EST = D(ABBA_TO_INCLUDE_SUM, BABA_TO_INCLUDE_SUM)

		#calculate the sd of the jk D estimates
		BOOT_ERR = np.std(BOOT_D)

		#calculate Z score of the D estimate
		Z = D_EST / BOOT_ERR

		#print out the total ABBA, BABA sites, D estimates, jk error, and the Z score
		print " ".join(str(e) for e in [ABBA_TO_INCLUDE_SUM, BABA_TO_INCLUDE_SUM, D_EST, BOOT_ERR, Z])

	if TEST == 5:

		for BLOCK in ABBA_TO_INCLUDE:

			ABABA_SUBSET = []

			BAABA_SUBSET = []

			#1000 bootstrap replicates
			for I in range(0,1000):

				RANDOM_INDEX = random.randrange(0,BOOT_NUM)

				ABABA_SUBSET.append(ABABA_TO_INCLUDE[RANDOM_INDEX])

				BAABA_SUBSET.append(BAABA_TO_INCLUDE[RANDOM_INDEX])

			BOOT_D.append(D(sum(ABABA_SUBSET), sum(BAABA_SUBSET) ))

			COUNTER += 1

		D_EST = D(ABABA_TO_INCLUDE_SUM, BAABA_TO_INCLUDE_SUM)

		BOOT_ERR = np.std(BOOT_D)

		Z = D_EST / BOOT_ERR

		print " ".join(str(e) for e in [ABABA_TO_INCLUDE_SUM, BAABA_TO_INCLUDE_SUM, D_EST, BOOT_ERR, Z])

def GET_ABBA_BABA(SPLINE_IN,H1_IN,H2_IN,H3_IN,NABBA_IN, NBABA_IN,H4_IN="none"):

	" Function for taking in a .haplo row and determining if ABBA/BABA count should be added to"
	ANCESTRAL_IN = SPLINE_IN[3]

	DERIVED_IN = "".join(SPLINE_IN[4:]).replace(ANCESTRAL_IN,"").replace("N","")[0]

	if H4_IN == "none":

		#add to the current ABBA/BABA count when appropriate
		if SPLINE_IN[H1_IN] == ANCESTRAL_IN and SPLINE_IN[H2_IN] == DERIVED_IN and SPLINE_IN[H2_IN] == SPLINE_IN[H3_IN]:

			NABBA_IN += 1

		elif SPLINE_IN[H2_IN] == ANCESTRAL_IN and SPLINE_IN[H1_IN] == DERIVED_IN and SPLINE_IN[H1_IN] == SPLINE_IN[H3_IN]:

			NBABA_IN += 1

	else:

		if SPLINE_IN[H1_IN] == ANCESTRAL_IN and SPLINE_IN[H2_IN] == DERIVED_IN and SPLINE_IN[H1_IN] == SPLINE_IN[H3_IN] and SPLINE_IN[H2_IN] == SPLINE_IN[H4_IN]:

			#NABBA is a stand in for ABABA
			NABBA_IN += 1

		if SPLINE_IN[H2_IN] == ANCESTRAL_IN and SPLINE_IN[H1_IN] == DERIVED_IN and SPLINE_IN[H2_IN] == SPLINE_IN[H3_IN] and SPLINE_IN[H1_IN] == SPLINE_IN[H4_IN]:

			#again, NBABA is a stand in for BAABA
			NBABA_IN += 1

	return(NABBA_IN, NBABA_IN)

def D(ABBA, BABA):
	" Function for calculating D statistic "
	return( (ABBA - BABA) / (ABBA + BABA))

main()
