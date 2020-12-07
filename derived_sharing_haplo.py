from __future__ import division

import sys

import numpy as np

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

	#jackknifing
	#Going with  5Mb JKs
	#assuming bootstrap regions are found at /home/kdaly/bootstrap.list
	#create a list storing the bootstrap region information

	JACKKNIFES = []
	with open("/home/kdaly/bootstrap.list") as FILE:

		for LINE in FILE:

			LINE = LINE.rstrip("\n")

			REGION = [0,0,0,0,0]

			REGION[0] = int(LINE.split(":")[0])

			REGION[1] = int(LINE.split(":")[1].split("-")[0])

			REGION[2] = int(LINE.split(":")[1].split("-")[1])

			JACKKNIFES.append(REGION)

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

			#keep jack-knifes estimates
			#create a counter to loop through each block
			COUNTER = 0

			#for each line, interate over the jackknife regions
			for IND_JACKKNIFE in JACKKNIFES:

				#if the position matches the current jk region
				if CHR == IND_JACKKNIFE[0] and POS >= IND_JACKKNIFE[1] and POS <= IND_JACKKNIFE[2]:

					#if we are doing a 4 pop test:
					if TEST == 4:

						#add to the current ABBA BABA count for the current jk region
						JACKKNIFES[COUNTER][3],JACKKNIFES[COUNTER][4] = GET_ABBA_BABA(SPLINE,H1,H2,H3,JACKKNIFES[COUNTER][3],JACKKNIFES[COUNTER][4])

						#also add to the ABBA BABA lists
						ABBA[COUNTER] =  JACKKNIFES[COUNTER][3]

						BABA[COUNTER] =  JACKKNIFES[COUNTER][4]

					elif TEST == 5:

						#add to ABABA BAABA counter
						JACKKNIFES[COUNTER][3],JACKKNIFES[COUNTER][4] = GET_ABBA_BABA(SPLINE,H1,H2,H3,JACKKNIFES[COUNTER][3],JACKKNIFES[COUNTER][4],H4)

						ABABA[COUNTER] = JACKKNIFES[COUNTER][3]

						BAABA[COUNTER] = JACKKNIFES[COUNTER][4]

					continue

				#update the counter to keep track of the current jk region
				COUNTER += 1

	#new counter, for tracking jk regions as we filter out those with no data
	COUNTER = 0

	#we want to exclude jk regions with no information
	TO_INCLUDE = []

	for BLOCK in JACKKNIFES:

		#if missing both ABBA or BABA information, skip
		if BLOCK[3] != 0 or BLOCK[4] != 0:

			TO_INCLUDE.append(COUNTER)

		#update counter for new jk region
		COUNTER += 1

	#strip out ABBA / BABA information from jk regions we are keeping; also sum

	ABBA_TO_INCLUDE = list(map(ABBA.__getitem__,TO_INCLUDE))

	ABBA_TO_INCLUDE_SUM = sum(ABBA_TO_INCLUDE)

	BABA_TO_INCLUDE = list(map(BABA.__getitem__,TO_INCLUDE))

	BABA_TO_INCLUDE_SUM = sum(BABA_TO_INCLUDE)

	#if we are doing a 5 pop test, do the same for ABABA BAABA
	if TEST == 5:

		ABABA_TO_INCLUDE = list(map(ABABA.__getitem__,TO_INCLUDE))

		ABABA_TO_INCLUDE_SUM = sum(ABABA_TO_INCLUDE)

		BAABA_TO_INCLUDE = list(map(BAABA.__getitem__,TO_INCLUDE))

		BAABA_TO_INCLUDE_SUM = sum(BAABA_TO_INCLUDE)

	#new counter for jk regions we are keeping
	COUNTER = 0

	#list for each jk estimate of D
	JACK_D = []

	#if we are doing a four pop test
	if TEST == 4:

		#for each jk block
		for BLOCK in ABBA_TO_INCLUDE:

			#"remove" a given region by adding its negative to the list, which will be summed
			ABBA_SUBSET = ABBA_TO_INCLUDE + [-ABBA_TO_INCLUDE[COUNTER]]

			BABA_SUBSET = BABA_TO_INCLUDE + [-BABA_TO_INCLUDE[COUNTER]]

			#add the jk D estimate to the jk D list
			JACK_D.append(D(sum(ABBA_SUBSET), sum(BABA_SUBSET) ))

			#update counter for each jk region
			COUNTER += 1

		#get the total D estimate
		#print ABBA_TO_INCLUDE_SUM

		#print BABA_TO_INCLUDE_SUM

		D_EST = D(ABBA_TO_INCLUDE_SUM, BABA_TO_INCLUDE_SUM)

		#calculate the sd of the jk D estimates
		JACK_ERR = np.std(JACK_D)

		#calculate Z score of the D estimate
		Z = D_EST / JACK_ERR

		#print out the total ABBA, BABA sites, D estimates, jk error, and the Z score
		print " ".join(str(e) for e in [ABBA_TO_INCLUDE_SUM, BABA_TO_INCLUDE_SUM, D_EST, JACK_ERR, Z])

	if TEST == 5:

		for BLOCK in ABBA_TO_INCLUDE:

			ABABA_SUBSET = ABABA_TO_INCLUDE + [-ABABA_TO_INCLUDE[COUNTER]]

			BAABA_SUBSET = BAABA_TO_INCLUDE + [-BAABA_TO_INCLUDE[COUNTER]]

			JACK_D.append(D(sum(ABABA_SUBSET), sum(BAABA_SUBSET) ))

			COUNTER += 1

		D_EST = D(ABABA_TO_INCLUDE_SUM, BAABA_TO_INCLUDE_SUM)

		JACK_ERR = np.std(JACK_D)

		Z = D_EST / JACK_ERR

		print " ".join(str(e) for e in [ABABA_TO_INCLUDE_SUM, BAABA_TO_INCLUDE_SUM, D_EST, JACK_ERR, Z])

def GET_ABBA_BABA(SPLINE_IN,H1_IN,H2_IN,H3_IN,NABBA_IN, NBABA_IN,H4_IN="none"):

	#print [SPLINE_IN,H1_IN,H2_IN,H3_IN,NABBA_IN, NBABA_IN,H4_IN]

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
