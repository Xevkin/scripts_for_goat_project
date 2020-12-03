from __future__ import division

import sys

import numpy as np

def main():

	#take samples to test
	#going by order in sample file, outgroup in position one
	#H1,H2,H3
	SAMPLES = sys.argv[2].rstrip("\n").split(",")

	H1 = int(SAMPLES[0])+2
	H2 = int(SAMPLES[1])+2
	H3 = int(SAMPLES[2])+2

	NABBA = 0
	NBABA = 0

	D_LIST = []

	#jackknifing
	#Going with  5Mb JKs
	#assuming bootstrap regions are found at /home/kdaly/bootstrap.list

	JACKKNIFES = []

	with open("/home/kdaly/bootstrap.list") as FILE:

		for LINE in FILE:

			LINE = LINE.rstrip("\n")

			REGION = [0,0,0,0,0]

			REGION[0] = int(LINE.split(":")[0])

			REGION[1] = int(LINE.split(":")[1].split("-")[0])

			REGION[2] = int(LINE.split(":")[1].split("-")[1])

			JACKKNIFES.append(REGION)

	with open(sys.argv[1]) as FILE:

		for LINE in FILE:

			LINE = LINE.rstrip("\n")

			SPLINE = LINE.split("\t")

			if LINE.startswith("chr"):

				print "H1 is " + SPLINE[H1]

				print "H2 is " + SPLINE[H2]

				print "H3 is " + SPLINE[H3]

				continue

			#skip if missing ancestral information
			#ancestral genome MUST be in position 1
			if SPLINE[3] == "N":

				continue

			#get ancestral and derived allele
			ANCESTRAL = SPLINE[3]

			DERIVED = "".join(SPLINE[4:]).replace(ANCESTRAL,"").replace("N","")[0]

			CHR = int(SPLINE[0])

			POS = int(SPLINE[1])

			#get ABBA/BABA sites
			NABBA, NBABA = GET_ABBA_BABA(SPLINE,ANCESTRAL,DERIVED,H1,H2,H3,NABBA, NBABA)

			for IND_JACKKNIFE in JACKKNIFES:

				if not(CHR == IND_JACKKNIFE[0] and POS >= IND_JACKKNIFE[1] and POS <= IND_JACKKNIFE[2]):

					IND_JACKKNIFE[3], IND_JACKKNIFE[4] = GET_ABBA_BABA(SPLINE,ANCESTRAL,DERIVED,H1,H2,H3,IND_JACKKNIFE[3], IND_JACKKNIFE[4])

				#print "EXCLUDE" + " " + str(CHR) + " " + str(IND_JACKKNIFE[1]) + " " + str(IND_JACKKNIFE[2])

	print str(NABBA) + " " + str(NBABA) + str(D(NABAB, NBABA))

	#get standard deviation of estimates and print
	for m in JACKKNIFES:

		 D_LIST.append( D(m[3], m[4]))

	print np.std(D_LIST)

def GET_ABBA_BABA(SPLINE_IN,ANCESTRAL_IN,DERIVED_IN,H1_IN,H2_IN,H3_IN,NABBA_IN, NBABA_IN):

	if SPLINE_IN[H1_IN] == ANCESTRAL_IN and SPLINE_IN[H2_IN] == DERIVED_IN and SPLINE_IN[H2_IN] == SPLINE_IN[H3_IN]:

		NABBA_IN += 1

	if SPLINE_IN[H2_IN] == ANCESTRAL_IN and SPLINE_IN[H1_IN] == DERIVED_IN and SPLINE_IN[H1_IN] == SPLINE_IN[H3_IN]:

		NBABA_IN += 1

	return(NABBA_IN, NBABA_IN)

def D(ABBA, BABA):

	return( (ABBA - BABA) / (ABBA + BABA))

main()
