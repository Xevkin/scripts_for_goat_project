from __future__ import division
import sys

#input is tped with ONE individual

WINDOW_SIZE = int(sys.argv[2].rstrip("\n"))

STEP_SIZE = int(WINDOW_SIZE) / 2

HET_COUNT = 0

CALL_COUNT = 0

CUR_CHR = 0

#switch for start of chromosomes
CHR_SWITCH = 0

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		LINE = LINE.rstrip("\n")

		SPLINE = LINE.split(" ")

		#current position and chromosome

		CHR = int(SPLINE[0])

		POS = int(SPLINE[3])

                #get genotypes of each individual
                IND_VAR = [SPLINE[4], SPLINE[5]]

		#update start of window if at new window

		if (CALL_COUNT == 0) and (IND_VAR != ["0", "0"]):

			WIND_START = POS

		#if we are on a new chromosome and there is data, update
		if (IND_VAR != ["0", "0"]) and (CHR > CUR_CHR ):

			#print the last window
			if (CHR != 1):

				#get the midpoint of window
	                        WIND_END = POS

        	                MIDPOINT = (WIND_START + WIND_END) / 2

                	        #get the het rate for current window
	                        HET_RATE = float((HET_COUNT + PREV_HET_COUNT) / (STEP_SIZE + CALL_COUNT) )

        	                print " ".join(str(x) for x in [CUR_CHR, MIDPOINT, HET_RATE])

                        #reset the windows
                        PREV_HET_COUNT = HET_COUNT

                        HET_COUNT = 0

                        CALL_COUNT = 0

			CUR_CHR = CHR

			WIND_START = POS

		#miissing site
		if (IND_VAR == ["0","0"]):

			continue

		#het site
		elif (IND_VAR[0] != IND_VAR[1] ):

			HET_COUNT += 1

			CALL_COUNT += 1

		else:

			CALL_COUNT += 1

		if CALL_COUNT == STEP_SIZE:

			if (CHR_SWITCH == 0):

				CHR_SWITCH = 1

				#reset the windows
				PREV_HET_COUNT = HET_COUNT

				HET_COUNT = 0

				CALL_COUNT = 0

				continue

			#get the midpoint of window
			WIND_END = POS

			MIDPOINT = (WIND_START + WIND_END) / 2
			#get the het rate for current window
			HET_RATE = float((HET_COUNT + PREV_HET_COUNT) / WINDOW_SIZE)

			print " ".join(str(x) for x in [CUR_CHR, MIDPOINT, HET_RATE])

			#reset the windows
			PREV_HET_COUNT = HET_COUNT

			HET_COUNT = 0

			CALL_COUNT = 0





