import sys

CURRENT_SAMPLE = ""

PREVIOUS_ROH_STATE = "A_not_roh"

PREVIOUS_CHR = 0

PREVIOUS_WINDOW_START = -2

PREVIOUS_WINDOW_END = -2

ROH_COUNTER = 0

CURRENT_SAMPLE = "N\A"

SMUDGE = int(sys.argv[2])

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		SPLINE = LINE.rstrip("\n").split(" ")

		CURRENT_ROH_STATE = SPLINE[10]

		CURRENT_CHR = int(SPLINE[0])

		CURRENT_WINDOW_START = int(SPLINE[1])

		CURRENT_WINDOW_END = int(SPLINE[2])

		if SPLINE[6] != CURRENT_SAMPLE:

			if CURRENT_SAMPLE != "N\A":

				if PREVIOUS_ROH_STATE == "B_roh":

					ROH_COUNTER += 1

				print CURRENT_SAMPLE + " " + str(ROH_COUNTER)

			 #if we are at a new sample, update and move on

			CURRENT_SAMPLE = SPLINE[6]

			PREVIOUS_ROH_STATE = CURRENT_ROH_STATE

			PREVIOUS_CHR = CURRENT_CHR

			PREVIVOUS_WINDOW_END = CURRENT_WINDOW_END

			PREVIOUS_WINDOW_START = CURRENT_WINDOW_START

			ROH_COUNTER = 0

			continue

		if (CURRENT_ROH_STATE == "B_roh"):

			#check in case we are not in a new ROH due to the consequetive break

			if CURRENT_CHR != PREVIOUS_CHR:

				ROH_COUNTER += 1

			if (CURRENT_CHR == PREVIOUS_CHR) and (CURRENT_WINDOW_START != PREVIVOUS_WINDOW_END+1) :

				#taking into account the SMUDGE factor of 1Mb
				if (CURRENT_WINDOW_START > (PREVIVOUS_WINDOW_END + SMUDGE + 1)):

					ROH_COUNTER += 1

		#if we are not in a ROH
		else:
			#check in case we are no longer in a ROH
			if PREVIOUS_ROH_STATE == "B_roh":

				ROH_COUNTER += 1

		PREVIOUS_ROH_STATE = CURRENT_ROH_STATE

		PREVIOUS_CHR = CURRENT_CHR

		PREVIVOUS_WINDOW_END = CURRENT_WINDOW_END

		PREVIOUS_WINDOW_START = CURRENT_WINDOW_START

		CURRENT_SAMPLE = SPLINE[6]


#catch any final ROH
if PREVIOUS_ROH_STATE == "B_roh":

	ROH_COUNTER += 1


#print the last sample
print CURRENT_SAMPLE + " " + str(ROH_COUNTER)



