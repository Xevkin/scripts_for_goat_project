import sys

import random

def most_common(lst):
    return max(set(lst), key=lst.count)

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		if LINE.startswith("#"):

			continue

		LINE = LINE.rstrip("\n")

		SPLINE = LINE.split("\t")

		REF = SPLINE[2]

		#skip sites with structural variants

		if bool([ele for ele in [">","<","+","*","[","]"] if(ele in SPLINE[4])]):

			continue

		# remove positional info
		SPLINE[4]=SPLINE[4].replace("^","").replace("$","")

		# randomly sample a base if one or two bases

		if int(SPLINE[3]) < 3:

			RANDBASE = SPLINE[4][random.randrange(0,len(SPLINE[4]))]

			if RANDBASE in [".",","]:

				BASE = REF

			else:

				BASE = RANDBASE.upper()

			print "\t".join([SPLINE[0],SPLINE[1],SPLINE[2],BASE])

		else:

			BASES = SPLINE[4].upper()

			BASES = BASES.replace(".",REF)

			BASES = BASES.replace(",",REF)

			print "\t".join([SPLINE[0],SPLINE[1],SPLINE[2],most_common(list(BASES))])
