from __future__ import division

import sys

import gzip

IBS_FILE = sys.argv[1]

BOOT_FILE = sys.argv[2]

BOOT_LIST = []

#list to store metrics for each pairwise
PAIRS = []

#create a list with the boot strap regions
with open(BOOT_FILE) as FILE:

	for LINE in FILE:

		BOOT_LIST.append(LINE.rstrip("\n").split(" "))

with gzip.open(IBS_FILE, "r") as FILE:

	for LINE in FILE:

		#counter to keep track of which pair we are looking at
		COUNTER1 = 0

		SPLINE = LINE.split("\t")

		#check if header line, build a list of pairwise comparisons
		if LINE.startswith("chr"):

			#get the number of innds from the last sample name
			IND_NUM = int(SPLINE[-1].replace("ind", "")) + 1

			#for all unique pairwise comparisons
			for IND in range(0,IND_NUM):

				for IND2 in range(IND, IND_NUM):

					if not (IND2 == IND):

						#create a list storing the ind #s for the pair members, and holders for #sites and #matches
						PAIRS.append([IND ,IND2, 0, 0])

		else:
			#check if variant is in bootstrap windows
			CHECK = "NO"

			for BOOT in BOOT_LIST:

				if (SPLINE[0] == BOOT[0]) and ( int(SPLINE[1]) >= int(BOOT[1]) ) and ( int(SPLINE[1]) <= int(BOOT[2]) ):

					CHECK = "YES"

					continue

			if CHECK == "YES":

				#for each non-header line, go through each pair
				for PAIR in PAIRS:

					PAIR1 = PAIR[0]

					PAIR2 = PAIR[1]

					# if both samples have data
					if (SPLINE[PAIR1 + 4] != "N") or ( SPLINE[PAIR2 + 4] != "N"):

						# add to the site count
						PAIRS[COUNTER1][2]  +=  int(BOOT[3])

						#add to the match count if they are matching
						if (SPLINE[PAIR1 + 4] == SPLINE[PAIR2 + 4] ):

							PAIRS[COUNTER1][3] += int(BOOT[3])

					#update the pair index
					COUNTER1 += 1

#list with the output statistic we care about
PAIRWISE_LIST = []

for PAIRWISE in PAIRS:

	PAIRWISE_LIST.append(PAIRWISE[3] / PAIRWISE[2])

#print all the uniq pairs
OUTPUT_LIST = []

#for each individual
for IND in range(0,IND_NUM):

	#create a new list in the output list
	OUTPUT_LIST.append([])

	#for this individual, pull out the correct number of pairwise comparisons
	OUTPUT_LIST[IND] = PAIRWISE_LIST[0:(IND_NUM - IND -1)]

	#remove those pairwise comparisons
	PAIRWISE_LIST = PAIRWISE_LIST[(IND_NUM - IND -1):]

	#add in self comparison for this individual
	OUTPUT_LIST[IND].insert(0,0)

	#if this isn't the first individual
	if IND > 0:

		#take the previous individual, reverse it, and iteratively add to current
		for PREVIOUS in OUTPUT_LIST[0:IND][::-1]:

			OUTPUT_LIST[IND].insert(0, PREVIOUS[IND])


for IND in OUTPUT_LIST:

	print IND
