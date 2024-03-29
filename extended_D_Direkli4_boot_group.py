#!/usr/bin/python
"""
Update 3/7/22 : adding 7th variable to turn on outgroup ascertainment



I apologize for the current state of this script, it was not envisioned as being shared but in hindsight that was very unwise. It is very much not-optimized, and can be slow to run given the
size of input files

The script is currently hardcoded for my own study organism and current reference genome i.e. goat, ARS1. I will be working on making it
compliant with user-defined files and improving comments in the coming weeks - please contact me @ the email below if I have not down this or if you have questions.

This script takes a haplo.gz file as its main input, as output by ANGSD -doHaploCall.

The next five inputs define the H1, H2, H3, H4, and "true" outgroup groups/genomes.

Within each group, sample names (as in the haplo.gz header) are comma-seperated.

Within H4 I allow multiple groups to be defined; sites must be covered at least once by each H4 group and also be fixed for ancestral before a H3 derived allele is examined

H4 groups are delimited by the underscore _ character. Within each H4 group, genomes are comma seperated.

A sixth variable is required to turn on "single" mode. All analyses reported in the Direkli paper did not use this option, so should be set to "no"

A seventh variable turns on outgroup (H5) ascertainment ("yes" or "Yes" to turn on)

Two example inputs are provided below:


python extended_D_Direkli4_boot_group.py H1-Genome1,H1-Genome2 H2-Genome H3-Genome H4-Group1-Genome1,H4-Group1-Genome2_H4-Group2-Genome1_H4-Group3-Genome1 Outgroup-Genome1,Outgroup-Genome2 no no


python extended_D_Direkli4_boot_group.py Ganjdareh3,Ganjdareh22 Blagotin3 Direkli4 Direkli1-2,Direkli5_Ibex1,Falconeri1,Sibirica1,Tur1,Tur2_IranBezoar1,ZagrosBezoar1 Sheep1,Sheep2,Sheep3 no no


Typically I will run a parallelized loop across H2-Genomes.

The regions used for bootstrap analyses are currently hardcoded into the script. The bootstrap region file is space seperated: chr block_start block_end

Note that the 11 "cutoff" thresholds refer to the frequency of the H3 derived allele among all H4 groups: <= 0%, <= 10%, <= 20%.
The calculated statistics for each of these cutoffs are outputted from the script. If only interested in those at <=0%, examine the first row of the output.

The output columns are the extended D value, standard error, Z score, number of ABBAA sites,  number of nBABAA sites, and the same statistics but transversions only.

Kevin G. Daly, dalyk1@tcd.ie
"""
from __future__ import division

import sys

import gzip

import random

import numpy

#input is gzipped haplo file, list of H1, H2 group, H3 ind, list of H4s groups, grouped with _, the H5 outgroup groups, and if single mode

#create the bootstrap regions
BOOTS=[]

with open("/raid_md0/kdaly/papers/direkli4_paper2021/analyses/extended_D/boots.list") as FILE:

	for LINE in FILE:

		SPLINE = LINE.rstrip("\n").split(" ")

		BOOTS.append([int(X) for X in SPLINE] + [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ])

HAPLO_FILE = sys.argv[1]

H1_GROUP = sys.argv[2].split(",")

H1_INDICES = []

H2_GROUP = sys.argv[3].split(",")

H2_INDICES = []

H3 = sys.argv[4]

H4_GROUPS = [X.split(",") for X in sys.argv[5].split("_")]

H4 = sys.argv[5].replace("_",",").split(",")

H5 = sys.argv[6].rstrip("\n").split(",")

H4_INDICES = []

H5_INDICES = []

H4_GROUP_INDICES = []

ALL_INDS = []

ALL_INDS.extend(H1_GROUP)

ALL_INDS.extend(H2_GROUP)

ALL_INDS.extend([H3])

ALL_INDS.extend(H5)

ALL_INDS.extend(H4)

SINGLE = sys.argv[7]

OUTGROUP_ASCERT = sys.argv[8]

SINGLE_TOTAL_COUNT = 0

SINGLE_B_COUNT = 0

SINGLE_B_COUNT_TRANSV = 0

with gzip.open(HAPLO_FILE) as FILE:

	for LINE in FILE:

		LINE = LINE.rstrip("\n")

		SPLINE = LINE.split("\t")

		if LINE.startswith("chr"):

			for IND in ALL_INDS:

				if IND not in SPLINE:

					print SPLINE

					print IND + " is not in the header. Exiting..."

					raise SystemExit

			for H1_INDIV in H1_GROUP:

				H1_INDICES.append(SPLINE.index(H1_INDIV))

			for H4_INDIV in H4:

				H4_INDICES.append(SPLINE.index(H4_INDIV))

			for H2_INDIV in H2_GROUP:

				H2_INDICES.append(SPLINE.index(H2_INDIV))

			for GROUP in H4_GROUPS:

				H4_GROUP_INDICES.append([SPLINE.index(X) for X in GROUP])

			H3_INDEX = SPLINE.index(H3)

			for H5_INDIV in H5:

				H5_INDICES.append(SPLINE.index(H5_INDIV))

		#now start the actual ABAAAA BAAAA etc calculations
		else:

			if  ([SPLINE[INDEX] for INDEX in [H3_INDEX]].count("N") == 0)  and ([SPLINE[INDEX] for INDEX in H5_INDICES ].count("N") != len(H5_INDICES))  and ([SPLINE[INDEX] for INDEX in H2_INDICES ].count("N") != len(H2_INDICES)):

				if SINGLE not in ["yes", "Yes", "Single"]:

					if not ([SPLINE[INDEX] for INDEX in H1_INDICES ].count("N") != len(H1_INDICES)):

						continue

				POSITIONS = [int(X) for X in SPLINE[0:2]]

				#filter on removing derived alleles present in our H4s
				H4_BASES = [SPLINE[INDEX] for INDEX in H4_INDICES]

				H5_BASES = [SPLINE[INDEX] for INDEX in H5_INDICES]

				H2_BASES = [SPLINE[INDEX] for INDEX in H2_INDICES]

				DER_ALLELE = SPLINE[H3_INDEX]

				if (len(''.join(set(H5_BASES)).replace("N","")) == 1):

					ANC_ALLELE = ''.join(set(H5_BASES)).replace("N","")

				else:

					if OUTGROUP_ASCERT not in ["yes", "Yes"]:

						continue

					else:

						ANC_ALLELE = ''.join(set(H5_BASES)).replace("N","").replace(DER_ALLELE,"")

						if (DER_ALLELE not in H5_BASES):

							continue

				if ANC_ALLELE == DER_ALLELE:

					continue

				H1_BASES = [SPLINE[INDEX] for INDEX in H1_INDICES]

				#randomly sample a H1 base
				H1_BASE = "N"

				if SINGLE not in ["yes", "Yes", "Single", "single"]:

					while H1_BASE == "N":

						H1_BASE = H1_BASES[random.randint(0,len(H1_BASES)-1)]

				H2_BASE = "N"

				while H2_BASE == "N":

					 H2_BASE = H2_BASES[random.randint(0,len(H2_BASES)-1)]

				#get the derived allele frequency of present bases
				DER_FREQ = H4_BASES.count(DER_ALLELE) /  len(list(filter(lambda a: a != "N", H4) ) )

				if SINGLE in ["single", "Yes", "yes", "Single"]:

					MISSING = "NO"

					for GROUP in H4_GROUP_INDICES:

						GROUP_BASES = [SPLINE[INDEX] for INDEX in GROUP]

						if GROUP_BASES.count("N") == len(GROUP):

							MISSING = "YES"

					if ( H4_BASES.count("N") != len(H4) )  and (MISSING == "NO"):

						SINGLE_TOTAL_COUNT += 1

						if  (H2_BASE == DER_ALLELE) and (DER_FREQ == 0):

							#print [H2_BASE, DER_ALLELE, H4_BASES, ANC_ALLELE]

							SINGLE_B_COUNT += 1

							if  ["GA", "AG", "CT", "TC"].count("".join([ANC_ALLELE,DER_ALLELE])) == 0:

								SINGLE_B_COUNT_TRANSV += 1

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

							for GROUP in H4_GROUP_INDICES:

								GROUP_BASES = [SPLINE[INDEX] for INDEX in GROUP]

								#if any of the groups are missing data, skip
								if GROUP_BASES.count("N") == len(GROUP):

									MISSING = "YES"

							if ( H4_BASES.count("N") != len(H4) ) and ( DER_FREQ <= DER_CUTOFF ) and (MISSING == "NO"):

								# AB sites
								if ( H1_BASE == ANC_ALLELE) and (H1_BASE != H2_BASE) and (H2_BASE == DER_ALLELE) :

									BOOTS[BOOT_COUNTER][3][COUNTER] += 1

									if ["GA", "AG", "CT", "TC"].count("".join([ANC_ALLELE,DER_ALLELE])) == 0:

										 BOOTS[BOOT_COUNTER][4][COUNTER]  += 1

								if ( H2_BASE == ANC_ALLELE) and (H1_BASE != H2_BASE) and (H1_BASE == DER_ALLELE):

									BOOTS[BOOT_COUNTER][5][COUNTER]  += 1

									if ["GA", "AG", "CT", "TC"].count("".join([ANC_ALLELE,DER_ALLELE])) == 0:

										 BOOTS[BOOT_COUNTER][6][COUNTER]  += 1

if SINGLE in ["yes","single", "Single"]:

	print " _".join(H2_GROUP) + " " +  " ".join([str(X) for X in [SINGLE_B_COUNT, SINGLE_B_COUNT_TRANSV, SINGLE_TOTAL_COUNT]])

	exit()

#list of the 11 different cutoffs and extended D estimates
OUTPUT = [[] , [], [], [], [] , [], [], [], [] , [], []]

TRANSV_OUTPUT = OUTPUT

#calculate the final extended D too
DE_OUT = [ [0, 0, 0, 0] , [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0] , [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0] , [0, 0, 0, 0],[0, 0, 0, 0]]

DE_FINAL = []

#for each bootstrap region
for BLOCK in BOOTS:

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
	for SUBSAMPLE in range(0, len(BOOTS)):

		#randomly sample with replacement
		SUBSAMPLE_INDEX = random.randrange(0, len(BOOTS))

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

			EXTENDED_D_TRANSV = "NA"

			OUTPUT[CUTOFF].append(EXTENDED_D)

			OUTPUT[CUTOFF].append(EXTENDED_D_TRANSV)

		else:

			EXTENDED_D = (AB_TOTAL[CUTOFF] - BA_TOTAL[CUTOFF]) / (AB_TOTAL[CUTOFF] + BA_TOTAL[CUTOFF])

			EXTENDED_D_TRANSV = (AB_TRANSV_TOTAL[CUTOFF] - BA_TRANSV_TOTAL[CUTOFF]) / (AB_TRANSV_TOTAL[CUTOFF] + BA_TRANSV_TOTAL[CUTOFF])

			OUTPUT[CUTOFF].append(EXTENDED_D)

			OUTPUT[CUTOFF].append(EXTENDED_D_TRANSV)

TO_PRINT = [sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]]

for CUTOFF in range(0, 11):

	if "NA" in OUTPUT[CUTOFF]:

		TO_PRINT.extend(["NA_" + str(CUTOFF)])

		continue

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

print " ".join(TO_PRINT)
