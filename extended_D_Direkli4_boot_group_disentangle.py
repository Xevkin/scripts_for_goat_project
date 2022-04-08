#!/usr/bin/python
"""
Please see the docstring of extended_D_Direkli4_boot_group.py for an explanation on the inputs.

This script will print, for each current_genome in H4, the number of times a H3-derived allele was observed to be 

i) ancestral in the outgroup and H1
ii) derived in H2
iii) at >0%, <=10% in across all of H4 and
iv) derived in the current_genome (in H4)

All-variant and transversion only values are output.
"""
from __future__ import division

import sys

import gzip

import random

import numpy

#input is gzipped haplo file, list of H1, H2 group, H3 ind, list of H4s groups, grouped with _, the H5 outgroup groups, and if single mode
HAPLO_FILE = sys.argv[1]

H1_GROUP = sys.argv[2].split(",")

H1_INDICES = []

H2_GROUP = sys.argv[3].split(",")

H2_INDICES = []

H3 = sys.argv[4]

H4_GROUPS = [X.split(",") for X in sys.argv[5].split("_")]

H4_INDS = sys.argv[5].rstrip("\n").replace("_",",").split(",")

H4 = sys.argv[5].replace("_",",").split(",")

H5 = sys.argv[6].rstrip("\n").split(",")

H4_INDICES = []

H5_INDICES = []

H4_GROUP_INDICES = []

H4_DICT = {}

H4_TRANSV_DICT = {}

for H4_IND in H4_INDS:

	if H4_IND not in H4_DICT:

		H4_DICT[H4_IND] = 0

		H4_TRANSV_DICT[H4_IND] = 0

ALL_INDS = []

ALL_INDS.extend(H1_GROUP)

ALL_INDS.extend(H2_GROUP)

ALL_INDS.extend([H3])

ALL_INDS.extend(H5)

ALL_INDS.extend(H4)

SINGLE = sys.argv[7]

DER_OUTLIST = []

DER_OUTLIST_TRANSV = []


with gzip.open(HAPLO_FILE) as FILE:

	for LINE in FILE:

		LINE = LINE.rstrip("\n")

		SPLINE = LINE.split("\t")

		if LINE.startswith("chr"):

			for IND in ALL_INDS:

				if IND not in SPLINE:

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

				if (len(''.join(set(H5_BASES)).replace("N","")) == 1):

					ANC_ALLELE = ''.join(set(H5_BASES)).replace("N","")

				else:

					continue

				DER_ALLELE = SPLINE[H3_INDEX]

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

				if DER_FREQ > 0 and DER_FREQ <= 0.1:

						#this is so we can filter on no missing sites in our defined groups
						MISSING = "NO"

						for GROUP in H4_GROUP_INDICES:

							GROUP_BASES = [SPLINE[INDEX] for INDEX in GROUP]

							#if any of the groups are missing data, skip
							if GROUP_BASES.count("N") == len(GROUP):

									MISSING = "YES"

						if ( H4_BASES.count("N") != len(H4) ) and (MISSING == "NO"):

							# AB sites
							if ( H1_BASE == ANC_ALLELE) and (H1_BASE != H2_BASE) and (H2_BASE == DER_ALLELE) :

								for H4_IND in [H4_INDS[Y] for Y in  [i for i, x in enumerate(H4_BASES) if x==DER_ALLELE]]:

									H4_DICT[H4_IND] += 1

								if ["GA", "AG", "CT", "TC"].count("".join([ANC_ALLELE,DER_ALLELE])) == 0:

									for H4_IND in [H4_INDS[Y] for Y in  [i for i, x in enumerate(H4_BASES) if x==DER_ALLELE]]:

										H4_TRANSV_DICT[H4_IND] += 1


for H4_IND in H4_INDS:

	print H4_IND + " " +  str(H4_DICT[H4_IND]) + " " + str(H4_TRANSV_DICT[H4_IND])



