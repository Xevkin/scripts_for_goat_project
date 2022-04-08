#!/usr/bin/python

from __future__ import division

import sys

import gzip

import random

import numpy

#input is gzipped haplo file, list of H1, H2 group, H3 ind, list of H4s groups, grouped with _, the H5 outgroup

#create the bootstrap regions
BOOTS=[]

with open("boots.list") as FILE:

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

H5 = sys.argv[6].rstrip("\n")

H4_INDICES = []

H4_GROUP_INDICES = []

ALL_INDS = H1_GROUP

ALL_INDS.extend(H2_GROUP)

ALL_INDS = [H3, H5]

ALL_INDS.extend(H4)

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

			H5_INDEX = SPLINE.index(H5)

		#now start the actual ABAAAA BAAAA etc calculations
		else:

			if  ([SPLINE[INDEX] for INDEX in [H3_INDEX, H5_INDEX ]].count("N") == 0) and ([SPLINE[INDEX] for INDEX in H1_INDICES ].count("N") != len(H1_INDICES)) and ([SPLINE[INDEX] for INDEX in H2_INDICES ].count("N") != len(H2_INDICES)):

				POSITIONS = [int(X) for X in SPLINE[0:2]]

				#filter on removing derived alleles present in our H4s
				H4_BASES = [SPLINE[INDEX] for INDEX in H4_INDICES]

				H2_BASES = [SPLINE[INDEX] for INDEX in H2_INDICES]

				ANC_ALLELE = SPLINE[H5_INDEX]

				DER_ALLELE = SPLINE[H3_INDEX]

				if ANC_ALLELE == DER_ALLELE:

					continue

				H1_BASES = [SPLINE[INDEX] for INDEX in H1_INDICES]

				#randomly sample a H1 base
				H1_BASE = "N"

				while H1_BASE == "N":

					H1_BASE = H1_BASES[random.randint(0,len(H1_BASES)-1)]

				H2_BASE = "N"

				while H2_BASE == "N":

					 H2_BASE = H2_BASES[random.randint(0,len(H2_BASES)-1)]

				#get the derived allele frequency of present bases
				DER_FREQ = H4_BASES.count(DER_ALLELE) /  len(list(filter(lambda a: a != "N", H4) ) )

				BOOT_COUNTER = -1

				#for each bootstrap region
				for BOOT in BOOTS:

					BOOT_COUNTER += 1

					BOOT_SPLIT = BOOT

					if not ( (BOOT[0] == POSITIONS[0]) and ( POSITIONS[1] >= BOOT[1]) and (POSITIONS[1] <= BOOT[2]) ):

						continue

					else:
						#this is so we can filter on no missing sites in our defined groups
						MISSING = "NO"

						for GROUP in H4_GROUP_INDICES:

							GROUP_BASES = [SPLINE[INDEX] for INDEX in GROUP]

							#if any of the groups are missing data, skip
							if GROUP_BASES.count("N") == len(GROUP):

								MISSING = "YES"

						if ( H4_BASES.count("N") != len(H4) ) and ( DER_FREQ <= 0.1 ) and (MISSING == "NO") and (DER_FREQ>= 0):

							COUNTER = 0

							DERIVED_H4 = []
							# AB sites
							for H4_BASE in H4_BASES:

								COUNTER += 1

								if H4_BASE == DER_ALLELE:

									DERIVED_H4.append(H4[COUNTER])

							print " ".join(SPLINE[0:2] + [ANC_ALLELE] + [DER_ALLELE] + H1_GROUP + H2_GROUP + DERIVED_H4)

