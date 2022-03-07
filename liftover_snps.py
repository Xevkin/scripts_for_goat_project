from __future__ import print_function

import sys

import pyfaidx

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def main():

	#give absolute path
	REF_TO_BE_LIFTED_TO=sys.argv[1]

	REF=pyfaidx.Fasta(REF_TO_BE_LIFTED_TO)

	#lifted snps should be a bim file
	LIFTED_SNP_FILE=sys.argv[2]

	#counters for the various possibilies
	MAJOR_REF_MATCH_COUNT = 0

	MINOR_REF_MATCH_COUNT = 0

	MAJOR_REF_FLIP_MATCH_COUNT = 0

	MINOR_REF_FLIP_MATCH_COUNT = 0

	FLIP_COUNT = 0

	FLIP_ISSUES = 0

	#for GC and TA variants
	PROBLEM_SNP_COUNT = 0

	with open(LIFTED_SNP_FILE) as FILE:

		for LINE in FILE:

			LINE=LINE.rstrip("\n")

			SPLINE = LINE.split("\t")

			ALLELE1=LINE.split("\t")[4].upper()

			ALLELE2=LINE.split("\t")[5].upper()

			#remove problem SNPS
			if (ALLELE1 == "G" and ALLELE2 == "C") or (ALLELE1 == "C" and ALLELE2 == "G") or (ALLELE1 == "T" and ALLELE2 == "A") or (ALLELE1 == "A" and ALLELE2 == "T"):

				PROBLEM_SNP_COUNT += 1

				#eprint(LINE + "\tGC don't print")

				continue

			CHR=LINE.split("\t")[0].replace("chr","")

			if CHR[0] in ["0","1","2","3","4","5","6","7","8","9"]  :

				CHR = int(CHR)

			#let's not use variants which have changed chromosome
			if str(CHR) != SPLINE[1].split(":")[0]:

				continue

			#position and in bim is base 1. This needs to accounted for by -1 when using pyfaidx
			POS=int(LINE.split("\t")[3])

			if isinstance(CHR,str):

				BASE_IN_REF=str(REF[CHR][POS-1:POS].seq).upper().replace("\n","").rstrip("\n")

			else:

				BASE_IN_REF=str(REF[CHR-1][POS-1:POS].seq).upper().replace("\n","").rstrip("\n")

			if ALLELE1 == BASE_IN_REF:

				MAJOR_REF_MATCH_COUNT += 1

				#print(LINE)

				print(SPLINE[1])

			elif ALLELE2 == BASE_IN_REF:

				MINOR_REF_MATCH_COUNT += 1

				print(SPLINE[1])

				#print("\t".join(SPLINE[0:4]) + "\t" + SPLINE[5] + "\t" + SPLINE[4])

			else:

				FLIP_COUNT += 1

				NEW_ALLELE_1 = flip_base(ALLELE1)

				NEW_ALLELE_2 = flip_base(ALLELE2)

				if NEW_ALLELE_1 == BASE_IN_REF:

					MAJOR_REF_FLIP_MATCH_COUNT += 1

					print(SPLINE[1])

					eprint(SPLINE[1])

					#print("\t".join(SPLINE[0:4]) + "\t" + NEW_ALLELE_1 + "\t" + NEW_ALLELE_2)

				elif NEW_ALLELE_2 == BASE_IN_REF:

					MINOR_REF_FLIP_MATCH_COUNT += 1

					print(SPLINE[1])

					eprint(SPLINE[1])

					#print("\t".join(SPLINE[0:4]) + "\t" + NEW_ALLELE_2 + "\t" + NEW_ALLELE_1 + "\tMINOR_FLIP")

				else:

					FLIP_ISSUES += 1

					#eprint(LINE + "\tFlip Issues")

	#print "Major allele matches reference " + str(MAJOR_REF_MATCH_COUNT) + " times."
	#print "Minor allele matches reference " + str(MINOR_REF_MATCH_COUNT) + " times."
	#print "Major allele flipped matches reference " + str(MAJOR_REF_FLIP_MATCH_COUNT) + " times."
	#print "Minor allele flipped matches reference " + str(MINOR_REF_FLIP_MATCH_COUNT) + " times."
	#print "Flipped alleles do not match " + str(FLIP_ISSUES) + " times."
	#print "GC|AT sites occur " + str(PROBLEM_SNP_COUNT) + " times."

def flip_base(BASE):

	BASE = BASE.upper()

	if BASE == "C":

		NEW_BASE = "G"

	elif BASE == "G":

		NEW_BASE = "C"

	elif BASE == "T":

		NEW_BASE = "A"

	else:

		NEW_BASE = "T"

	return(NEW_BASE)


main()
