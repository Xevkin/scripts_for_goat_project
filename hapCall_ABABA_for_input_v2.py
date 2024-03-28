import sys
import gzip
import multiprocessing
import subprocess

#take the root of the haplo file names
INPUT_HAP_ROOT = sys.argv[1]

#should be numbers according to bamlist
#expect the first sample to be the outgroup ie sheep
ANC_POP_SET_INPUT = sys.argv[2].rstrip("\n").split(",")

#offset number to fit .haplo format
ANC_POP_SET_INPUT = [int(i)+2 for i in ANC_POP_SET_INPUT]

DER_POP_SET_INPUT = sys.argv[3].rstrip("\n").split(",")
DER_POP_SET_INPUT = [int(i)+2 for i in DER_POP_SET_INPUT]

#remove the output file if it exists already
subprocess.call("rm " + INPUT_HAP_ROOT + ".out ",shell=True)

#function for getting counts of derived allele sharing
def calc_X(INPUT_HAP_FILE_GZ, ANC_POP_SET=ANC_POP_SET_INPUT, DER_POP_SET=DER_POP_SET_INPUT):

	#create count variables
	COUNT = 0

	COUNT_TRANSV = 0

	#open the current file
	with gzip.open(INPUT_HAP_FILE_GZ) as FILE:

		#create a shallow copy of the ancestral allele indices
		POP_SET = ANC_POP_SET[:]

		#add the derived alleles
		POP_SET.extend(DER_POP_SET)

		#scan each line, skip if header
		for LINE in FILE:

			if LINE.startswith("chr"):

				continue

			LINE = LINE.rstrip("\n")

			SPLINE = LINE.split("\t")

			#store the ancestral alleles
			ANC_ALLELE = SPLINE[2]

			#store the current bases
			CURRENT_BASES = [SPLINE[i] for i in POP_SET]

			#if no missing data, store the bases for the ancestral and derived groups
			if "N" not in CURRENT_BASES:

				ANC_BASES = [SPLINE[i] for i in ANC_POP_SET]

				DER_BASES = [SPLINE[i] for i in DER_POP_SET]

				SPLINE_POP = [SPLINE[i] for i in POP_SET]

				#if the sheep sample matches the major ancestral allele
				#and all ancestral samples have the same allele
				#and all derived samples have the same allele
				#and derived-allele alleles are different
				if ANC_BASES[0] == ANC_ALLELE and \
				ANC_BASES.count(ANC_BASES[0]) == len(ANC_BASES) and \
				DER_BASES.count(DER_BASES[0]) == len(DER_BASES) and \
				ANC_BASES[0] != DER_BASES[0]:

					#add to the count
					COUNT += 1

					#if the variant is not a transition, add to the transversion count
					if not ( (ANC_BASES[0] == "G" and  DER_BASES[0] == "A") or (ANC_BASES[0] == "A" and  DER_BASES[0] == "G") or \
					(ANC_BASES[0] == "C" and  DER_BASES[0] == "T") or (ANC_BASES[0] == "T" and  DER_BASES[0] == "C") ):

						COUNT_TRANSV += 1

		#print the counts
		subprocess.call("echo " + INPUT_HAP_FILE_GZ + " " + str(COUNT) + " " + str(COUNT_TRANSV) + " >> " + INPUT_HAP_ROOT + ".out ",shell=True)

#input file variable
INPUT_FILES = []

#iterate through goat chromosomes
for CHR in range(1,30):

	INPUT_FILES.append(INPUT_HAP_ROOT + "_chr" + str(CHR) + ".haplo.gz")

#multithread over 8 threads
p = multiprocessing.Pool(8)

p.map(calc_X,INPUT_FILES)


#create final count variables
FINAL_COUNT = 0

FINAL_TRANSV_COUNT = 0

#run over the output file and total the final (all and transversion) counts
with open(INPUT_HAP_ROOT + ".out") as f:

	for line in f:

		split_f = f.rstrip("\n").split(" ")

		FINAL_COUNT += int(split_f[2])

		FINAL_TRANSV_COUNT += int(split_f[3])

print str(FINAL_COUNT) + " " + str(FINAL_TRANSV_COUNT)
