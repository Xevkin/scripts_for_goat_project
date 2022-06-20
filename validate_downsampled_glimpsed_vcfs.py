import gzip
import sys

#truth set should be vcf.gz
#sample set should be vcf

TRUTH_VARIANT_COUNT = 0

SAMPLE_VARIANT_COUNT = 0

#define transtions
TRANSITION_SET = [("T", "C"), ("G","A")]

#transitions then transverions
HET_MATCH_TRANS = [0, 0, 0, 0, 0, 0, 0]

HET_MATCH_TRNV = [0, 0, 0, 0, 0, 0, 0]

HET_MISMATCH_TRANS = [0, 0, 0, 0, 0, 0, 0]

HET_MISMATCH_TRNV = [0, 0, 0, 0, 0, 0, 0]

HOMO_ALT_MATCH_TRANS = [0, 0, 0, 0, 0, 0, 0]

HOMO_ALT_MATCH_TRNV = [0, 0, 0, 0, 0, 0, 0]

HOMO_REF_MATCH_TRANS = [0, 0, 0, 0, 0, 0, 0]

HOMO_REF_MATCH_TRNV = [0, 0, 0, 0, 0, 0, 0]

MISMATCH_TRANS = [0, 0, 0, 0, 0, 0, 0]

MISMATCH_TRNV = [0, 0, 0, 0, 0, 0, 0]

with gzip.open(sys.argv[1]) as TRUTH_GVCF:

	for TRUTH_LINE in TRUTH_GVCF:

		TRUTH_LINE = TRUTH_LINE.rstrip("\n")

		if TRUTH_LINE.startswith("#"):

			continue

		#add to variant count in truth data set
		TRUTH_VARIANT_COUNT += 1

		TRUTH_SPLINE = TRUTH_LINE.split()

		with open(sys.argv[2]) as SAMPLE_VCF:

			for SAMPLE_LINE in SAMPLE_VCF:

				SAMPLE_LINE = SAMPLE_LINE.rstrip("\n")

				if SAMPLE_LINE.startswith("#"):

					continue

				SAMPLE_SPLINE = SAMPLE_LINE.split()

				# skip to the matching variant and break if we go past
				if int(TRUTH_SPLINE[1]) < int(SAMPLE_SPLINE[1]):

					break

				if int(TRUTH_SPLINE[1]) > int(SAMPLE_SPLINE[1]):

					continue

				#check if missing data in sample
				if SAMPLE_SPLINE[9].startswith('./.'):

					break

				if TRUTH_SPLINE[1] == SAMPLE_SPLINE[1]:

					SAMPLE_VARIANT_COUNT += 1

					n = -1

					for MAF_LIMIT in [0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.10]:

						n += 1

						RAF=float(SAMPLE_SPLINE[7].split(";")[0].split("=")[1])

						if RAF > MAF_LIMIT:

							REF = SAMPLE_SPLINE[3]

							ALT = SAMPLE_SPLINE[4]

							SAMPLE_CALL = SAMPLE_SPLINE[9].split(":")[0]

							TRUTH_CALL = TRUTH_SPLINE[9].split(":")[0]

							if [SAMPLE_CALL,TRUTH_CALL] == ["0/0","0/1"] or [SAMPLE_CALL,TRUTH_CALL] == ["1/1","0/1"] :

								if (REF,ALT) in TRANSITION_SET:

									HET_MISMATCH_TRANS[n] += 1

								else:

									HET_MISMATCH_TRNV[n] += 1

							if [SAMPLE_CALL,TRUTH_CALL] == ["0/1","0/1"] :

								if (REF,ALT) in TRANSITION_SET:

									HET_MATCH_TRANS[n] += 1

								else:

									HET_MATCH_TRNV[n] += 1

							elif [SAMPLE_CALL,TRUTH_CALL] == ["1/1","1/1"]:

								if (REF,ALT) in TRANSITION_SET:

									HOMO_ALT_MATCH_TRANS[n] += 1

								else:

									HOMO_ALT_MATCH_TRNV[n] += 1

							elif [SAMPLE_CALL,TRUTH_CALL] == ["0/0","0/0"]:

								if (REF,ALT) in TRANSITION_SET:

									HOMO_REF_MATCH_TRANS[n] += 1

								else:

									HOMO_REF_MATCH_TRNV[n] += 1

							else:

								if (REF,ALT) in TRANSITION_SET:

									MISMATCH_TRANS[n] += 1

								else:

									MISMATCH_TRNV[n] += 1

					break

A= "truth_variant_count sample_variant_count"
B="homo_ref_match_trans_maf0 homo_ref_match_trans_maf0-1 homo_ref_match_trans_maf0-5 homo_ref_match_trans_maf1 homo_ref_match_trans_maf2-5 homo_ref_match_trans_maf5 homo_ref_match_trans_maf10"
C="homo_ref_match_trnv_maf0 homo_ref_match_trnv_maf0-1 homo_ref_match_trnv_maf0-5 homo_ref_match_trnv_maf1 homo_ref_match_trnv_maf2-5 homo_ref_match_trnv_maf5 homo_ref_match_trnv_maf10"
D="het_match_trans_maf0 het_match_trans_maf0-1 het_match_trans_maf0-5 het_match_trans_maf1 het_match_trans_maf2-5 het_match_trans_maf5 het_match_trans_maf10"
E="het_match_trnv_maf0 het_match_trnv_maf0-1 het_match_trnv_maf0-5 het_match_trnv_maf1 het_match_trnv_maf2-5 het_match_trnv_maf5 het_match_trnv_maf10"
F = B.replace("_ref_","_alt_")
G = C.replace("_ref_","_alt_")
H = D.replace("het_", "mis")
I = E.replace("het_", "mis")
J = D.replace("match","mismatch")
K = J.replace("trans","trvn")

print "downsampled " +  " ".join([A, B, C, D, E, F, G, H, I, J, K])

print sys.argv[2] + " " + " ".join([str(x) for x in [TRUTH_VARIANT_COUNT, SAMPLE_VARIANT_COUNT, " ".join([str(x) for x in HOMO_REF_MATCH_TRANS]), " ".join([str(x) for x in HOMO_REF_MATCH_TRNV]), " ".join([str(x) for x in HET_MATCH_TRANS]), " ".join([str(x) for x in HET_MATCH_TRNV]), " ".join([str(x) for x in HOMO_ALT_MATCH_TRANS]), " ".join([str(x) for x in HOMO_ALT_MATCH_TRNV]), " ".join([str(x) for x in MISMATCH_TRANS]) , " ".join([str(x) for x in MISMATCH_TRNV]), " ".join([str(x) for x in HET_MISMATCH_TRANS]), " ".join([str(x) for x in HET_MISMATCH_TRNV]) ] ] )
