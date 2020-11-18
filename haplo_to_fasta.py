import sys

from subprocess import call

SAMPLES = []

TEMP_FILE = "temp.fa"

OUT_FILE = sys.argv[1].split(".")[0] +  ".fa"

call("rm " + TEMP_FILE,shell=True)

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		LINE = LINE.rstrip("\n")

		SPLINE = LINE.split("\t")

		if LINE.startswith("chr"):

			SAMPLE_NUM = len(SPLINE[3:])

			for CURRENT_SAMPLE in range(3, SAMPLE_NUM+3):

				call ("echo -e \"\\\\n>\""  + str(SPLINE[CURRENT_SAMPLE]) + " >> " + TEMP_FILE, shell=True)

				call("cut -f " + str(CURRENT_SAMPLE + 1) + " " + sys.argv[1] + " | tail -n +2 | tr \'\n\' \' \' | sed -e \"s/ //g\" >> " + TEMP_FILE,shell=True)

			break

call("tail -n +2 " + TEMP_FILE + " > tmp; mv tmp " + TEMP_FILE,shell=True)

call("awk -F \'\' \'{if ($0 ~ /^>/) print; else for(i=1;i<length;i+=70) print substr($0,i,70) }\' "  +  TEMP_FILE + " > " + OUT_FILE,shell=True)
