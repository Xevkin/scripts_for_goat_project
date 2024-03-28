import gzip

import sys

input_file = sys.argv[1]
H1 = int(sys.argv[2]) + 2
H2 = int(sys.argv[3]) + 2
H3 = int(sys.argv[4]) + 2

ABBA_count = 0

BABA_count = 0

with gzip.open(input_file) as file:

	for line in file:

		if (line.startswith("chr")):

			continue

		if not line.startswith("1") and not line.startswith("2") and not line.startswith("3") and not line.startswith("4") and not line.startswith("5") and not line.startswith("6") and not line.startswith("7") and not line.startswith("8") and not line.startswith("9"):

			break

		#remove Ns - FOR COMPARISON
		if any("N" in s for s in line.split()[2:3]+line.split()[H1:H1+1]+line.split()[H2:H2+1]+line.split()[H3:H3+1] ):

			continue

		else:

			#remove transitions
			if (line.split()[2:3] == "C" and line.split()[H3:H3+1] == "T")  or (line.split()[2:3] == "G" and line.split()[H3:H3+1] == "A"):

				continue

		if ((line.split()[2:3] ==  line.split()[H1:H1+1]) and  (line.split()[H2:H2+1] ==  line.split()[H3:H3+1]) and (line.split()[2:3] !=  line.split()[H2:H2+1])):

			ABBA_count += 1


		elif ((line.split()[2:3] ==  line.split()[H2:H2+1]) and  (line.split()[H1:H1+1] ==  line.split()[H3:H3+1]) and (line.split()[2:3] !=  line.split()[H1:H1+1])):
				
			BABA_count += 1

		else:

			continue


print str((ABBA_count - BABA_count)/ (ABBA_count + BABA_count))


