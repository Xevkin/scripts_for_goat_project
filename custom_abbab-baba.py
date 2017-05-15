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

		if (line.split()[0:1] not in range(1,30) ):

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


print str(ABBA_count )

print str(BABA_count)
