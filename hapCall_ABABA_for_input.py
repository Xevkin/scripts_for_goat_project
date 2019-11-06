from __future__ import division
import sys

#input the sample numbers as in the bamlist

#expect "yes" or "no" for "transversions only"

#last input is the hapCall file

#unadmixed
H1_index = int(sys.argv[1]) + 2

#admixed?
H2_index = int(sys.argv[2]) + 2

#H3 will be the admixer
H3_index = int(sys.argv[3]) + 2

#confounder
H4_index = int(sys.argv[4]) + 2

outgroup_index = int(sys.argv[5]) + 2

transversions_only = sys.argv[6]

nABABA = 0

nBAABA = 0

nABBBA = 0

nBABBA = 0

nABBAA = 0

nBABAA = 0

nBBBAA = 0

nAAABA = 0

nBBABA = 0

nAABAA = 0

nABAAA = 0

nBAAAA = 0

nAAAAA = 0

nBBBBA = 0

with open(sys.argv[7]) as file:

	for line in file:

		if line.startswith("chr"):

			continue

		line = line.rstrip("\n")

		split_line = line.split("\t")

		observed_bases = [split_line[H1_index], split_line[H2_index], split_line[H3_index], split_line[H4_index], split_line[outgroup_index]]

		if "N" in observed_bases:

			continue

		observed_bases_unique = list( dict.fromkeys(observed_bases) )

		if (len(observed_bases_unique) != 2):

			continue

		if (  ("C" in  observed_bases_unique) and ("T" in  observed_bases_unique) ) or ( ("G" in  observed_bases_unique) and ("A" in  observed_bases_unique) ) and (transversions_only == "yes"):

			continue

		#extract and count ABABA sites
		if ( observed_bases[4] ==  observed_bases[0]) and ( observed_bases[4] ==  observed_bases[2]) and ( observed_bases[1] ==  observed_bases[3]) and (observed_bases[3] !=  observed_bases[4]):

			#print str(observed_bases) + " ABABA"

			nABABA += 1

			continue

		#BAABA sites
                if ( observed_bases[4] ==  observed_bases[1]) and ( observed_bases[4] ==  observed_bases[2]) and ( observed_bases[0] ==  observed_bases[3]) and (observed_bases[0] !=  observed_bases[4]):

			#print str(observed_bases) + " BAABA"

                        nBAABA += 1

                        continue

		#ABBBA sites
		if ( observed_bases[4] ==  observed_bases[0]) and ( observed_bases[1] ==  observed_bases[2]) and ( observed_bases[1] ==  observed_bases[3]) and (observed_bases[3] !=  observed_bases[4]):

			#print str(observed_bases) + " ABBBA"

			nABBBA += 1

			continue


	        #BABBA sites
                if ( observed_bases[4] ==  observed_bases[1]) and ( observed_bases[0] ==  observed_bases[2]) and ( observed_bases[0] ==  observed_bases[3]) and (observed_bases[4] !=  observed_bases[0]):

			#print str(observed_bases) + " BABBA"

                        nBABBA += 1

                        continue


		#ABBAA sites
		if ( observed_bases[4] ==  observed_bases[0]) and ( observed_bases[1] ==  observed_bases[2]) and  (observed_bases[3] ==  observed_bases[4]):

			#print str(observed_bases) + " ABBAA"

			nABBAA += 1

			continue

               #BABAA sites
                if ( observed_bases[4] ==  observed_bases[1]) and ( observed_bases[0] ==  observed_bases[2]) and ( observed_bases[1] ==  observed_bases[3]) and  (observed_bases[0] !=  observed_bases[4]):

			#print str(observed_bases) + " BABAA"

                        nBABAA += 1

                        continue


		#BBBAA sites
		if ( observed_bases[4] ==  observed_bases[3]) and ( observed_bases[0] ==  observed_bases[2]) and ( observed_bases[1] ==  observed_bases[0]) and  (observed_bases[0] !=  observed_bases[4]):

			nBBBAA += 1

			continue

		#nAAABA sites
		if ( observed_bases[4] ==  observed_bases[0]) and ( observed_bases[0] ==  observed_bases[1]) and ( observed_bases[1] ==  observed_bases[2]) and  (observed_bases[3] !=  observed_bases[4]):

			nAAABA += 1

			continue

		#BBABA sites
		if ( observed_bases[4] ==  observed_bases[2]) and ( observed_bases[0] ==  observed_bases[1]) and ( observed_bases[1] ==  observed_bases[3]) and  (observed_bases[0] !=  observed_bases[4]):

			nBBABA += 1

			continue

		#AABAA sites
		if ( observed_bases[4] ==  observed_bases[0]) and ( observed_bases[0] ==  observed_bases[1]) and ( observed_bases[1] ==  observed_bases[3]) and  (observed_bases[2] !=  observed_bases[4]):

			nAABAA += 1

			continue

		#ABAAA sites
		if ( observed_bases[4] ==  observed_bases[0]) and ( observed_bases[0] ==  observed_bases[2]) and ( observed_bases[2] ==  observed_bases[3]) and  (observed_bases[1] !=  observed_bases[4]):

			nABAAA += 1

			continue

                #BAAAA sites
                if ( observed_bases[4] ==  observed_bases[1]) and ( observed_bases[1] ==  observed_bases[2]) and ( observed_bases[2] ==  observed_bases[3]) and  (observed_bases[0] !=  observed_bases[4]):

                        nBAAAA += 1

                        continue

		if ( observed_bases[4] ==  observed_bases[1]) and ( observed_bases[1] ==  observed_bases[2]) and ( observed_bases[2] ==  observed_bases[3]) and  (observed_bases[0] ==  observed_bases[4]):

			nAAAAA += 1

		if ( observed_bases[4] !=  observed_bases[0]) and ( observed_bases[1] ==  observed_bases[2]) and ( observed_bases[2] ==  observed_bases[3]) and  (observed_bases[0] ==  observed_bases[1]):

			nBBBBA += 1

print "dummy\tdummy\t" + str(nAAAAA) + "\t" + str(nAAABA) + "\t"+ str(nAABAA) + "\t" + str(nAABBA) + "\t" + str(nABAAA) + "\t" + str(nABABA) + "\t" + str(nABBAA) + "\t" + str(nABBBA) + "\t" + str(nBAAAA) + "\t" + str(nBAABA) + "\t" + str(nBABAA) + "\t" + str(nBABBA) + "\t" + str(nBBAAA) + "\t" + str(nBBABA) + "\t" + str(nBBBAA) + "\t" + str(nBBBBA)

print "nABABA: " + str(nABABA)
print "nBAABA: " + str(nBAABA)
print "nABBBA: " + str(nABBBA)
print "nBABBA: " + str(nBABBA)
print "nABBAA: " + str(nABBAA)
print "nBABAA: " + str(nBABAA)
print "nBBBAA: " + str(nBBBAA)
print "nAAABA: " + str(nAAABA)
print "nBBABA: " + str(nBBABA)
print "nAABAA: " + str(nAABAA)
print "nABAAA: " + str(nABAAA)
print "nBAAAA: " + str(nBAAAA)
