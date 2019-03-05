import sys

from random import randint

vcf=sys.argv[1]

with open(vcf) as file:

	for line in file:

		line=line.rstrip("\n")

		if line.startswith("#"):

			print line

		else:

			splitline=line.split("\t")

			newline=[]

			newline.extend(splitline[0:9])
			#account for cases where the format is GT:PL:DP:SP:ADF:ADR:AD:GP:GQ
			if (splitline[8].count(":") == 8):

				for sample in splitline[9:len(splitline)]:

					splitsample=sample.split(":")

					refreads=int(splitsample[6].split(",")[0])

					altreads=int(splitsample[6].split(",")[1])

					if (refreads == 0) and (altreads == 0):

						newline.append("./.:" + ":".join(splitsample[1:len(splitsample)]))

					elif sample.startswith("0/0") or sample.startswith("0/1") or sample.startswith("1/1"):

						total = refreads+altreads

						random_draw = randint(0,total-1)

						if random_draw >= refreads:

							newline.append("1/1:" + ":".join(splitsample[1:len(splitsample)]))

						else:

							newline.append("0/0:"+ ":".join(splitsample[1:len(splitsample)]))

				print "\t".join(newline)


			elif (splitline[8].count(":") == 5):
				#these only appear to be reference or missing, so can just print the line
				print line.rstrip("\n")

