import sys

chromosomes = []

BLOCK = 5000000

with open("/home/kdaly/chromosomes.length") as file:

	for line in file:

		spline = line.rstrip("\n").split("\t")

		chromosomes.append(spline)


for chromosome in chromosomes:

	CHR_END = int(chromosome[1])

	CHR = chromosome[0]

	START = 1

	END = START + BLOCK - 1

	while END < CHR_END:

		print CHR + ":" + str(START) + "-" + str(END)

		START = START + BLOCK

		END = END + BLOCK

	print CHR + ":" + str(CHR_END - BLOCK) + "-" + str(CHR_END)
