import os

import sys

import subprocess

import numpy

import random

random_number = str(random.random())

seed = random_number.split(".")[1]

#remove the previous parallel file if it had the same seed
subprocess.call("rm parallel_" + seed + ".txt",shell=True)

bootstrap_size = 5000000

chromosomes = []

with open("/home/kdaly/chromosomes.length") as file1:

	for line in file1:

		split_line = line.rstrip("\n").split("\t")

		to_add = [split_line[0], split_line[1]]

		chromosomes.append(to_add)

#create an output file
output_file = "tmp-"+ seed +  ".out"

subprocess.call("rm " + output_file,shell=True)

#for each chromosome
for chr in chromosomes:

	#set the length of the chromosome
	chr_len_remain = int(chr[1])

	#set a start position that will move with each iteration
	start_pos = 1

	#so we can parallelize
	i=1

	#this is the bootstrap step

	#keep going until we have no block big enough
	while (chr_len_remain >= bootstrap_size):

		#initialize some files for the parallel step
		current_file = "tmp-"+ seed + "-chr" + chr[0] + "-" + "boot" + str(i) + ".txt"

		f=open(current_file,"w+")

		current_limit = start_pos + bootstrap_size

		#take the input file
		with open("test.haplo") as file2:

			for line2 in file2:

				split_line2 = line2.rstrip("\n").split("\t")

				#now check if the current line has the same chromosome and falls within the current bootstrap range
				#if it does, do not print the line
				if not ( (split_line2[0] == chr[0]) and (int(split_line2[1]) >= start_pos) and (int(split_line2[1]) <= current_limit) ):

					f.write(line2)


		#update the start position
		start_pos += bootstrap_size

		#decrease the amount of chromosome remaining
		chr_len_remain -= bootstrap_size

		f.close()

		#update the counter
		#files which have no data in the bootstrap chunk will not be included in the parallel command
                i += 1

		#if file we've made has the same number of lines i.e. their was no data in the chunk remove, skip this chunk
		if os.stat(sys.argv[1]).st_size != os.stat(current_file).st_size:

			#now, run the python script using the bootstrap file

			subprocess.call("echo python /home/kdaly/programs/scripts_for_goat_project/biallelic_lineage_specific.py " + sys.argv[2] + " " + sys.argv[3] + " " + sys.argv[4] + "  " + current_file + " \"|\"  wc -l \">>\" " + output_file + " >> parallel_" + seed + ".txt",shell=True)


subprocess.call("parallel -a parallel_" + seed + ".txt -j 24",shell=True)

subprocess.call("rm tmp-"+ seed + "*.txt",shell=True)


#intialize a list for the bootstrap estimates
bootstraps = []

#now gather the output and print the variance
with open(output_file) as file:

	for line in file:

		bootstraps.append(float(line.rstrip("\n")))

	print "The number of entries is: " + str(len(bootstraps))

	variance = numpy.var(bootstraps)

	print "The variance is: " +  str(variance)
