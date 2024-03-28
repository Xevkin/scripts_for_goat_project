##converts two input files to the individual Q-matrix file
#first input is the working file - the Q file with the first col being population ID
#second input is a "individual key population" fil, tab deliminated, which links individual IDs to both population
#and the population key

import sys

#individual_Q_file
ind_Qfile = sys.argv[1]

indiv_key_pop = sys.argv[2]

indiv_key_list = []

with open(indiv_key_pop) as file:

	for line in file:

		indiv_key_list.append(line.rstrip("\n"))


with open(ind_Qfile) as file:
	#cycle through the pop_ind_key file
	#extract data
	#assuming they are ordered the same as in the working file, print the indivdual Q file in the correct format
	i=0
	for line in file:
		
		split_line = line.split(" ")
		
		cur_indiv = indiv_key_list[i].split("\t")
		
		population = cur_indiv[0]
		individual_ID = cur_indiv[1]
		pop_key = cur_indiv[2] 
		
		#col 2 is code/ID for individual
		#col 4 is population code (e.g. 47) for that individual
		#cols 6 onwards show membership coefficient
 
		start_cols = [pop_key, individual_ID, pop_key, pop_key, pop_key]

		new_line = start_cols + split_line
		
		print ' '.join(new_line).rstrip("\n") 
		

		i = i+1
