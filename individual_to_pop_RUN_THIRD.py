#Inputs should be indiv Q file, number of pops, pop_ind_key file, and the name of output file

from __future__ import division

from operator import add

import sys

import csv

#Q file
Q_file = sys.argv[1]

master_list = []

number_of_pops = int(sys.argv[2])

#get data from the pop_ind_key file so Q file can be labelled by population code
pop_ind_key = sys.argv[3]

#save this in a list
indiv_key_list = []

with open(pop_ind_key) as file:

        for line in file:

                indiv_key_list.append(line.rstrip("\n"))



with open(Q_file) as Q:
	
	#for each population in the study
	#i variable tracks the population code
	for i in range(1,(number_of_pops+1)):
				
		#create a list to put individuals of a given pop into
		pop = []
	
		#tracking variable so script knows which individual is being processed
		m = 0 
		
		#for each individual in the Q file
		for line in Q:
						
			split_line = line.rstrip("\n").split(" ")
			
			#get the information on that individual, turn it into a list
			info = indiv_key_list[m]
			split_info = info.split("\t")

			#get population identifier of that individual
			pop_of_ind = split_info[2]
			
			#inserts element into an index, pushing later indexes  
			split_line.insert(0, pop_of_ind)
			
			#add the pop ID to the starting position of the line
			line = " ".join(split_line)
			
			#now check if the population we're interested contains this individual
			#if it does, add them to the pop list
			if line.startswith(str(i) + " "):

				pop.append(line)  
			#move on to the next individual
			m = m + 1
		
		#go back to the start of the Q file so each individual can be looped through again
		Q.seek(0)
		
		#Make a variable for the ancestry components
		#we add all the proportions up and then average them against the number of individuals in the pop
		#variable starts with [0] as we will dump the pop ID into it
		proc_pop = [0]

		#we change the size of this variable depending on K
		for b in range(1, len(split_line)):

			proc_pop.append(0)
		
		for indiv in pop:

			split_indiv = indiv.split(" ")
			
			#converts all strings in this list to floats
			split_indiv = map(float, split_indiv)

			#add the ancestry components to the intialized variable
 			proc_pop = map(add, proc_pop, split_indiv)

		pop_count = len(pop)
		
		#for each totalled ancestry component
		#divide that total by the number of individuals in that population
		averaged = ([ round(x / pop_count, 3) for x in proc_pop])

		
		averaged[0] = str(int(split_indiv[0])) + ":"
		
		averaged.append(pop_count)
		master_list.append(averaged) 
	
#need to write to file				
				
	Q.seek(0)	
	
with open(sys.argv[4].rstrip("\n"), "w") as f:

	writer = csv.writer(f, delimiter=' ')
	writer.writerows(master_list)
