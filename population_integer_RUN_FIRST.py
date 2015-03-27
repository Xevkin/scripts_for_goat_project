import sys
#input files are:
#1. a bed file containing breed/population in column 1 [0], individual ID in column 2 [1]
#2. a file containing each breed/population - can just "cut -f1 -d' ' | sort | uniq"
#output file will be the "Population	Individual	Pop_Key" file

ID_list = sys.argv[1]

uniq_pops = sys.argv[2]

pops = []

#make a list of all pops
with open(uniq_pops) as file:
	for line in file:
		pops.append(line.rstrip("\n"))

#go through file that has both individual ID and breed
#once the breed matches, print appropriate number

with open(ID_list) as q:
	
	for individual in q:
			
		split_individual = individual.split(" ")
		
		i = 0
		a = 0
		
		while (i != 255):
			 
			if (split_individual[0].rstrip("\n") == pops[i]):

				i = 255

				print split_individual[0].rstrip("\n") + "\t" + split_individual[1].rstrip("\n") + "\t"+ str(a + 1)

			else:

				i = i + 1

				a = a + 1
