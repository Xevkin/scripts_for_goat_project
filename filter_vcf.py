import sys

with open(sys.argv[1]) as file:

	for line in file:

		if not (line.startswith("#")):

			if not (line.startswith("\n")):
	
				split_line = line.rstrip("\n").split("\t")
				
				split_format = split_line[9].split(":")

				if (int(split_format[2]) >= 10) and (int(split_format[3]) >= 30):

					print line.rstrip("\n")  

		
