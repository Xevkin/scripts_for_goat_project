#!/usr/bin/python2.7

import sys

#script.py <input file> <column to modify> <value to add> <column deliminator eg \t,\s

#adds a single value to every value of a particular column

#file must be tab seperated

#input should be the file with columns that will be printed with a single column modified


in_file = sys.argv[1]

column_to_add = int(sys.argv[2].rstrip("\n"))

to_add = float(sys.argv[3].rstrip("\n"))

with open(in_file) as file:

	for line in file:

		line_columns = line.split("\t")

		line_columns[column_to_add] = str(float(line_columns[column_to_add]) + to_add)

		new_line = "\t".join(line_columns)

		print new_line.rstrip("\n")
