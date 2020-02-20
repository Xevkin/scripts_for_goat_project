import sys

#tab seperated
regions_file = sys.argv[1]

chr = 0

SWITCH = 0

with open(regions_file) as regions:

	for line in regions:

		split_line = line.rstrip("\n").split("\t")
		current_start = split_line[1]
		current_end = split_line[2]
		if split_line[0].startswith("N"):

			current_chr = split_line[0]

			if (current_chr != chr):

				if SWITCH != 0:

					print str(chr) + "\t" + last_start + "\t" + last_end

				chr = current_chr
				last_start = split_line[1]
				last_end = split_line[2]
				SWITCH = 1
				continue

			if (int(last_end) <= int(current_start)):

				print str(chr) + "\t" + last_start + "\t" + last_end
				last_start = current_start
				last_end = current_end

	                else:

		                last_end =  current_end


		else:

			current_chr = int(split_line[0])

			if current_chr > chr:

				if not (chr == 0):

					print str(chr) + "\t" + last_start + "\t" + last_end
					chr = current_chr
					last_start = split_line[1]
					last_end = split_line[2]
					continue

			if (int(last_end) <= int(current_start)):

				print str(chr) + "\t" + last_start + "\t" + last_end

				last_start = current_start

				last_end = current_end

			else:

				last_end =  current_end




print str(chr) + "\t" + last_start + "\t" + last_end

