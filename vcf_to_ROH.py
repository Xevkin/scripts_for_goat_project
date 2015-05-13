#takes input vcf file, and returns output in the format required for HHn
#Chromosome ObsHeterozPos_bp NextHeterozPos_bp LengthHomozRoh_bp
#I will assume \s seperation for output

from os import sys

def main(input_file):
	"""Required to produce input for HHn. vcf>ROH file"""

	with open(input_file) as vcf_file:

		#initialize some variables
		prev_pos = "0"
		prev_chrom = "0"

		#initialize a list that will take individual entries
		ROH_list = [["Chromosome" ,"ObsHeterozPos_bp", "NextHeterozPos_bp","LengthHomozRoh_bp"]]

		for line in vcf_file:
	
			#skip header lines
			if not line.startswith(r"#"):

				
				split_line = line.split("\t")
				
				#take chromosome and position information
				chrom = split_line[0] 

				NextHetPos = split_line[1]
				
				#if the script has moved onto a new chromosome, set the previous position back to "0"
				if (int(prev_chrom) < int(chrom)):

					ObsHetPos = "0" 
				
				#calculate the length of the ROH
				ROH_length = str(int(NextHetPos) - int(ObsHetPos))
		
				#add entry to ROH_list if it ISN'T a stretch starting from 0
				if (ObsHetPos != "0"):
	
					ROH_list.append([chrom, ObsHetPos, NextHetPos, ROH_length])
				
				ObsHetPos = NextHetPos

				prev_chrom = chrom


		return(ROH_list)

#run input vcf into the main function
#save output of main function in a list
input_vcf = sys.argv[1].rstrip("\n") 

output_list = main(input_vcf)

#print out the output in the correct format
for list in output_list:

	print ' '.join(list)
