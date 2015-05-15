#takes input vcf file, and returns output in the format required for HHn
#Chromosome ObsHeterozPos_bp NextHeterozPos_bp LengthHomozRoh_bp
#I will assume \s seperation for output
#Going to ignore the X and Y chromosomes,as X het calls are errors in males, only relevant in females

from os import sys

def main(input_file):
	"""Required to produce input for HHn. vcf>ROH file"""

	with open(input_file) as vcf_file:

		#initialize some variables
		ObsHetPos = "0"
		prev_chrom = "0"
		
		#initialize a list that will take individual entries
		ROH_list = [["Chromosome" ,"ObsHeterozPos_bp", "NextHeterozPos_bp","LengthHomozRoh_bp"]]

		for line in vcf_file:
			
			#skip header lines
			if not (line.startswith(r"#") or line.startswith("chrX") or line.startswith("chrY")):

				
				split_line = line.split("\t")
				
				#skip this variant if it is not a SNP
				if  (split_line[4] not in ["A", "G", "T", "C"]):

					continue

				#skip variant if it is not a heterozygous position
				if (split_line[9].split(":")[0] not in ["1/0", "0/1", "1|0", "0|1"]):

					continue
				#take chromosome and position information
				chrom = split_line[0] 

				NextHetPos = split_line[1]
								
				if (int(ObsHetPos) > int(NextHetPos)):

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
