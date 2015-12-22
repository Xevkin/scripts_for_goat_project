#!/usr/bin/env python
'''
Script is run in a directory with bam files to be aligned to either goat, wild goat or sheep genome
Reads will be trimmed, aligned, read groups added, duplicates removed, then 

Also need to supply a date for the Hiseq run - will automatically make results file
python script <date_of_hiseq> <meyer> <species> <trim> <align> <read group file> <directory in which to place output dir/>

read group file needs to be in the follow format (\t are actual tab characters)
FASTQ_FILE.GZ\t@RG\tID:X\tSM:X\tPL:X\tLB:X\tLANE\tSAMPLE_NAME
ID should be in the format <sample_name>-<macrogen_index_number>-<lane_number>-<hiseq_number>
LB should refer to PCR: <sample>-<lab index>-<macrogen-index>-<PCR_number>



'''

#import csv module to easily write the list of list used as a .csv file
import csv

#apparently this is a better way to make system calls using python, rather than "os.system"
import subprocess 
from subprocess import call

#need to import sys anyways to access input file (a list)
import sys

#we will use this later to check if the input files actually exist
import os

#dictionary of species names and genome paths
nuclear_genomes = {

	"goat" : "/kendrick/reference_genomes/goat_CHIR1_0/goat_CHIR1_0.fasta",

	"sheep" : "/kendrick/reference_genomes/sheep_oviAri3/oviAri3.fa",

	"bezoar" : "/kendrick/reference_genomes/bezoar_CapAeg_1_0/CapAeg_renamed.fa"
}


def main(date_of_hiseq, meyer, species, trim, align, process, RG_file, output_dir):
	
	#run the set up function.#set up will create some output directories
	#and return variables that will be used in the rest of the script
	
	files, reference, out_dir, cut_adapt, alignment_option, fastq_list = set_up(date_of_hiseq, meyer, species, RG_file, output_dir, trim) 
	
	#sample is the file root
	#trim fastq files and produce fastqc files - if we want to
	#the masterlist will change each time so it needs be equated to the function
	
	if (trim == "yes"):
 
		for fastq in fastq_list:
		
			trim_fastq(fastq, cut_adapt, out_dir)
	
	#at this stage we have our fastq files with adaptors trimmed
	#We can now move on to the next step: alignment
	
	#going to align to CHIR1.0, as that what was used for AdaptMap
	#this step will align each fastq and produce raw bams with RGs
	#they should have the same stem of the initial bam


	if (align == "yes"):
 
		map(lambda fastq : align(fastq, RG_file, alignment_option, reference, trim), fastq_list)

	#Remove the sai files that we don't care about
	call("rm *sai",shell=True)

	#sort and remove duplicates from each bam	
	if (process == "yes"):

		map(process_bam, fastq_list)

        #now merge the lanes for each sample, and rmdup the merged bams
	merged_bam_list = merge_lanes_and_sample(RG_file)
	
	realigned_bam_list = []

	#indel realignment
	for i in merged_bam_list:

		print "Indel realignment on sample " + i

		realigned_bam_list.append(indel_realignment(i,reference))

	#now run the script to process the realigned bams
	#this function will invoke a function that rescales the bams

	for i in realigned_bam_list:

		process_realigned_bams(i,reference,output_dir)
	
		
	#clean up files
	
	call("mkdir trimmed_fastq_files_and_logs",shell=True)
	call("mv *trimmed* trimmed_fastq_files_and_logs/",shell=True)
	
	call("gzip *",shell=True)
	
	for sample in fastq_list:
			
		#going to make an output directory for each sample
		#then move all produced files to this directory
		call("mkdir " + out_dir + sample,shell=True)
		print "mv *" + sample + "*.bam* "+ sample + "*.idx* "+ sample + "*flagstat* " + out_dir + sample
		call("mv *" + sample + "*.bam* "+ sample + "*.idx* "+ sample + "*flagstat* " + out_dir + sample,shell=True)

	
	call("gzip trimmed_fastq_files_and_logs/*",shell=True)
	call("mv trimmed_fastq_files_and_logs/ " + out_dir,shell=True)

	output_summary = date_of_hiseq + "_summary.table"
	
	call("mv *bam* *flagstat* *log" + out_dir,shell=True)

def set_up(date_of_hiseq, meyer, species, RG_file, output_dir, trim):
	#take all .fastq.gz files in current directory; print them
	files = []

	files = [file for file in os.listdir(".") if file.endswith(".fastq.gz")] 
	
	print "fastq.gz files in current directory:"
	
	print map(lambda x : x ,files)


	#however, if trim option is not yes, then we use fastq files
	if (trim != "yes"):

		files = [file for file in os.listdir(".") if file.endswith("_trimmed.fastq")]

        	print "trimmed fastq files in current directory:"

	        print map(lambda x : x ,files)

	#variables will be initialized here so they can be modified by options 


	#reference genome to be used
	reference = nuclear_genomes[species] 

	#define default cut_adapt

	cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 "

        print "Species selected is " + species

        print "Path to reference genome is " + reference


	#Prepare output directory

	out_dir = output_dir + date_of_hiseq +  "/"

	call("mkdir " + out_dir, shell=True)

	#allow meyer option to be used
	meyer_input = meyer.rstrip("\n").lower()

	alignment_option = "bwa aln -t 5 -l 1024 "  

	if (meyer_input == "meyer"):
		
		print "Meyer option selected."
		
		alignment_option = "bwa aln -t 5 -l 1000 -n 0.01 -o 2 " 

	#variable for RG file
	RG_file = RG_file.rstrip("\n")

	fastq_list = []
	
	#cycle through each line in the input file, gunzip
	for file in files:

		#unzip fastq
		call("gunzip " + file, shell=True)
		current_file = file.split(".")[0] 

		fastq_list.append(current_file.rstrip("\n"))
	
	return files, reference, out_dir, cut_adapt, alignment_option, fastq_list


def trim_fastq(current_sample, cut_adapt, out_dir):
	
	print "Current sample is: " + current_sample
	
	unzipped_fastq = current_sample + ".fastq"
	
	#Get number of lines (and from that reads - divide by four) from raw fastq
	trimmed_fastq = current_sample + "_trimmed" + ".fastq" 
	
	#cut raw fastq files
	call(cut_adapt + unzipped_fastq + " > " + trimmed_fastq + " 2> " + trimmed_fastq + ".log", shell=True)
	

def align(sample, RG_file, alignment_option, reference, trim):

    trimmed_fastq = sample + "_trimmed.fastq"

    if (trim != "yes"):

    	trimmed_fastq = sample + ".fastq"

	sample = "_".join(sample.split("_")[:-1]) 

    print(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai")
    call(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai",shell=True)
    
    with open(RG_file) as file:
	print sample
	for line in file:
		
		split_line = line.split("\t")
		print split_line	
		if (sample == split_line[0].split(".")[0]):

			RG = split_line[1].rstrip("\n")
				
			#check if RG is an empty string
			if not RG:
		
				print "No RGs were detected for this sample - please check sample names in fastq files and in RG file agree" 
				#should probably do something here is there are no read groups
		                break
		        else:
	
				print "Reads groups being used are:"
                               	print RG
	file.seek(0)

	#Print the current sample and RG
        print sample

	print RG                        		
	print "bwa samse -r \'" + RG.rstrip("\n") + "\' " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam"
        call("bwa samse -r \'" + RG.rstrip("\n") + "\' " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)
	
def process_bam(sample_name):

	sample_name = "_".join(sample_name.split("_")[:-1])		

	#sort this bam
	print "samtools sort " + sample_name + ".bam " + sample_name + "_sort"
	call("samtools sort " + sample_name + ".bam " + sample_name + "_sort",shell=True)
	
	print "samtools rmdup -s " + sample_name + "_sort.bam " + sample_name + "_rmdup.bam"
	#remove duplicates from the sorted bam
	call("samtools rmdup -s " + sample_name + "_sort.bam " + sample_name + "_rmdup.bam", shell=True)
	
	#remove the "sorted with duplicates" bam
	call("rm " + sample_name + "_sort.bam",shell=True)

	#make a copy of the samtools flagstat
	call("samtools flagstat " + sample_name + "_rmdup.bam > " + sample_name + "_flagstat.txt",shell=True)

	#remove dupicates from this bam
	#call("samtools view -b -F 4 " + sample_name + "_rmdup.bam > " + sample_name + "_rmdup_only_aligned.bam",shell=True)

	#produce q25 bams
	#call("samtools view -b -q25 " + sample_name + "_rmdup_only_aligned.bam >" + sample_name + "_q25.bam",shell=True)

    	#make a copy of the samtools flagstat
      	#call("samtools flagstat " + sample_name + "_q25.bam > " + sample_name + "_q25.txt",shell=True)
      	
      	#index the q25 bam
	#call("samtools index "+ sample_name + "_q25.bam",shell=True)
	
	#get idx stats
	#call("samtools idxstats "+ sample_name + "_q25.bam > "  + sample_name + ".idx",shell=True)



def merge_lanes_and_sample(RG_file):

	#get sample list from the RG file
	
	sample_list = []
	
        with open(RG_file) as r:

                for line in r:
			
                        sample = line.split("\t")[3].rstrip("\n")
			
			if [sample] not in sample_list:
                        
				sample_list.append([sample])

	print sample_list

	#then cycle through the RG file and associate each lane with the correct sample
        for sample in sample_list:

		lane_list = []

		with open(RG_file) as r:

			for line in r:

				if sample[0] == line.split("\t")[3].rstrip("\n"):
				
					lane = line.split("\t")[2]
				
					if lane not in lane_list:

						lane_list.append(lane)

		sample.append(lane_list)
		print sample


	#create a list of final merged,rmdup bams that will be returned
	merged_bam_list = []
	
	#for each sample, go through each lane for that sample and merge each fastq for that sample
	for sample in sample_list:
				
		merged_lane_list = []
		
		for lane in sample[1]:
		
			files_in_lane = []

			merge_cmd = "java -Xmx20g -jar /research/picard-tools-1.119/MergeSamFiles.jar VALIDATION_STRINGENCY=SILENT "

			with open(RG_file) as r:

				for line in r:

					if (lane == line.split("\t")[2] ) and (sample[0] == line.split("\t")[3].rstrip("\n")):
						print lane
						print sample[0]	
						print line				
						files_in_lane.append(line.split("\t")[0].split(".")[0] + "_rmdup.bam")

			for bam in files_in_lane:

				merge_cmd = merge_cmd + "INPUT=" + bam + " "
				
				if not "-" in bam:
				
					sample_name = bam.split("_")[0]
 
				else:
			
					sample_name = bam.split("-").split("_")[0]


			sample_lane = sample_name + "_"	+ lane + "_merged"		 

			merge_cmd = merge_cmd + "OUTPUT=" + sample_lane + ".bam 2>" + sample_lane + ".log"

			print merge_cmd

			#now merge each bam file associated with a given lane
			call(merge_cmd,shell=True)

			#flagstat the merged bam
			call("samtools flagstat "+ sample_lane + ".bam >"+ sample_lane + "_flagstat.txt" ,shell=True)
		
			#remove duplicates from the merged lane bam
			cmd = "samtools rmdup -s " + sample_lane + ".bam " + sample_lane + "_rmdup.bam 2> " + sample_lane + "_rmdup.log"  

			call(cmd,shell=True)

			#flagstat the rmdup_merged file
			call("samtools flagstat "+ sample_lane + "_rmdup.bam >"+ sample_lane + "_rmdup_flagstat.txt" ,shell=True)

			merged_lane_list.append(sample_lane + "_rmdup.bam")
		
		#each lane has been merged
		#now, merge each lane for a given sample

		merge_cmd = "java -Xmx20g -jar /research/picard-tools-1.119/MergeSamFiles.jar VALIDATION_STRINGENCY=SILENT "
	
		for bam in merged_lane_list:

			merge_cmd = merge_cmd + "INPUT=" + bam + " "
		
		merge_cmd = merge_cmd + "OUTPUT=" + sample[0] + "_merged.bam 2>" + sample[0] + "_merged.log"

		print merge_cmd

		call(merge_cmd,shell=True)

		call("samtools flagstat " + sample[0] + "_merged.bam >" + sample[0]+ "_merged_flagstat.txt",shell=True)

		#now rmdup the sample bam
		call("samtools rmdup -s " + sample[0] + "_merged.bam " + sample[0] + "_merged_rmdup.bam",shell=True)

		call("samtools flagstat " + sample[0] + "_merged_rmdup.bam >" + sample[0]+ "_merged_rmdup_flagstat.txt",shell=True)
		
		merged_bam_list.append(sample[0] + "_merged_rmdup.bam")

	return merged_bam_list

def indel_realignment(rmdup_bam, reference_genome):
	
	print "starting realignment for sample "+rmdup_bam
	
	call("samtools index " + rmdup_bam,shell=True)
	cmd = "java -Xmx20g -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T RealignerTargetCreator -R " + reference_genome + " -I " + rmdup_bam + " -o forIndelRealigner" + rmdup_bam.split(".")[0] + ".intervals 2> " + rmdup_bam.split(".")[0] + "_intervals.log"

	call(cmd, shell=True)
	
	call("java -Xmx20g -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T IndelRealigner -R " + reference_genome + " -I " + rmdup_bam + " -targetIntervals forIndelRealigner" + rmdup_bam.split(".")[0] + ".intervals -o " +  rmdup_bam.split(".")[0] + "_realigned.bam 2> " + rmdup_bam.split(".")[0] + "_realignment.log",shell=True)

	return rmdup_bam.split(".")[0] + "_realigned.bam"


def mapDamage_rescale(bam,reference_genome, out_dir):

	call("mapDamage -i " + bam + " -r " + reference_genome + " --rescale --verbose",shell=True)
	
	call("mv results_" + bam.split(".")[0] + "/" +bam.split(".")[0] + ".rescaled.bam ./ ; mv " + bam.split(".")[0] + ".rescaled.bam " + bam.split(".")[0] + "_rescaled.bam",shell=True)

	call("mkdir "+out_dir+bam.split(".") + "; mv " + bam.split(".")[0] + "_rescaled.bam " + out_dir+bam.split(".")+"/",shell=True)


	print bam.split(".")[0] + "_rescaled.bam"

	return (bam.split(".")[0] + "_rescaled.bam")


def process_realigned_bams(realigned_bam, reference_genome,output_dir):
	
	print realigned_bam
	print realigned_bam.split(".")[0]
	print "samtools rmdup -s " + realigned_bam + " " + realigned_bam.split(".")[0] + "_rmdup.bam && samtools flagstat " + realigned_bam.split(".")[0] + "_rmdup.bam 2> " + realigned_bam.split(".")[0] + "_rmdup_flagstat.txt"
	call("samtools rmdup -s " + realigned_bam + " " + realigned_bam.split(".")[0] + "_rmdup.bam && samtools flagstat " + realigned_bam.split(".")[0] + "_rmdup.bam 2> " + realigned_bam.split(".")[0] + "_rmdup_flagstat.txt",shell=True)	

	call("samtools view -b -F4 " + realigned_bam.split(".")[0] + "_rmdup.bam > " + realigned_bam.split(".")[0] + "_rmdup_F4.bam && samtools view -q25 -b " + realigned_bam.split(".")[0] + "_rmdup_F4.bam > " + realigned_bam.split(".")[0] + "_rmdup_q25.bam",shell=True)

	call("samtools flagstat " + realigned_bam.split(".")[0] + "_rmdup_q25.bam > " + realigned_bam.split(".")[0] + "_rmdup_q25_flagstat.txt",shell=True)

	call("samtools index " + realigned_bam.split(".")[0] + "_rmdup_q25.bam",shell=True)

	cmd="java -Xmx20g -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T DepthOfCoverage -R " + reference_genome + " -o DoC_" + realigned_bam.split(".")[0] + " -I " + realigned_bam.split(".")[0] + "_rmdup_q25.bam --omitDepthOutputAtEachBase"

	call(cmd,shell=True)
	
	q25_bam=realigned_bam.split(".")[0] + "_rmdup_q25.bam"
	
	mapDamage_rescale(q25_bam,reference_genome,output_dir)
	
	call("mv DoC* " + output_dir,shell=True)
try:
	date_of_hiseq  = sys.argv[1]
	meyer = sys.argv[2]
	species = sys.argv[3]
	trim = sys.argv[4]
	align = sys.argv[5]
	process = sys.argv[6]
	RG_file  = sys.argv[7]
	output_dir = sys.argv[8]

except IndexError:
	print "Incorrect number of variables have been provided"
	print "Input variables are date_of_hiseq, meyer, species, trim, process, RG_file, and the directory to put output directories/files"
	print "Exiting program..."
	sys.exit()


if not output_dir[-1] == "/":

	output_dir = output_dir + "/"

main(date_of_hiseq, meyer, species, trim, align, process, RG_file, output_dir)
