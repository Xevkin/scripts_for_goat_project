#!/usr/bin/env python
'''
Script is run in a directory with bam files to be aligned to either goat, wild goat or sheep genome
Reads will be trimmed, aligned, read groups added, duplicates removed, then 

Also need to supply a date for the Hiseq run - will automatically make results file
python script <date_of_hiseq> <meyer> <species> <fastq_screen> <read group file> <directory in which to place output dir/>

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


def main(date_of_hiseq, meyer, species, RG_file, output_dir):
	
	#run the set up function.#set up will create some output directories
	#and return variables that will be used in the rest of the script
	
	files, reference, out_dir, cut_adapt, alignment_option, master_list, fastq_list = set_up(date_of_hiseq, meyer, species, RG_file, output_dir) 
	
	#sample is the file root
	#trim fastq files and produce fastqc files
	#the masterlist will change each time so it needs be equated to the function
	
	for fastq in fastq_list:
		
		trim_fastq(fastq, cut_adapt, out_dir)
	
	#at this stage we have our fastq files with adaptors trimmed, fastqc and fastq screen run
	#we can now move on to the next step: alignment
	
	#going to align to CHIR1.0, as that what was used for AdaptMap
	#this step will align each fastq and produce raw bams with RGs
	#they should have the same stem of the initial bam

	##########map(lambda fastq : align(fastq, RG_file, alignment_option, reference), fastq_list)


	#sort and remove duplicates from each bam	
	map(process_bam, fastq_list)

        for sample in fastq_list:

                master_list = get_summary_info(master_list, sample)

	
	#get sample list from the RG file
	
	sample_list = []
	
        with open(RG_file) as r:

                for line in r:
			
                        sample = line.split("\t")[3].rstrip("\n")
			
                        if not sample in sample_list:

                                sample_list.append([sample])

	print sample_list

	#then cycle through the RG file and associate each lane with the correct sample
        for sample in sample_list:

		lane_list = []

		with open(RG_file) as r:

			for line in r:

				if sample[0] == line.split("\t")[3].rstrip("\n"):
				
					lane = line.split("\t")[2]
				
					if not lane in lane_list:

						lane_list.append(lane)

		sample.append(lane_list)
	print lane_list
	#for each sample, go through each lane for that sample and merge each fastq for that sample
	for sample in sample_list:
		
		merged_lane_list = []
		
		for lane in sample[1]:
		
			files_in_lane = []

			merge_cmd = "java -Xmx20g -jar /research/picard-tools-1.119/MergeSamFiles.jar VALIDATION_STRINGENCY=SILENT "

			with open(RG_file) as r:

				for line in r:
	
					if lane == line.split("\t")[2]:

						files_in_lane.append(line.split("\t")[0] + "_rmdup.bam")

			for bam in files_in_lane:

				merge_cmd = merge_cmd + "INPUT=" + bam + " "

				sample_name = bam.split("-")[1].split("_")[0]
 
			sample_lane = sample_name + "_"	+ lane + "_merged"		 

			merge_cmd = merge_cmd + "OUTPUT=" + sample_lane + ".bam 2>" + sample_lane + ".log"
			print merge_cmd
			#now merge each bam file associated with a given lane
			call(merge_cmd,shell=True)

			#flagstat the merged bam
			call("samtools flagstat "+ sample_lane + ".bam >"+ sample_lane + "_flagstat.txt" ,shell=True)
		
			#remove duplicates from the merged lane bam
			cmd = "samtools rmdup -s " + sample_lane + ".bam " + sample_lane + "_rmdup.bam"

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

	
	
	call("mkdir trimmed_fastq_files_and_logs",shell=True)
	call("mv *trimmed* trimmed_fastq_files_and_logs/",shell=True)
		
	for sample in fastq_list:
	
		#clean up files
		call("gzip "+ sample + "*",shell=True)
	
		#going to make an output directory for each sample
		#then move all produced files to this directory
		call("mkdir " + out_dir + sample,shell=True)
		print "mv *" + sample + "*.bam* "+ sample + "*.idx* "+ sample + "*flagstat* " + out_dir + sample
		call("mv *" + sample + "*.bam* "+ sample + "*.idx* "+ sample + "*flagstat* " + out_dir + sample,shell=True)

	#remove all .sai files
	call("rm *sai*",shell=True)
	call("mv trimmed_fastq_files_and_logs/ " + out_dir,shell=True)

	output_summary = date_of_hiseq + "_summary.table"
	
	call("mv *bam " + out_dir,shell=True)

	#print summary stats
	#with open(output_summary, "w") as f:

	#	writer = csv.writer(f, delimiter='\t', lineterminator='\n')
	#	writer.writerows(master_list)

	#call("wc -l " + output_summary,shell=True)
	
	#number_of_samples = (int((subprocess.check_output("wc -l " + output_summary,shell=True).split(" ")[0])) - 1)
	
	#call("head -n1 " + output_summary + "> header.txt ",shell=True)
	
	#call("tail -n " + str(number_of_samples) + " " + output_summary + " | sort | cat header.txt - > tmp; mv tmp " + output_summary + ";rm header.txt tmp",shell=True)

	#call("mv " + output_summary + " " + out_dir,shell=True)


def set_up(date_of_hiseq, meyer, species, RG_file, output_dir):
	#take all .fastq.gz files in current directory; print them
	files = []

	files = [file for file in os.listdir(".") if file.endswith(".fastq.gz")] 
	
	print "fastq.gz files in current directory:"
	
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

	alignment_option = "bwa aln -t 5 -l 1000 "  

	if (meyer_input == "meyer"):
		
		print "Meyer option selected."
		
		alignment_option = "bwa aln -t 5 -l 1000 -n 0.01 -o 2 " 

	#variable for RG file
	RG_file = RG_file.rstrip("\n")

	#initialize a masterlist that will carry summary stats of each sample
	master_list = [["Sample", "read_count_raw", "trimmed_read_count","raw_reads_aligned", "raw %age endogenous", "rmdup_reads_remaining","rmdup_reads_aligned" ,"rmdup_alignment_percent"]]

	fastq_list = []
	
	#cycle through each line in the input file, gunzip
	for file in files:

		#unzip fastq
		call("gunzip " + file, shell=True)
		current_file = file.split(".")[0] 

		fastq_list.append(current_file.rstrip("\n"))
	
	for i in fastq_list:
	
		master_list.append([i])
	
	return files, reference, out_dir, cut_adapt, alignment_option, master_list, fastq_list


def trim_fastq(current_sample, cut_adapt, out_dir):
	
	print "Current sample is: " + current_sample
	
	unzipped_fastq = current_sample + ".fastq"
	
	#Get number of lines (and from that reads - divide by four) from raw fastq
	trimmed_fastq = current_sample + "_trimmed" + ".fastq" 
	
	#cut raw fastq files
	call(cut_adapt + unzipped_fastq + " > " + trimmed_fastq + " 2> " + trimmed_fastq + ".log", shell=True)
	

def align(sample, RG_file, alignment_option, reference):

    trimmed_fastq = sample + "_trimmed.fastq"

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
		
	#sort this bam
	call("samtools sort " + sample_name + ".bam " + sample_name + "_sort",shell=True)
	
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


def get_summary_info(master_list, current_sample):

	print "Initial master list is:"
	print master_list
	to_add = []

	unzipped_fastq = current_sample + ".fastq"
	
	trimmed_fastq = current_sample + "_trimmed.fastq"
	
	cmd = "wc -l " + unzipped_fastq + " | cut -f1 -d' '" 
	
	#raw reads
	file_length = subprocess.check_output(cmd,shell=True)
	
	raw_read_number = int(file_length) / 4
	
	to_add.append(raw_read_number)
	
	#grab summary statistics of trimmed file
	cmd = "wc -l " + trimmed_fastq + "| cut -f1 -d' '"

       	file_length = subprocess.check_output(cmd,shell=True)

	trimmed_read_number = int(file_length) / 4
       	
	to_add.append(str(trimmed_read_number).rstrip("\n"))

       	#get number of reads aligned without rmdup
	raw_reads_aligned = subprocess.check_output("samtools flagstat " + current_sample + ".bam |  grep 'mapped (' | cut -f1 -d' '",shell=True)
	to_add.append(raw_reads_aligned.rstrip("\n"))

	#get %age raw alignment
	raw_alignment_percentage = ((float(raw_reads_aligned)) * 100)/ float(trimmed_read_number)

	to_add.append(str(raw_alignment_percentage).rstrip("\n"))

	rmdup_reads_remaining = subprocess.check_output("more " + current_sample + "_flagstat.txt | head -n1 | cut -f1 -d' '",shell=True) 

	to_add.append(rmdup_reads_remaining.rstrip("\n"))
	#get reads that aligned following rmdup
	rmdup_reads_aligned = subprocess.check_output("more " + current_sample + "_flagstat.txt | grep 'mapped (' | cut -f1 -d' '",shell=True)
	to_add.append(rmdup_reads_aligned.rstrip("\n"))

  	#capture the alignment percentage of the flagstat file, both no q and q30
	raw_alignment = subprocess.check_output("more " + current_sample + "_flagstat.txt | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	to_add.append(raw_alignment.rstrip("\n"))
	
	for i in master_list:
		
		print i		
		if (i[0] == current_sample):
			i.extend(to_add)		
			break
		
	return master_list


def merge_lanes()


try:
	date_of_hiseq  = sys.argv[1]
	meyer = sys.argv[2]
	species = sys.argv[3]
	RG_file  = sys.argv[4]
	output_dir = sys.argv[5]

except IndexError:
	print "Incorrect number of variables have been provided"
	print "Input variables are date_of_hiseq, meyer, species, RG_file, and the directory to put output directories/files"
	print "Exiting program..."
	sys.exit()


main(date_of_hiseq, meyer, species, RG_file, output_dir)
