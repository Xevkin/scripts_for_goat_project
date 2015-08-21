#!/usr/bin/python
"""
Script is run in a directory with bam files to be aligned to the goat genome
Bams will be trimmed, quality metrics obtained, screened using fastq screen

Also need to supply a date for the Miseq run - will automatically make results file
python script <date_of_miseq> <meyer> <option> <read group file>

read group file needs to be in the follow format: 
Sample_name\tRGs_to_add (in following format- ID:X\tSM:X\tPL:X\tLB:X)



"""

#import csv module to easily write the list of list used as a .csv file
import csv

#apparently this is a better way to make system calls using python, rather than "os.system"
import subprocess 
from subprocess import call

#need to import sys anyways to access input file (a list)
import sys

#we will use this later to check if the input files actually exist
import os

def main(date_of_miseq, meyer, option, RG_file):
	
	#run the set up function.#set up will create some output directories
	#and return variables that will be used in the rest of the script
	
	files, reference, out_dir, cut_adapt, alignment_option, master_list, sample_list, fastq_screen_option = set_up(date_of_miseq, meyer, option) 
	
	#sample is the file root
	#trim fastq files and produce fastqc files
	#the masterlist will change each time so it needs be equated to the function
	
	for sample in sample_list:
		
		trim_fastq(sample, cut_adapt, out_dir)
	
	#run fastq screen on the samples

	map(lambda sample: run_fastq_screen(sample, out_dir, fastq_screen_option), sample_list)
	
	#at this stage we have our fastq files with adaptors trimmed, fastqc and fastq screen run
	#we can now move on to the next step: alignment
	
	#going to align to CHIR1.0, as that what was used for AdaptMap
	
	map(lambda sample :align(sample,ＲＧ＿ｆｉｌｅ， alignment_option, reference), sample_list)
	
	#testing a function here, to process a bam to a q25 version
	map(process_bam, sample_list)
		
	ｆｏｒ　ｓａｍｐｌｅ　ｉｎ　ｓａｍｐｌｅｓ：
	
		ｍａｓｔｅｒ＿ｌｉｓｔ　＝　ｇｅｔ＿ｓｕｍｍａｒｙ＿ｉｎｆｏ（ｍａｓｔｅｒ＿ｌｉｓｔ，　ｓａｍｐｌｅ）
	
	call("mkdir trimmed_fastq_files_and_logs",shell=True)
	call("mv *trimmed* trimmed_fastq_files_and_logs/",shell=True)
		
	ｆｏｒ　ｓａｍｐｌｅ　ｉｎ　ｓａｍｐｌｅ＿ｌｉｓｔ：
	
		#clean up files
		call("gzip "+ sample + "*",shell=True)
	
		#going to make an output directory for each sample
		#then move all produced files to this directory
		call("mkdir " + out_dir + sample,shell=True)
		call("mv *" + sample +＂*.bam*　＂+ sample + "*.idx* "+ sample＋　"*flagstat* " + out_dir + sample,shell=True)

	#remove all .sai files
	call("rm *sai*",shell=True)
	call("mv trimmed_fastq_files_and_logs/ " + out_dir,shell=True)

	output_summary = date_of_miseq + "_summary.table"
	
	#print summary stats
	with open(output_summary, "w") as f:

		writer = csv.writer(f, delimiter='\t', lineterminator='\n')
		writer.writerows(master_list)
	
	number_of_samples = int(subprocess.check_output("wc -l " + output_summary,shell=True)) - 1

	call("head -n1 " + output_summary "> header.txt ",shell=True)
	
	call("tail -n " + str(number_of_samples) "output_summary | sort | cat header.txt - > tmp; mv tmp " + output_summary ";rm header.txt tmp",shell=True)

	call("mv " + output_summary + " " + out_dir,shell=True)


def set_up(date_of_miseq, meyer, option):
	#take all .fastq.gz files in current directory; print them
	files = []

	files = [file for file in os.listdir(".") if file.endswith(".fastq.gz")] 
	
	print "fastq.gz files in current directory:"
	
	print map(lambda x : x ,files)

	#variables will be initialized here so they can be modified by options 

	#reference genomes
	reference = "~/goat/miseq/data/reference_genomes/goat_CHIR1_0/goat_CHIR1_0.fasta"

	#Prepare output directory

	out_dir = "~/goat/miseq/results/" + date_of_miseq +  "/"

	call("mkdir " + out_dir, shell=True)

	#allow meyer option to be used
	meyer_input = meyer.rstrip("\n").lower()

	alignment_option = "bwa aln -t 5 -l 1000 "  

	if (meyer_input == "meyer"):
		print "Meyer option selected."
		alignment_option = "bwa aln -t 5 -l 1000 -n 0.01 -o 2 " 

	#turn the reference genome path to the sheep if sheep option is selected
	option = option.rstrip("\n").lower()

	if (option == "sheep"):
		print "Sheep alignment selected."
		reference = "~/goat/miseq/data/reference_genomes/sheep_oviAri3/oviAri3.fa" 

	#define default cut_adapt and fastq screen

	cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 "

	fastq_screen = "~/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --outdir ./"
	
	#if option is mit, then several changes need to occur

	if (option == "mit"):

		print "Mitochondrial alignment selected."

		cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 25 "

		reference = "~/goat/miseq/data/mit_genomes_fastq_screen/goat_mit/goat_mit.fasta"

		fastq_screen = "~/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --conf ~/goat/miseq/src/fastq_screen_v0.4.4/capture_config/capture.conf --outdir ./"

	#variable for RG file
	RG_file = RG_file.rstrip("\n")

	#initialize a masterlist that will carry summary stats of each sample
	master_list = []
	master_list.append([["Sample", "wc-l", "read_count_raw", "wc-l_trimmed", "trimmed_read_count","raw_reads_aligned", "raw %age endogenous", "rmdup_reads_aligned" ,"rmdup_alignment_percent", "reads_aligned_q25", "percentage_reads_aligned_q25"]])

	sample_list = []
	#cycle through each line in the input file, gunzip
	for file in files:

		#unzip fastq
		call("gunzip " + current_file, shell=True)
		current_file = current_file.split(".")[0] + ".fastq"

		#rename the file
		#also save the current sample as a variable to be used later
		split_file = current_file.split(".")[0].split("_")
		
		call("mv " + current_file + " " + split_file[0] + ".fastq",shell=True)
	
		#current_file = split_file[0] + ".fastq"
	
		sample_list.append(split_file[0].rstrip("\n"))
		
		for i in sample_list:
		
			master_list.append([i])

	return files, reference, out_dir, cut_adapt, alignment_option, master_list, sample_list, fastq_screen


def trim_fastq(current_sample, cut_adapt, out_dir):
	
	print "Current samples is: " + sample
	
	unzipped_fastq = current_sample＋".fastq"
	
	#Get number of lines (and from that reads - divide by four) from raw fastq
	trimmed_fastq = current_sample + "_trimmed" + ".fastq" 
	cmd = "wc -l " + unzipped_fastq + " | cut -f1 -d' '" 
	
	
	#cut raw fastq files
	call(cut_adapt + unzipped_fastq + " > " + trimmed_fastq + " 2> " + trimmed_fastq + ".log", shell=True)
	
       	#run fastqc on both the un/trimmed fastq files
	#first we want to create an output directory if there is none to begin with
	call("mkdir " +  out_dir + "fastqc/", shell=True)
		
	call("fastqc " + unzipped_fastq + " -o " + out_dir + "fastqc/", shell=True)
	call("fastqc " + trimmed_fastq + " -o " + out_dir + "fastqc/", shell=True)
       	
       	return master_list
       	
def run_fastq_screen(current_sample, out_dir, fastq_screen_option):

	call("mkdir " + out_dir + "fastq_screen/",shell=True)
	
	call("mkdir " + out_dir + "fastq_screen/" + current_sample, shell=True)
	
	call(fastq_screen + current_sample + " " + current_sample + "_trimmed.fastq --outdir " + out_dir + "fastq_screen/" + sample, shell=True)


def align(sample, RG_file, alignment_option, reference):

    trimmed_fastq = sample + "_trimmed.fastq"

    print(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai")
    call(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai",shell=True)
    
    with open(RG_file) as file:

	for line in file:

		split_line = line.split("\t")
			
			if (sample == split_line[0]):

				RG = split_line[1].rstrip("\n")
				print RG
				
				#check if RG is an empty string
				if not RG:
		
					print "No RGs were detected for this sample - please check sample names in fastq files and in RG file agree" 
			                          #should probably do something here is there are no read groups
			                break
			        else:
		
					print "Reads groups being used are:"
                                	print RG
                                		
        #produce bam with all reads
	#at this stage we also add read groups
	call("bwa samse -r \'@RG\t" + RG + "\t\' " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)


def process_bam(sample_name):
		
	#sort this bam
	call("samtools sort " + sample_name + ".bam " + sample_name + "_sort",shell=True)
	
	#remove duplicates from the sorted bam
	call("samtools rmdup -s " + sample_name + "_sort.bam " + sample_name + "_rmdup.bam", shell=True)
	
	#remove the "sorted with duplicates" bam
	call("rm " + sample_name + "_sort.bam",shell=True)

	#make a copy of the samtools flagstat
	call("samtools flagstat " + sample_name + "_rmdup.bam > " + sample_name + "_flagstat.txt",shell=True)

	#remove unaligned reads from this bam
	call("samtools view -b -F 4 " + sample_name + "_rmdup.bam > " + sample_name + "_rmdup_only_aligned.bam",shell=True)

	#produce q25 bams
	call("samtools view -b -q25 " + sample_name + "_rmdup_only_aligned.bam >" + sample_name + "_q25.bam",shell=True)

	#sort this bam
       	call("samtools sort " + sample_name + "_q25.bam " + sample_name + "_q25_sort",shell=True)

       	#remove duplicates from the sorted bam
       	call("samtools rmdup -s " + sample_name + "_q25_sort.bam " + sample_name + "_q25_rmdup.bam", shell=True)

       	#remove the "sorted with duplicates" bam
       	call("rm " + sample_name + "_q25_sort.bam",shell=True)

      	#make a copy of the samtools flagstat
      	call("samtools flagstat " + sample_name + "_q25_rmdup.bam > " + sample_name + "_q25_flagstat.txt",shell=True)
      	
      	#index the q25 bam
	call("samtools index "+ sample_name + "_q25_rmdup.bam",shell=True)
	
	#get idx stats
	call("samtools idxstats "+ sample_name + "_q25_rmdup.bam > "  + sample_name + ".idx",shell=True)

def　get_summary_info(master_list, current_sample)：

	to_add　＝［］
	
	#raw reads
	raw_read_number = int(subprocess.check_output(cmd,shell=True)) / 4
	
	to_add．ａｐｐｅｎｄ（raw_read_number）
	
	#grab summary statistics of trimmed file
	cmd = "wc -l " + trimmed_fastq + "| cut -f1 -d' '"
       	trimmed_read_number = int(subprocess.check_output(cmd,shell=True)) / 4
       	raw_read_number.append(subprocess.check_output(ｓｔｒ（cmd）,shell=True))
       	
       	#get number of reads aligned without rmdup
	raw_reads_aligned = subprocess.check_output("samtools flagstat " + sample + ".bam |  grep 'mapped (' | cut -f1 -d' '",shell=True)
	to_add.append(raw_reads_aligned)

	#get %age raw alignment
	raw_alignment_percentage = ((float(raw_reads_aligned)) * 100)/ float(trimmed_read_number)
	to_add.append(ｓｔｒ（raw_alignment_percentage)）

	#get reads that aligned following rmdup
	rmdup_reads_aligned = subprocess.check_output("more " + sample + "_flagstat.txt | grep 'mapped (' | cut -f1 -d' '",shell=True)
	to_add.append(rmdup_reads_aligned)

  	#capture the alignment percentage of the flagstat file, both no q and q30
	raw_alignment = subprocess.check_output("more " + sample + "_flagstat.txt | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	to_add.append(raw_alignment)

	#get q25 reads aligned
	q25_reads_aligned = subprocess.check_output("more " + sample + "_q25_flagstat.txt | grep 'mapped (' | cut -f1 -d' '",shell=True)
       	to_add.append(q25_reads_aligned)

	#q25_percent_aligned = subprocess.check_output("more " + sample + "_q25_flagstat.txt | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	fixed_percentage = ｓｔｒ（((float(q25_reads_aligned)) * 100)/ float(trimmed_read_number)）
	to_add.append(fixed_percentage)

	
	for i in master_list:
		
		if (i[0] == current_sample):
			
			i.ｅｘｔｅｎｄ(ｔｏ＿ａｄｄ)
	
			break
		
	ｒｅｔｕｒｎ　ｍａｓｔｅｒ＿ｌｉｓｔ

	
main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
