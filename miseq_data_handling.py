#!/usr/bin/env python
'''
Script is run in a directory with bam files to be aligned to the goat genome
Bams will be trimmed, quality metrics obtained, screened using fastq screen

Also need to supply a date for the Miseq run - will automatically make results file
python script <date_of_miseq> <meyer> <species> <mit> <fastq_screen> <read group file>

read group file needs to be in the follow format: 
Sample_name\tRGs_to_add (in following format- ID:X\tSM:X\tPL:X\tLB:X)


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

	"goat" : "~/goat/miseq/data/reference_genomes/goat_CHIR1_0/goat_CHIR1_0.fasta",

	"sheep" : "~/goat/miseq/data/reference_genomes/sheep_oviAri3/oviAri3.fa",

	"cow" : "~/goat/miseq/data/reference_genomes/cow_bosTau8/bosTau8.fa",

	"dog" : "~/goat/miseq/data/reference_genomes/canFam3.fa"
}

mitochondrial_genomes = {

	"goat" : "~/goat/miseq/data/mit_genomes_fastq_screen/goat_mit/goat_mit.fasta",

	"sheep" : "~/goat/miseq/data/mit_genomes_fastq_screen/sheep_mit/sheep_mit.fasta"

}

def main(date_of_miseq, meyer, species, mit, fastq_screen, RG_file):
	
	#run the set up function.#set up will create some output directories
	#and return variables that will be used in the rest of the script
	
	files, reference, out_dir, cut_adapt, alignment_option, master_list, sample_list, fastq_screen_option = set_up(date_of_miseq, meyer, species, mit, RG_file) 
	
	#sample is the file root
	#trim fastq files and produce fastqc files
	#the masterlist will change each time so it needs be equated to the function
	
	for sample in sample_list:
		
		trim_fastq(sample, cut_adapt, out_dir)
	
	#run fastq screen on the samples if input variable for fastq_screen is "yes"
	if (fastq_screen.lower() == "yes"):

		map(lambda sample: run_fastq_screen(sample, out_dir, fastq_screen_option), sample_list)
	
	#at this stage we have our fastq files with adaptors trimmed, fastqc and fastq screen run
	#we can now move on to the next step: alignment
	
	#going to align to CHIR1.0, as that what was used for AdaptMap

	map(lambda sample : align(sample, RG_file, alignment_option, reference), sample_list)
	
	#testing a function here, to process a bam to a q25 version
	map(process_bam, sample_list)
	
	for sample in sample_list:

		master_list = get_summary_info(master_list, sample)
	
	call("mkdir trimmed_fastq_files_and_logs",shell=True)
	call("mv *trimmed* trimmed_fastq_files_and_logs/",shell=True)
		
	for sample in sample_list:
	
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

	output_summary = date_of_miseq + "_summary.table"
	
	#print summary stats
	with open(output_summary, "w") as f:

		writer = csv.writer(f, delimiter='\t', lineterminator='\n')
		writer.writerows(master_list)

	call("wc -l " + output_summary,shell=True)
	
	number_of_samples = (int((subprocess.check_output("wc -l " + output_summary,shell=True).split(" ")[0])) - 1)
	
	call("head -n1 " + output_summary + "> header.txt ",shell=True)
	
	call("tail -n " + str(number_of_samples) + " " + output_summary + " | sort | cat header.txt - > tmp; mv tmp " + output_summary + ";rm header.txt tmp",shell=True)

	call("mv " + output_summary + " " + out_dir,shell=True)


def set_up(date_of_miseq, meyer, species, mit, RG_file):
	#take all .fastq.gz files in current directory; print them
	files = []

	files = [file for file in os.listdir(".") if file.endswith(".fastq.gz")] 
	
	print "fastq.gz files in current directory:"
	
	print map(lambda x : x ,files)

	#variables will be initialized here so they can be modified by options 


	if (mit != "yes"): 
		
		#reference genome to be used
		reference = nuclear_genomes[species] 

		#define default cut_adapt and fastq screen

	        cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 "

        	fastq_screen = "~/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --outdir ./"

	else:

		print "Mitochondrial alignment selected"

		reference = mitochondrial_genomes[species]
	
		cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 25 "

                fastq_screen = "~/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --conf ~/goat/miseq/src/fastq_screen_v0.4.4/capture_config/capture.conf --outdir ./"


	print "Species selected is " + species

        print "Path to reference genome is " + reference


	#Prepare output directory

	out_dir = "~/goat/miseq/results/" + date_of_miseq +  "/"

	if (species == "sheep"):

		out_dir = "~/sheep/results/" + date_of_miseq +  "/" 

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
	master_list = [["Sample", "read_count_raw", "trimmed_read_count","raw_reads_aligned", "raw %age endogenous", "rmdup_reads_remaining","rmdup_reads_aligned" ,"rmdup_alignment_percent", "reads_aligned_q25", "percentage_reads_aligned_q25"]]

	sample_list = []
	#cycle through each line in the input file, gunzip
	for file in files:

		#unzip fastq
		call("gunzip " + file, shell=True)
		current_file = file.split(".")[0] + ".fastq"

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
	
	print "Current sample is: " + current_sample
	
	unzipped_fastq = current_sample + ".fastq"
	
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
       	
       	
def run_fastq_screen(current_sample, out_dir, fastq_screen_option):

	call("mkdir " + out_dir + "fastq_screen/",shell=True)
	
	call("mkdir " + out_dir + "fastq_screen/" + current_sample, shell=True)
	
	call(fastq_screen_option + current_sample + " " + current_sample + "_trimmed.fastq --outdir " + out_dir + "fastq_screen/" + current_sample, shell=True)


def align(sample, RG_file, alignment_option, reference):

    trimmed_fastq = sample + "_trimmed.fastq"

    print(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai")
    call(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai",shell=True)
    
    with open(RG_file) as file:
	print sample
	for line in file:
		
		split_line = line.split("\t")
		print split_line	
		if (sample == split_line[0]):

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

        print sample
	print RG                        		
        call("bwa samse "  + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)
	call("cp " + sample + ".bam " +sample+"_cp.bam",shell=True)
	#add the read groups
	print "java -jar /research/picard-tools-1.119/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=SILENT  INPUT=" + sample + ".bam OUTPUT=" + sample + "_RG.bam " + RG.rstrip("\n")
	call("java -jar /research/picard-tools-1.119/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=SILENT  INPUT=" + sample + ".bam OUTPUT=" + sample + "_RG.bam " + RG.rstrip("\n"), shell=True)
	
	print "mv " + sample + "_RG.bam " +  sample + ".bam "
	call("mv " + sample + "_RG.bam " +  sample + ".bam ",shell=True)	

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

	#get q25 reads aligned
	q25_reads_aligned = subprocess.check_output("more " + current_sample + "_q25_flagstat.txt | grep 'mapped (' | cut -f1 -d' '",shell=True)
       	to_add.append(q25_reads_aligned.rstrip("\n"))

	#q25_percent_aligned = subprocess.check_output("more " + sample + "_q25_flagstat.txt | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	fixed_percentage = str(((float(q25_reads_aligned)) * 100)/ float(trimmed_read_number))
	to_add.append(fixed_percentage.rstrip("\n"))

	
	for i in master_list:
		
		print i		
		if (i[0] == current_sample):
			i.extend(to_add)		
			break
		
	return master_list

try:
	date_of_miseq  = sys.argv[1]
	meyer = sys.argv[2]
	species = sys.argv[3]
	mit = sys.argv[4]
	fastq_screen = sys.argv[5]
	RG_file  = sys.argv[6]
except IndexError:
	print "Incorrect number of variables have been provided"
	print "Input variables are date_of_miseq, meyer, species, mit, fastq_screen, RG_file"
	print "Exiting program..."
	sys.exit()


main(date_of_miseq, meyer, species, mit, fastq_screen, RG_file)
