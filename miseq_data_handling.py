#!/bin/python


#also need to supply a date for the Miseq run - will automatically make results file
#python script <date_of_miseq> <meyer> <sheep> <read group file>

#read group file needs to be in the follow format: 
#Sample_name\tRGs_to_add (in following format- ID:X\tSM:X\tPL:X\tLB:X)

#import csv module to easily write the list of list used as a .csv file
import csv

#apparently this is a better way to ake system calls using python, rather than "os.system"
import subprocess 
from subprocess import call

#need to import sys anyways to access intput file (a list)
import sys

#we will use this later to check if the input files actually exist
import os

#take all .fastq.gz files in current directory; print them

files = []

for file in os.listdir("."):

        if file.endswith(".fastq.gz"):

                files.append(file)
print "fastq.gz files in current directory:"
for file in files:
	print file


#reference genomes
reference = "~/goat/miseq/data/reference_genomes/goat_CHIR1_0/goat_CHIR1_0.fasta"

#Prepare output directory

miseq_date = sys.argv[1]

out_dir = "~/goat/miseq/results/" + miseq_date + "/"

call("mkdir " + out_dir, shell=True)

#allow meyer option to be used
meyer_input = sys.argv[2].rstrip("\n").lower()

alignment_option = "bwa aln -l 1000 "  

if (meyer_input == "meyer"):
	print "Meyer option selected."
	alignment_option = "bwa aln -l 1000 -n 0.01 -o 2 " 

#turn the reference genome path to the sheep if sheep option is selected
sheep_option = sys.argv[3].rstrip("\n").lower()

if (sheep_option == "sheep"):
	print "Sheep alignment selected."
	reference = "~/goat/miseq/data/reference_genomes/sheep_oviAri3/oviAri3.fa" 

#variable for RG file
RG_file = sys.argv[4].rstrip("\n")

#initialize a masterlist that will carry summary stats of each sample
master_list = []
master_list.append(["Sample", "wc-l", "read_count_raw", "wc-l_trimmed", "trimmed_read_count","raw_reads_aligned","rmdup_reads_aligned" ,"rmdup_alignment_percent", "reads_aligned_q30", "percentage_reads_aligned_q30"])

#create list that will carry any lines in the file which fail at any stage
failures = []

#cycle through each line in the input file
for file in files:
	
	current_file = file.strip()
	#if the line in the file is not a .fastq.gz, report an error and skip to next one
	#record those which were not successful

	if not current_file.endswith(".fastq.gz"):

		print current_file + " is not a fastq.gz file."
	
		failures.append(current_file)
	
		continue
	
	#check if file actually exists in this directory
	#if not print error message and record it
	if not os.path.isfile(current_file):

		print current_file + " is not in the current directory."
		failures.append(current_file)

		continue
	#using call([]) to invoke system
	#cutadapt with the adaptor already defined
	#minimum overlap before trimming (-O 1) and minimum sequence length before trimming (-m 30) are also preset
	#note that the -m 30 option is not appropriate for mitochondrial captures
	
	#unzip fastq
	#print "gunzip " + current_file
	call("gunzip " + current_file, shell=True)
	
	#seperate fastq.gz file name so the file names generated during the process can be used
	split_file = current_file.split(".")
	
	sample = split_file[0].rstrip("\n")
	print sample
	
	#initialize variable to carry summary statistics
	summary = []
	summary.append(sample)
	#wc -l raw fastq
	unzipped_fastq = sample + "." + split_file[1]
	trimmed_fastq = sample + "_trimmed" + "." + split_file[1]
	
	cmd = "wc -l " + unzipped_fastq + " | cut -f1 -d' '" 
	summary.append(subprocess.check_output(cmd,shell=True))
	
	call(cmd,shell=True)
	#raw reads
	raw_read_number = int(subprocess.check_output(cmd,shell=True)) / 4
	summary.append(raw_read_number)

	call("cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 " + unzipped_fastq + " > " + trimmed_fastq + " 2> " + trimmed_fastq + ".log", shell=True)
	#print "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 " + unzipped_fastq + " > " + trimmed_fastq
	
	#rezip original fastq files

	call("gzip " + unzipped_fastq, shell=True)
	
	#grab summary statistics of trimmed file
	cmd = "wc -l " + trimmed_fastq + "| cut -f1 -d' '"
       	summary.append(subprocess.check_output(cmd,shell=True))
       	trimmed_read_number = int(subprocess.check_output(cmd,shell=True)) / 4
       	summary.append(str(trimmed_read_number))


	#run fastqc on trimmed fastq file
	#first we want to create an output directory if there is none to begin with
			
	call("mkdir " +  out_dir + "fastqc/", shell=True)
	call("fastqc " + trimmed_fastq + " -o " + out_dir + "fastqc/", shell=True)

	call("mkdir " + out_dir + "fastq_screen/" + sample, shell=True)
	call("mkdir ./" + sample, shell=True)
	
	call("~/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --outdir ./" + sample + " " + trimmed_fastq, shell=True)
	print "Output of fastq_screen:"
	call("ls ./" + sample,shell=True)
	#call("rsync --remove-source-files ./" + sample + " " + out_dir + "fastq_screen/" + sample,shell=True)
	call("mv ./" + sample + " " +  out_dir + "fastq_screen/",shell=True)
	call("rm -r ./" + sample,shell=True)

	#at this stage we have our fastq files with adaptors trimmed
	#we can now move on to the next step: alignment
	#going to align to CHIR1.0, as that what was used for AdaptMap
	print(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai")
	call(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai",shell=True)	
	
	#Obtain the appropriate read group from the supplied read group file
	with open(RG_file) as file:

		for line in file:

			split_line = line.split("\t")
			print sample
			print split_line
			if (sample == split_line[0]):

				RG = split_line[1].rstrip("\n")

		#check if RG is an empty string
		if not RG:
		
			print "No RGs were detected for this sample - please check sample names in fastq files and in RG file agree" 
	
	#produce bam with all reads
	#at this stage we also add read groups
	call("bwa samse -r \'@RG\t" + RG + "\t\' " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)
	#print("bwa samse " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam")
	#call("bwa samse " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)
		
	#sort this bam
	call("samtools sort " + sample + ".bam " + sample + "_sort",shell=True)
	
	#remove duplicates from the sorted bam
	call("samtools rmdup -s " + sample + "_sort.bam " + sample + "_rmdup.bam", shell=True)
	
	#remove the "sorted with duplicates" bam
	call("rm " + sample + "_sort.bam",shell=True)

	#make a copy of the samtools flagstat
	call("samtools flagstat " + sample + "_rmdup.bam > " + sample + "_flagstat.txt",shell=True)

	#remove unaligned reads from this bam
	call("samtools view -b -F 4 " + sample + "_rmdup.bam > " + sample + "_rmdup_only_aligned.bam",shell=True)

	#produce q30 bams
	call("samtools view -b -q30 " + sample + ".bam >" + sample + "_q30.bam",shell=True)

	#sort this bam
       	call("samtools sort " + sample + "_q30.bam " + sample + "_q30_sort",shell=True)

       	#remove duplicates from the sorted bam
       	call("samtools rmdup -s " + sample + "_q30_sort.bam " + sample + "_q30_rmdup.bam", shell=True)

       	#remove the "sorted with duplicates" bam
       	call("rm " + sample + "_q30_sort.bam",shell=True)

       	#make a copy of the samtools flagstat
       	call("samtools flagstat " + sample + "_q30_rmdup.bam > " + sample + "_q30_flagstat.txt",shell=True)

	#get number of reads aligned without rmdup
	raw_reads_aligned = subprocess.check_output("samtools flagstat " + sample + ".bam |  grep 'mapped (' | cut -f1 -d' '",shell=True)
	summary.append(raw_reads_aligned)

	#get reads that aligned following rmdup
	rmdup_reads_aligned = subprocess.check_output("more " + sample + "_flagstat.txt | grep 'mapped (' | cut -f1 -d' '",shell=True)
	summary.append(rmdup_reads_aligned)

  	#capture the alignment percentage of the flagstat file, both no q and q30
	raw_alignment = subprocess.check_output("more " + sample + "_flagstat.txt | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	summary.append(raw_alignment)

	#get q30 reads aligned
	q30_reads_aligned = subprocess.check_output("more " + sample + "_q30_flagstat.txt | grep 'mapped (' | cut -f1 -d' '",shell=True)
       	summary.append(q30_reads_aligned)

	#q30_percent_aligned = subprocess.check_output("more " + sample + "_q30_flagstat.txt | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	fixed_percentage = (int(q30_reads_aligned)) / int(trimmed_read_number)
	summary.append(fixed_percentage)
	
	#remove unaligned reads from this bam
       	call("samtools view -b -F 4 " + sample + "_q30_rmdup.bam > " + sample + "_q30_rmdup_only_aligned.bam",shell=True)
	
	#index the q30 bam
	call("samtools index "+ sample + "_q30_rmdup_only_aligned.bam",shell=True)
	
	#get idx stats
	call("samtools idxstats "+ sample + "_q30_rmdup_only_aligned.bam > "  + sample + ".idx",shell=True)
	
	#clean up files
	call("gzip "+ sample + "*",shell=True)
	
	#going to make an output directory for each sample
	#then move all produced files to this directory
	call("mkdir " + out_dir + sample,shell=True)

	call("mv *.bam* *.idx* *trimmed.fastq *.txt* " + out_dir + sample,shell=True)
	
	#add sample summary to the masterlist
	fixed_summary = []
	for entry in summary:
	
		entry = str(entry).rstrip("\n")
		fixed_summary.append(entry)
	
	print fixed_summary

	master_list.append(fixed_summary)

#Move all cutadapt logs to a cutadapt log directory - if there is no such dir, create it
call("mkdir " + out_dir + "cutadapt_logs/", shell=True)
call("mv *.log " + out_dir + "cutadapt_logs/", shell=True)

#remove all .sai files
call("rm *sai*",shell=True)

print "Here a list of the files which failed:"

print failures


output_summary = miseq_date + "_summary.table"
	
#print summary stats
with open(output_summary, "w") as f:

	writer = csv.writer(f, delimiter='\t', lineterminator='\n')
	writer.writerows(master_list)

call("mv " + output_summary + " " + out_dir,shell=True)