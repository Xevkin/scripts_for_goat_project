#!/usr/bin/python


#also need to supply a date for the Miseq run - will automatically make results file
#python script <date_of_miseq> <meyer> <option> <read group file>

#read group file needs to be in the follow format: 
#Sample_name\tRGs_to_add (in following format- ID:X\tSM:X\tPL:X\tLB:X)

#import csv module to easily write the list of list used as a .csv file
import csv

#apparently this is a better way to make system calls using python, rather than "os.system"
import subprocess 
from subprocess import call

#need to import sys anyways to access input file (a list)
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

#variables will be initialized here so they can be modified by options 


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
option = sys.argv[3].rstrip("\n").lower()

if (option == "sheep"):
	print "Sheep alignment selected."
	reference = "~/goat/miseq/data/reference_genomes/sheep_oviAri3/oviAri3.fa" 

#if option is mit, then several changes need to occur

cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 "

fastq_screen = "~/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --outdir ./"

#now change variables if the mitochondrial option has been selected

if (option == "mit"):

	print "Mitochondrial alignment selected."

	cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 25 "

	reference = "~/goat/miseq/data/mit_genomes_fastq_screen/goat_mit/goat_mit.fasta"

	fastq_screen = "~/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --conf ~/goat/miseq/src/fastq_screen_v0.4.4/capture_config/capture.conf --outdir ./"

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
	call("gunzip " + current_file, shell=True)
	current_file = current_file.split(".")[0] + ".fastq"

	#rename the file
	#also save the current sample as a variable to be used later
	split_file = current_file.split(".")[0].split("_")
	print split_file[0]
	call("mv " + current_file + " " + split_file[0] + ".fastq",shell=True)
	
	current_file = split_file[0] + ".fastq"
	
	sample = split_file[0].rstrip("\n")
	
	print "Current samples is: " + sample
	
	
	#initialize variable to carry summary statistics
	summary = []
	summary.append(sample)
	
	#Get number of lines (and from that reads - divide by four) from raw fastq
	unzipped_fastq = current_file
	trimmed_fastq = sample + "_trimmed" + ".fastq" 
	
	cmd = "wc -l " + unzipped_fastq + " | cut -f1 -d' '" 
	summary.append(subprocess.check_output(cmd,shell=True))
	
	call(cmd,shell=True)
	
	#raw reads
	raw_read_number = int(subprocess.check_output(cmd,shell=True)) / 4
	summary.append(raw_read_number)
	
	#cut raw fastq files
	call(cut_adapt + unzipped_fastq + " > " + trimmed_fastq + " 2> " + trimmed_fastq + ".log", shell=True)
	
	#grab summary statistics of trimmed file
	cmd = "wc -l " + trimmed_fastq + "| cut -f1 -d' '"
       	summary.append(subprocess.check_output(cmd,shell=True))
       	trimmed_read_number = int(subprocess.check_output(cmd,shell=True)) / 4
       	summary.append(str(trimmed_read_number))


	#run fastqc on trimmed fastq file
	#first we want to create an output directory if there is none to begin with
	call("mkdir " +  out_dir + "fastqc/", shell=True)
	
	call("fastqc " + unzipped_fastq + " -o " + out_dir + "fastqc/", shell=True)
	call("fastqc " + trimmed_fastq + " -o " + out_dir + "fastqc/", shell=True)

	#run fastq screen
	call("mkdir " + out_dir + "fastq_screen/",shell=True)
	call("mkdir " + out_dir + "fastq_screen/" + sample, shell=True)
	
	call(fastq_screen + sample + " " + trimmed_fastq + " --outdir " + out_dir + "fastq_screen/" + sample, shell=True)
	
	#at this stage we have our fastq files with adaptors trimmed
	#we can now move on to the next step: alignment
	#going to align to CHIR1.0, as that what was used for AdaptMap
	print(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai")
	call(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai",shell=True)	
	
	#Obtain the appropriate read group from the supplied read group file
	with open(RG_file) as file:

		for line in file:

			split_line = line.split("\t")
			
			if (sample == split_line[0]):

				RG = split_line[1].rstrip("\n")
				print RG
				#check if RG is an empty string
				if not RG:
		
					print "No RGs were detected for this sample - please check sample names in fastq files and in RG file agree" 
					failures.append(current_file)
	
				else:
		
					print "Reads groups being used are:"
                                	print RG

	#produce bam with all reads
	#at this stage we also add read groups
	call("bwa samse -r \'@RG\t" + RG + "\t\' " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)
			
	#Run mapDamage on  bam to check degradation patterns 
	#mapdmg_out = out_dir + "mapDamage/" + sample
	#call("mkdir "+ mapdmg_out,shell=True)	
	#call("mapDamage map -i " + sample + ".bam -d " + mapdmg_out +" -r " +reference +" -t sample",shell=True)
	
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

	#produce q25 bams
        call("samtools view -b -q25 " + sample + "_rmdup_only_aligned.bam >" + sample + "_q25.bam",shell=True)

        #sort this bam
        call("samtools sort " + sample + "_q25.bam " + sample + "_q25_sort",shell=True)

        #remove duplicates from the sorted bam
        call("samtools rmdup -s " + sample + "_q25_sort.bam " + sample + "_q25_rmdup.bam", shell=True)

        #remove the "sorted with duplicates" bam
        call("rm " + sample + "_q25_sort.bam",shell=True)
        
	#produce q30 bams
	call("samtools view -b -q30 " + sample + "_rmdup_only_aligned.bam >" + sample + "_q30.bam",shell=True)

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
	fixed_percentage = ((float(q30_reads_aligned)) * 100)/ float(trimmed_read_number)
	summary.append(fixed_percentage)
		
	#index the q30 bam
	call("samtools index "+ sample + "_q30_rmdup.bam",shell=True)
	
	#get idx stats
	call("samtools idxstats "+ sample + "_q30_rmdup.bam > "  + sample + ".idx",shell=True)
	
	#clean up files
	call("gzip "+ sample + "*",shell=True)
	
	#going to make an output directory for each sample
	#then move all produced files to this directory
	call("mkdir " + out_dir + sample,shell=True)
	call("mkdir trimmed_fastq_files_and_logs",shell=True)
	call("mv *trimmed* trimmed_fastq_files_and_logs/",shell=True)
	
	#add sample summary to the masterlist
	fixed_summary = []
	for entry in summary:
	
		entry = str(entry).rstrip("\n")
		fixed_summary.append(entry)
	
	print fixed_summary

	master_list.append(fixed_summary)

	call("mv *.bam* *.idx* *flagstat* " + out_dir + sample,shell=True)

#remove all .sai files
call("rm *sai*",shell=True)
call("mv trimmed_fastq_files_and_logs/ " + out_dir,shell=True)
print "Here a list of the files which failed:"

print failures


output_summary = miseq_date + "_summary.table"
	
#print summary stats
with open(output_summary, "w") as f:

	writer = csv.writer(f, delimiter='\t', lineterminator='\n')
	writer.writerows(master_list)

call("mv " + output_summary + " " + out_dir,shell=True)
