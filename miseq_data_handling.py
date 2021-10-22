#!/usr/bin/env python
'''
Script is run in a directory with bam files to be aligned to the goat genome
Bams will be trimmed, quality metrics obtained, screened using fastq screen

Also need to supply a date for the Miseq run - will automatically make results file
python script <date_of_miseq> <meyer> <species> <mit> <trim> <fastqc> <fastq_screen> <directory in which to place output dir/> <PE>


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

	"goat" : "/home/kdaly/ARS1/ARS1.fa",

	"sheep" : "/home/kdaly/st1/miseq/reference_genomes/oviAri3.fa",

	"cow" : "/Reference_Genomes/bos_2019/ARS-UCD1.2_Btau5.0.1Y.fa",

	"dog" : "/eno/reference_genomes/dog_canFam3/canFam3.fa",

	"horse" : "/Reference_Genomes/For_Fastq_Screen/horse.fa",

	"tur" : " ",

	"deer" : "/home/kdaly/st1/deer/deer_genome/CerEla1-0_mod.fa"
}

mitochondrial_genomes = {

	"goat" : "/eno/miseq/goat/miseq/data/mit_reference_genomes/goat/goat_mit_revised.fa",

	"sheep" : "/eno/miseq/data/mit_genomes_fastq_screen/sheep_mit/sheep_mit.fasta",

	"tur" : "/eno/miseq/goat/miseq/data/mit_reference_genomes/tur/west_caucus_tur.fasta"

}

def main(date_of_miseq, meyer, species, mit, fastq_screen,  output_dir, trim, fastqc, pair):

	#run the set up function.#set up will create some output directories
	#and return variables that will be used in the rest of the script

	files, reference, out_dir, cut_adapt, alignment_option, master_list, sample_list, fastq_screen_option,adaptor_removal = set_up(date_of_miseq, meyer, species, mit, output_dir, trim, pair)

	#sample is the file root
	#trim fastq files and produce fastqc files
	#the masterlist will change each time so it needs be equated to the function

	if (trim == "yes"):

		#make an output directory for fastqc
		call("mkdir " +  out_dir + "fastqc/", shell=True)

		for sample in sample_list:

			trim_fastq(sample, cut_adapt, out_dir, fastqc,adaptor_removal, pair)

	#run fastq screen on the samples if input variable for fastq_screen is "yes"
	if (fastq_screen.lower() == "yes"):

		map(lambda sample: run_fastq_screen(sample, out_dir, fastq_screen_option, pair), sample_list)

	#at this stage we have our fastq files with adaptors trimmed, fastqc and fastq screen run
	#we can now move on to the next step: alignment

	#going to align to ARS1

	map(lambda sample : align(sample,  alignment_option, reference, pair), sample_list)

	#testing a function here, to process a bam to a q30 version
	map(process_bam, sample_list)\

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

		writer = csv.writer(f, delimiter=' ', lineterminator='\n')
		writer.writerows(master_list)

	call("wc -l " + output_summary,shell=True)

	number_of_samples = (int((subprocess.check_output("wc -l " + output_summary,shell=True).split(" ")[0])) - 1)

	call("head -n1 " + output_summary + "> header.txt ",shell=True)

	call("tail -n " + str(number_of_samples) + " " + output_summary + " | sort | cat header.txt - > tmp; mv tmp " + output_summary + ";rm header.txt tmp",shell=True)

	call("mv " + output_summary + " " + out_dir,shell=True)


def set_up(date_of_miseq, meyer, species, mit,  output_dir, trim, pair):

	adaptor_removal = ""

	#clean up some file names
	call("for i in $(ls *fastq.gz); do I=`echo $i | sed -e \"s/_001././g\" | sed -e \"s/_/-/g\"`; mv $i $I; done",shell=True)

	files = []
	#if not trim, take all "trimmed.fastq" files:
	if (trim == "no"):

		if  (pair == "pair") or (pair == "yes"):

			 files = [file for file in os.listdir(".") if file.endswith("R1.trimmed.fastq.gz")]

		else:

			files = [file for file in os.listdir(".") if file.endswith("trimmed.fastq.gz")]
		#print the trimmed fastq files in current directory
		print "Trimmed fastq files in the curent directory:"
		print map(lambda x : x ,files)

	else:
		#take all .fastq.gz files in current directory; print them

		if  (pair == "pair") or (pair == "yes"):

			files = [file for file in os.listdir(".") if file.endswith("R1.fastq.gz")]

		else:

			files = [file for file in os.listdir(".") if file.endswith(".fastq.gz")]

		print "fastq.gz files in current directory:"

		print map(lambda x : x ,files)

	#variables will be initialized here so they can be modified by options


	if (mit != "yes"):

		#reference genome to be used
		reference = nuclear_genomes[species]

		#define default cut_adapt and fastq screen

	        cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 "

        	fastq_screen = "/Software/fastq_screen --aligner bowtie --force --outdir ./"

		adaptor_removal = "/home/kdaly/programs/adapterremoval-2.3.1/build/AdapterRemoval --threads 2 --collapse --minadapteroverlap 1 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minlength 30 --gzip --trimns --trimqualities "

	else:

		print "Mitochondrial alignment selected"

		reference = mitochondrial_genomes[species]

		cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 25 "

                fastq_screen = "/eno/miseq/goat/miseq/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --conf /eno/miseq/src/fastq_screen_v0.4.4/capture_config/capture.conf --outdir ./"


	print "Species selected is " + species

        print "Path to reference genome is " + reference


	#Prepare output directory

	out_dir = output_dir + date_of_miseq +  "/"

	call("mkdir " + out_dir, shell=True)

	#allow meyer option to be used
	meyer_input = meyer.rstrip("\n").lower()

	alignment_option = "bwa aln -t 8 -l 1000 "

	if (meyer_input == "meyer" or meyer_input == "yes"):

		print "Meyer option selected."

		alignment_option = "bwa aln -t 8 -l 1000 -n 0.01 -o 2 "

	#initialize a masterlist that will carry summary stats of each sample
	master_list = [["Sample", "read_count_raw", "trimmed_read_count","raw_reads_aligned", "raw %age endogenous", "rmdup_reads_remaining","rmdup_reads_aligned" ,"rmdup_alignment_percent", "reads_aligned_q30", "percentage_reads_aligned_q30"]]

	sample_list = []
	#cycle through each line in the input file, gunzip
	for file in files:

		#handle trimmed.fastq files differently

		if (trim == "no"):

			current_file = "_".join(file.split("_")[0:-1])
			sample_list.append(current_file.rstrip("\n"))

		else:

			#call("gunzip " + file, shell=True)
			current_file = file.split(".")[0] + ".fastq.gz"

			#rename the file
			#also save the current sample as a variable to be used later
			split_file = current_file.split(".")[0].split("_")

			call("mv " + current_file + " " + split_file[0] + ".fastq.gz",shell=True)

			#current_file = split_file[0] + ".fastq"

			sample_list.append(split_file[0].rstrip("\n"))

	for i in sample_list:

		master_list.append([i])

	return files, reference, out_dir, cut_adapt, alignment_option, master_list, sample_list, fastq_screen, adaptor_removal


def trim_fastq(current_sample, cut_adapt, out_dir ,fastqc, adaptor_removal, pair):

	print "Current sample is: " + current_sample

	fastq = current_sample + ".fastq.gz"

	#Get number of lines (and from that reads - divide by four) from raw fastq
	trimmed_fastq = current_sample + "_trimmed" + ".fastq.gz"

	cmd = "gunzip -c " + fastq + " | wc -l |cut -f1 -d' '"

	#cut raw fastq files
	if (pair == "pair") or (pair == "yes"):

		cmd = adaptor_removal + " --file1 " + fastq + " --file2 "  + fastq.replace("_R1","_R2").replace("-R1","-R2") + " --basename " + current_sample.replace("-R1","").replace("_R1","") + " 2> " + current_sample + ".trim.log && mv " + current_sample.replace("-R1","").replace("_R1","") + ".collapsed.gz " + current_sample.replace("-R1","").replace("_R1","") + "_trimmed.collapsed.fastq.gz  ; mv " +  current_sample + "_trimmed.pair1.truncated.gz " + current_sample + "_r1_trimmed.fastq.gz ; mv " + current_sample + "_trimmed.pair2.truncated.gz " +  current_sample + "_r2_trimmed.fastq.gz ; mv " +  current_sample + "_trimmed.singleton.truncated.gz " + current_sample + "_mate-discard_trimmed.fastq.gz "

		print(cmd)

		call(cmd,shell=True)

	else:

		call(cut_adapt + fastq + " -o " + current_sample + "_trimmed.fastq.gz  > " + trimmed_fastq + ".log", shell=True)

       	#run fastqc on both the un/trimmed fastq files
	#first we want to create an output directory if there is none to begin with
	if (fastqc == "yes"):

		call("~/programs/FastQC/fastqc " + fastq + " -o " + out_dir + "fastqc/", shell=True)

		if  (pair == "pair") or (pair == "yes"):

			call("~/programs/FastQC/fastqc " + current_sample  + "_trimmed.collapsed.gz -o " + out_dir + "fastqc/", shell=True)


def run_fastq_screen(current_sample, out_dir, fastq_screen_option,pair):

	call("mkdir " + out_dir + "fastq_screen/",shell=True)

	call("mkdir " + out_dir + "fastq_screen/" + current_sample, shell=True)

	if (pair == "pair") or (pair == "yes"):

		 call(fastq_screen_option + current_sample + " " + current_sample + "_trimmed.collapsed.fastq.gz --outdir " + out_dir + "fastq_screen/" + current_sample, shell=True)

	else:

		call(fastq_screen_option + current_sample + " " + current_sample + "_trimmed.fastq.gz --outdir " + out_dir + "fastq_screen/" + current_sample, shell=True)


def align(sample, alignment_option, reference, pair):

	if (pair == "pair") or (pair == "yes"):

		trimmed_fastq = sample.replace("-r1","").replace("-R1","") + "_trimmed.collapsed.fastq.gz"

	else:

		trimmed_fastq = sample + "_trimmed.fastq.gz"

	print(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai")
	call(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai",shell=True)

	print sample
	call("bwa samse "  + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)

	call("samtools flagstat " + sample + ".bam > " +  sample + ".flagstat",shell=True)


def process_bam(sample_name):

	#sort this bam
	call("~/programs/samtools-HL/samtools sort " + sample_name + ".bam " + sample_name + "_sort",shell=True)

	#remove duplicates from the sorted bam
	call("samtools rmdup -s " + sample_name + "_sort.bam " + sample_name + "_rmdup.bam", shell=True)

	#remove the "sorted with duplicates" bam
	call("rm " + sample_name + "_sort.bam",shell=True)

	#make a copy of the samtools flagstat
	call("samtools flagstat " + sample_name + "_rmdup.bam > " + sample_name + "_rmdup.flagstat",shell=True)

	#produce q30 bams
	call("samtools view -b -F 4 -q30 " + sample_name + "_rmdup.bam >" + sample_name + "_rmdup_q30.bam",shell=True)

      	#make a copy of the samtools flagstat
      	call("samtools flagstat " + sample_name + "_rmdup_q30.bam > " + sample_name + "_rmdup_q30.flagstat",shell=True)

      	#index the q30 bam
	call("samtools index "+ sample_name + "_rmdup_q30.bam",shell=True)

	#get idx stats
	call("samtools idxstats "+ sample_name + "_rmdup_q30.bam > "  + sample_name + ".idx",shell=True)


def get_summary_info(master_list, current_sample):

	print "Initial master list is:"
	print master_list
	to_add = []

	fastq = current_sample + ".fastq.gz"

	trimmed_fastq = current_sample + "_trimmed.fastq.gz"

	cmd = "gunzip -c " + fastq + " | wc -l | cut -f1 -d' '" 

	#raw reads
	file_length = subprocess.check_output(cmd,shell=True)
	raw_read_number = int(file_length) / 4

	to_add.append(raw_read_number)

	#grab summary statistics of trimmed file
	cmd = "gunzip -c  " + trimmed_fastq + " | wc -l| cut -f1 -d' '"
       	file_length = subprocess.check_output(cmd,shell=True)
	trimmed_read_number = int(file_length) / 4

	to_add.append(str(trimmed_read_number).rstrip("\n"))
       	#get number of reads aligned without rmdup
	raw_reads_aligned = subprocess.check_output("samtools flagstat " + current_sample + ".bam |  grep 'mapped (' | cut -f1 -d' '",shell=True)
	to_add.append(raw_reads_aligned.rstrip("\n"))

	#get %age raw alignment
	raw_alignment_percentage = ((float(raw_reads_aligned)) * 100)/ float(trimmed_read_number)
	to_add.append(str(raw_alignment_percentage).rstrip("\n"))


	rmdup_reads_remaining = subprocess.check_output("more " + current_sample + "_rmdup.flagstat | head -n1 | cut -f1 -d' '",shell=True) 
	to_add.append(rmdup_reads_remaining.rstrip("\n"))

	#get reads that aligned following rmdup
	rmdup_reads_aligned = subprocess.check_output("more " + current_sample + "_rmdup.flagstat | grep 'mapped (' | cut -f1 -d' '",shell=True)
	to_add.append(rmdup_reads_aligned.rstrip("\n"))

  	#capture the alignment percentage of the flagstat file, both no q and q30
	raw_alignment = subprocess.check_output("more " + current_sample + ".flagstat | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	to_add.append(raw_alignment.rstrip("\n"))

	#get q30 reads aligned
	q30_reads_aligned = subprocess.check_output("more " + current_sample + "_rmdup_q30.flagstat | grep 'mapped (' | cut -f1 -d' '",shell=True)
       	to_add.append(q30_reads_aligned.rstrip("\n"))

	#q30_percent_aligned = subprocess.check_output("more " + sample + "_q30_flagstat.txt | grep 'mapped (' | cut -f5 -d' ' | cut -f1 -d'%' | sed 's/(//'", shell=True)
	fixed_percentage = str(((float(q30_reads_aligned)) * 100)/ float(rmdup_reads_remaining))
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
	trim = sys.argv[5]
	fastqc = sys.argv[6]
	fastq_screen = sys.argv[7]
	output_dir = sys.argv[8]
	pair = sys.argv[9]

except IndexError:
	print "Incorrect number of variables have been provided"
	print "Input variab\es are date_of_miseq, meyer, species, mit, trim, fastqc, fastq_screen, output directory, and paired"
	print "Exiting program..."
	sys.exit()


main(date_of_miseq, meyer, species, mit, fastq_screen,  output_dir, trim, fastqc, pair)
