#!/bin/python


#this script takes a list of .fastq.gz files, and does intial processing of the files
#also need to supply a date for the Miseq run - will automatically make results file
#python script <input.list> <date_of_miseq>

#apparently this is a better way to ake system calls using python, rather than "os.system"
from subprocess import call

#need to import sys anyways to access intput file (a list)
import sys

#we will use this later to check if the input files actually exist
import os.path

#reference genomes
goat_ref = "~/goat/data/reference_genomes/goat_CHIR1_0/goat_CHIR1_0.fasta"

#Prepare output directory

miseq_date = sys.argv[2]

out_dir = "~/goat/results/" + miseq_date + "/"

call("mkdir " + out_dir, shell=True)

#allow meyer option to be used
meyer_input = sys.argv[3].rstrip("\n").lower()

alignment_option = "bwa aln -l 1000 "  

if meyer_input == "meyer":

	alignment_option = "bwa aln -l 1000 -n 0.01 -o 2 " 

	

#sys.argv[1] refers to the list input
#create list that will carry any lines in the file which fail at any stage
failures = []
f = open(sys.argv[1])
#cycle through each line in the input file
for lines in f:
	
	current_file = lines.strip()
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
	
	sample = split_file[0]
	print sample
	unzipped_fastq = sample + "." + split_file[1]

	trimmed_fastq = sample + "_trimmed" + "." + split_file[1]
	
	call("cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 " + unzipped_fastq + " > " + trimmed_fastq + " 2> " + trimmed_fastq + ".log", shell=True)
	#print "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 " + unzipped_fastq + " > " + trimmed_fastq
	
	#rezip original fastq files
	
	call("gzip " + unzipped_fastq, shell=True)
	
	#run fastqc on trimmed fastq file
	#first we want to create an output directory if there is none to begin with
	#output_dir = "fastqc_output/"
	#if not os.path.exists(output_dir):
	#	os.makedirs(output_dir)
	
	#print "fastqc " + trimmed_fastq + " -o fastqc_output/ "	
	#if not os.path.exists(out_dir + "fastqc/"):
	#	os.makedirs(out_dir + "fastqc/")
	call("mkdir " +  out_dir + "fastqc/", shell=True)
	call("fastqc " + trimmed_fastq + " -o " + out_dir + "fastqc/", shell=True)

	#run fastq_screen on each fastq, making a directory for each output
	#if not os.path.exists(out_dir + "fastq_screen/"):
	#	os.makedirs(out_dir + "fastq_screen/")
	
	#if not os.path.exists("./" + sample):
	#	os.makedirs("./" + sample)
	
	call("mkdir " + out_dir + "fastq_screen/" + sample, shell=True)
	call("mkdir ./" + sample, shell=True)
	
	call("~/goat/src/fastq_screen_v0.4.4/fastq_screen --aligner bowtie --outdir ./" + sample + " " + trimmed_fastq, shell=True)
	print "Output of fastq_screen:"
	call("ls ./" + sample,shell=True)
	#call("rsync --remove-source-files ./" + sample + " " + out_dir + "fastq_screen/" + sample,shell=True)
	call("mv ./" + sample + " " +  out_dir + "fastq_screen/",shell=True)
	call("rm -r ./" + sample,shell=True)

	#at this stage we have our fastq files with adaptors trimmed
	#we can now move on to the next step: alignment
	#going to align to CHIR1.0, as that what was used for AdaptMap
	print(alignment_option + goat_ref + " " + trimmed_fastq + " > " + sample + ".sai")
	call(alignment_option + goat_ref + " " + trimmed_fastq + " > " + sample + ".sai",shell=True)	
	
	#produce bam with all reads
	print("bwa samse " + goat_ref + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam")
	call("bwa samse " + goat_ref + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam",shell=True)
	#add RGs
	#call("java -jar /research/picard-tools-1.119/AddOrReplaceReadGroups.jar RGID=MG2b_MiSeq RGLB=" + sample + "RGPL=ILLUMINA_Miseq RGPU=index1 RGSM=MG2 INPUT=" + sample + "OUTPUT=" + sample + "_RG.bam",shell=True)
	
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
	call("samtools view -b -q30 " + sample + ".bam" > sample + "_q30.bam",shell=True)
	#process as normal
	#sort this bam
        call("samtools sort " + sample + "_q30.bam " + sample + "_q30_sort",shell=True)
        #remove duplicates from the sorted bam
        call("samtools rmdup -s " + sample + "_q30_sort.bam " + sample + "_q30_rmdup.bam", shell=True)
        #remove the "sorted with duplicates" bam
        call("rm " + sample + "_q30_sort.bam",shell=True)
        #make a copy of the samtools flagstat
        call("samtools flagstat " + sample + "_q30_rmdup.bam > " + sample + "_q30_flagstat.txt",shell=True)
        #remove unaligned reads from this bam
        call("samtools view -b -F 4 " + sample + "_30_rmdup.bam > " + sample + "_rmdup_only_aligned.bam",shell=True)

	#index the q30 bam
	call("samtools index "+ sample + "_q30_rmdup_only_aligned.bam",shell=True)
	
	#get idx stats
	call("samtools idxstats "+ sample + "_q30_rmdup_only_aligned.bam",shell=True)

	#gzip both files
	call("gzip " + sample + "_rmdup_only_aligned.bam",shell=True)
	call("gzip " + sample + "_rmdup.bam",shell=True)
	call("gzip " + sample + ".bam",shell=True)
        call("gzip " + sample + "_q30_rmdup_only_aligned.bam",shell=True)
        call("gzip " + sample + "_q30_rmdup.bam",shell=True)
        call("gzip " + sample + "_q30.bam",shell=True)

#Move all cutadapt logs to a cutadapt log directory - if there is no such dir, create it
call("mkdir " + out_dir + "cutadapt_logs/", shell=True)
call("mv *.log " + out_dir + "cutadapt_logs/", shell=True)

#remove all .sai files
call("rm *sai*",shell=True)
print "Here a list of the files which failed:"

print failures
			


