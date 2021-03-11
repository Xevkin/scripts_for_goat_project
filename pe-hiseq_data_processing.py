#!/usr/bin/env python
'''
Kevin Daly 2020 - for PE data

Script is run in a directory with bam files to be aligned to either goat, wild goat or sheep genome
Reads will be trimmed, aligned, read groups added, duplicates removed, indel realignment performed, and then softclipped

Also need to supply a date for the Hiseq run - will automatically make results file
python script <date_of_hiseq> <meyer> <species> <mit_file> <trim> <skip mit alignment> <align> <process_individual_bams> <merge+process> <clip> <read group file> <directory in which to place output dir/>

mit_file should be tab deliminated, col1 with name of reference, col2 with path

read group file needs to be in the follow format (\t are actual tab characters)
FASTQ_FILE_r1.GZ\t@RG+ID:X+SM:X+PL:X+LB:X\tSAMPLE_NAME
ID should be in the format <sample_name>-<macrogen_index_number>-<lane_number>-<hiseq_number>
LB should refer to PCR: <sample>-<lab index>-<macrogen-index>-<PCR_number>

fastq files should be in the format <sample>-<PCR number and letter, if any>_<underscores and numbers>, followed by "_r1"

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

	"ARS1" : "/home/kdaly/ARS1/ARS1.fa",

	"CHIR1" : "/home/kdaly/goat_CHIR1_0/goat_CHIR1_0.fasta",

	"Sheep" : "/lascelles/sheep_OarRambouillet/OarRambouillet.fa",

	"Deer" : "/home/kdaly/st1/deer/deer_genome/CerEla1-0_mod.fa"
}


def main(date_of_hiseq, meyer, threads, species, mit, skip_mit_align, trim, align, process, merge, clip, RG_file, output_dir):

	#run the set up function.#set up will create some output directories
	#and return variables that will be used in the rest of the script

	files, reference, mit_references, out_dir, adapter_removal, alignment_option, fastq_list = set_up(date_of_hiseq, meyer, threads,species, mit, RG_file, output_dir, trim)

	#make a folder for flagstat files
	call("mkdir flagstat_files",shell=True)

	#sample is the file root
	#trim fastq files and produce fastqc files - if we want to
	#the masterlist will change each time so it needs be equated to the function

	if (trim == "yes" or trim == "trim"):

	        #remove the trimming command script
	        #call("rm trim.sh",shell=True)

		#create the trimming command file
		for fastq in fastq_list:

			prepare_trim_fastq(fastq, adapter_removal, out_dir)

		#now actually trim in parallel, and remove the trim.sh file afterwards
		call("parallel -j " + threads + " -a trim.sh; rm trim.sh",shell=True)

	#at this stage we have our fastq files with adaptors trimmed
	#We can now move on to the next step: alignment

	#going to align to CHIR1.0, as that what was used for AdaptMap
	#this step will align each fastq and produce raw bams with RGs
	#they should have the same stem of the initial bam


	if (mit_references != "no" ):

		#align to every mitochondrial reference provided
		for mitochondrial_reference in mit_references:

			print "Doing mit alignment to " + mitochondrial_reference[0]

			if (skip_mit_align=="no"):

				map (lambda fastq : align_process_mit(fastq, RG_file, alignment_option, mitochondrial_reference, trim), fastq_list)

				merge_and_process_mit(RG_file, mitochondrial_reference, trim)

		clean_up_mit(mit,out_dir)

	print "alignment option is " + align

	if (align == "yes" or align == "align"):

		print fastq_list

		#remove the samse.sh file
		call("rm samse.sh",shell=True)

		map (lambda fastq : align_bam(fastq, RG_file, alignment_option, reference, trim, species), fastq_list)

		#run the samse.sh command file, then remove it
		call("parallel -j " +threads + " -a samse.sh; rm samse.sh",shell=True)

		#remove the sampe.sh file
		call("rm sampe.sh",shell=True)

		map (lambda fastq : align_bam_pe(fastq, RG_file, alignment_option, reference, trim, species), fastq_list)

		#run the sampe.sh command file, then remove it
                call("parallel -j " +threads + " -a sampe.sh; rm sampe.sh",shell=True)


	#sleep "$i"s & pids+=( $! ); done; wait "${pids[@]}
	#Remove the sai files that we don't care about
	call("rm *sai",shell=True)

	#sort and remove duplicates from each bam
	if (process == "yes" or process == "process"):

		for sample in fastq_list:

			process_bam(sample,species)

		#map(process_bam, fastq_list, species)

        	#do all markdup at same time
		call("echo Removing duplicates",shell=True)

		call("PIDS_list="";for i in $(ls *sort_q20.bam | grep -v \"_pe_\" | rev | cut -f3- -d'_' | rev ); do echo samtools markdups -r \"$i\"_sort_q20.bam \"$i\"_q20_markdup.bam \"2>\" \"$i\"_q20_markdup.log;  samtools markdups -r \"$i\"_sort_q20.bam \"$i\"_q20_markdup.bam 2> \"$i\"_q20_markdup.log & PIDS_list=`echo $PIDS_list $!`; done; for i in $(ls *resort.bam | grep \"_pe_\" | rev | cut -f3- -d'_' | rev ); do echo samtools markdup -r -@ 6 \"$i\"_sort_q20_fix_resort.bam \"$i\"_q20_markdup.bam \"2>\" \"$i\"_q20_markdup.log; samtools markdup  \"$i\"_sort_q20_fix_resort.bam \"$i\"_q20_markdup.bam 2> \"$i\"_q20_markdup.log & PIDS_list=`echo $PIDS_list $!`; done; for pid in $PIDS_list; do wait $pid; done",shell=True)

		call("for i in $(ls *markdup.bam | cut -f 1 -d'.'); do samtools flagstat -@ 12 ${i}.bam > ${i}.flagstat; samtools view -@ 12 -b -F 4 $i.bam > tmp.bam; mv tmp.bam $i.bam ;done; rm tmp.bam",shell=True)

		call("for i in $(ls *markdup.bam | cut -f 1 -d'.'); do samtools flagstat $i.bam > $i.flagstat; done",shell=True)

		call("for i in $(ls *_pe_*markdup.bam | cut -f 1 -d'.'); do samtools view -@ 12 -f 2 -b $i.bam > ${i}_ppair.bam && samtools flagstat -@ 12 ${i}_ppair.bam > ${i}_ppair.flagstat ; done",)

	#add an option here to kill the script if you do not want merging to occur
	if (merge == "no"):

		call("mkdir auxillary_files",shell=True)

		call("mv " + RG_file + " auxillary_files",shell=True)

		call("mv auxillary_files " + out_dir,shell=True)

		sys.exit("Script is terminated as no merging was desired")

    	#now merge the lanes for each sample
	#NOTE: if all the options up to process are "no", then this is where the script will start
	#expects markdup bams for each index in each lane, ungziped
	#output will be merged, no markdup
	merged_bam_list = merge_lanes_and_sample(RG_file, trim, species)

	merged_markdup_bam_list = []

	call("rm markdup.sh",shell=True)

	for bam in merged_bam_list:

		sample =  bam.split(".")[0]

		call("echo samtools markdups -r " + bam + " " + sample + "_markdup.bam \">\"  " + sample + "_markdup.log >> markdup.sh",shell=True)

		merged_markdup_bam_list.append(sample + "_markdup.bam")

	call("parallel -a markdup.sh -j " + threads, shell=True)

	realigned_bam_list = []

	#indel realignment
	for i in merged_markdup_bam_list:

		print "Indel realignment on sample " + i

		realigned_bam_list.append(indel_realignment(i,reference))

	#now run the script to process the realigned bams
	#this function will invoke a function that rescales the bams
	#duplicates will also be removed at this point
	for i in realigned_bam_list:

		process_realigned_bams(i,reference,clip,output_dir,species)

	clean_up(out_dir)


def set_up(date_of_hiseq, meyer, threads ,species, mit, RG_file, output_dir, trim):
	#take all .fastq.gz files in current directory; print them
	files = []

	if (trim == "yes"):

		files = [file for file in os.listdir(".") if file.endswith("r1.fastq.gz")]

		print "r1.fastq.gz files in current directory:"

		print map(lambda x : x ,files)

	#however, if trim option is not yes, then we use fastq files
	if (trim != "yes"):

		files = [file for file in os.listdir(".") if file.endswith("r1_trimmed.fastq.gz")]

        	print "trimmed fastq files in current directory:"

	        print map(lambda x : x ,files)

	#variables will be initialized here so they can be modified by options

	#reference genome to be used
	reference = nuclear_genomes[species]

	print "Species selected is " + species

	print "Path to reference genome is " + reference

	#if mit isn't no, pick a mitochondrial reference to use
	if not (mit.rstrip("\n") == "no"):

		mit_references = []

		with open(mit) as file:

			for line in file:

				reference_name = line.rstrip("\n").split("\t")[0]

				reference_path = line.rstrip("\n").split("\t")[1]

				mit_references.append([reference_name, reference_path])

				print "Path to " +  reference_name +" reference is " + reference_path


	if (mit.rstrip("\n") == "no"):

		mit_references = "no"

	#define default AdapterRemoval

	adapter_removal = "/home/kdaly/programs/adapterremoval-2.3.1/build/AdapterRemoval --threads 2 --collapse --minadapteroverlap 1 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minlength 30 --gzip --trimns --trimqualities "

	#Prepare output directory

	out_dir = output_dir + date_of_hiseq +  "/"

	call("mkdir " + out_dir, shell=True)

	#allow meyer option to be used
	meyer_input = meyer.rstrip("\n").lower()

	alignment_option = "bwa aln -l 1024 -t " + threads + " "

	if (meyer_input == "meyer" or meyer_input == "yes"):

		print "Meyer option selected."

		alignment_option = "bwa aln -l 1024 -n 0.01 -o 2 -t " + threads + " "

	#variable for RG file
	RG_file = RG_file.rstrip("\n")

	fastq_list = []

	#create a list of all fastq files
	for file in files:

		print file

		current_file = "_".join(file.split("_")[:-1])

		fastq_list.append(current_file.rstrip("\n"))

	return files, reference, mit_references, out_dir, adapter_removal, alignment_option, fastq_list


def prepare_trim_fastq(current_sample, adapter_removal, out_dir):

	print "Trimming: current sample is: " + current_sample

	zipped_fastq1 = current_sample + "_r1.fastq.gz"

	zipped_fastq2 = current_sample + "_r2.fastq.gz"

	#Get number of lines (and from that reads - divide by four) from raw fastq
	trimmed_fastq1 = current_sample + "_r1_trimmed" + ".fastq.gz"

	trimmed_fastq2 = current_sample + "_r2_trimmed" + ".fastq.gz"

	#cut raw fastq files
	call("echo " + adapter_removal  + "--file1 " + zipped_fastq1 + " --file2 " + zipped_fastq2 + " --basename " + current_sample + "_trimmed \">\" " + current_sample + "_trimmed.log \"&&\" mv " + current_sample + "_trimmed.pair1.truncated.gz " + current_sample + "_r1_trimmed.fastq.gz \"&&\" mv " + current_sample + "_trimmed.pair2.truncated.gz " + current_sample + "_r2_trimmed.fastq.gz \"&&\" mv " + current_sample + "_trimmed.singleton.truncated.gz " + current_sample + "_mate-discard_trimmed.fastq.gz \"&&\" mv " + current_sample + "_trimmed.collapsed.gz " + current_sample + "_trimmed.fastq.gz >> trim.sh" ,shell=True)


def align_process_mit(fastq, RG_file, alignment_option, reference, trim):

	reference_sequence = reference[0]

	reference_path = reference[1]

	sample =  fastq.split(".")[0].replace("_r1","")

	sample_and_ref = sample + "_" + reference_sequence

	trimmed_fastq = sample + "_trimmed.fastq.gz"

	trimmed_r1 = sample + "_r1_trimmed.fastq.gz"

	trimmed_r2 = sample + "_r2_trimmed.fastq.gz"

	print alignment_option + reference_path + " " + trimmed_fastq + " > " + sample_and_ref + "_mit.sai 2>>"+ sample_and_ref + "_mit_alignment.log"
	call(alignment_option + reference_path + " " + trimmed_fastq + " > " + sample_and_ref + "_mit.sai 2>>"+ sample_and_ref + "_mit_alignment.log",shell=True)

	print alignment_option + reference_path + " " + trimmed_r1 + " > " + sample_and_ref + "_r1_mit.sai 2>>"+ sample_and_ref + "_r1_mit_alignment.log"
	call(alignment_option + reference_path + " " + trimmed_r1 + " > " + sample_and_ref + "_r1_mit.sai 2>>"+ sample_and_ref + "_r1_mit_alignment.log",shell=True)

	print alignment_option + reference_path + " " + trimmed_r2 + " > " + sample_and_ref + "_r2_mit.sai 2>>"+ sample_and_ref + "_r2_mit_alignment.log"
	call(alignment_option + reference_path + " " + trimmed_r2 + " > " + sample_and_ref + "_r2_mit.sai 2>>"+ sample_and_ref + "_r2_mit_alignment.log",shell=True)


	with open(RG_file) as file:

		print "Looking for RG. Current sample is " + sample

		for line in file:

			split_line = line.split("\t")

	                if (sample.replace("_r1","") == split_line[0].split(".")[0].replace("_r1","")):

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
        print sample + " aligning to " + reference_path

        print RG
        print "bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference_path + " " + sample_and_ref + "_mit.sai " + trimmed_fastq + " | samtools view -Sb -F 4 - > " + sample_and_ref + "_mit_F4.bam + 2> " + sample_and_ref + "_mit_alignment.log"
        call("bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference_path + " " + sample_and_ref + "_mit.sai " + trimmed_fastq + " | samtools view -Sb -F 4 - > " + sample_and_ref +"_mit_F4.bam 2> " + sample_and_ref + "_mit_alignment.log", shell=True)

	print "bwa sampe -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference_path + " " + sample_and_ref + "_r1_mit.sai "  + " " + sample_and_ref + "_r2_mit.sai " + trimmed_r1 + " " + trimmed_r2 + " | samtools view -Sb -F 4 - > " + sample_and_ref + "_pe_mit_F4.bam + 2> " + sample_and_ref + "_pe_mit_alignment.log"
	call( "bwa sampe -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference_path + " " + sample_and_ref + "_r1_mit.sai "  + " " + sample_and_ref + "_r2_mit.sai " + trimmed_r1 + " " + trimmed_r2 +  " | samtools view -Sb -F 4 - > " + sample_and_ref + "_pe_mit_F4.bam + 2> " + sample_and_ref + "_pe_mit_alignment.log",shell=True)

	call("samtools flagstat " + sample_and_ref +"_mit_F4.bam > " +sample_and_ref + "_mit_F4.flagstat 2>> " + sample_and_ref + "_mit_alignment.log",shell=True)
	call("samtools flagstat " + sample_and_ref +"_pe_mit_F4.bam > " +sample_and_ref + "_pe_mit_F4.flagstat 2>> " + sample_and_ref + "_pe_mit_alignment.log",shell=True)

	call ("rm "+ sample_and_ref + "*_mit.sai ",shell=True)

	print "samtools sort -@ 24 "  + sample_and_ref +"_mit_F4.bam -O BAM -o " + sample_and_ref + "_mit_F4_sort.bam 2>>" + sample_and_ref + "_mit_alignment.log"
	call("samtools sort -@ 24 "  + sample_and_ref +"_mit_F4.bam -O BAM -o " + sample_and_ref + "_mit_F4_sort.bam 2>> " + sample_and_ref + "_mit_alignment.log",shell=True)

	print "samtools sort -n -@ 24 "  + sample_and_ref + "_pe_mit_F4.bam -O BAM -o " + sample_and_ref + "_pe_mit_F4_sort.bam 2>> " + sample_and_ref + "_pe_mit_alignment.log"
	call("samtools sort -n -@ 24 "  + sample_and_ref + "_pe_mit_F4.bam -O BAM -o " + sample_and_ref + "_pe_mit_F4_sort.bam 2>> " + sample_and_ref + "_pe_mit_alignment.log",shell=True)

	print "samtools fixmate -m -@ 24 "  + sample_and_ref + "_pe_mit_F4_sort.bam " + sample_and_ref + "_pe_mit_F4_sort_fixmate.bam"
	call("samtools fixmate -m -@ 24 "  + sample_and_ref + "_pe_mit_F4_sort.bam " + sample_and_ref + "_pe_mit_F4_sort_fixmate.bam",shell=True)

	print "samtools sort -@ 24 " + sample_and_ref + "_pe_mit_F4_sort_fixmate.bam " + sample_and_ref + "_pe_mit_F4_sort_fixmate_resort.bam"
	call("samtools sort -@ 24 " + sample_and_ref + "_pe_mit_F4_sort_fixmate.bam " + sample_and_ref + "_pe_mit_F4_sort_fixmate_resort.bam",shell=True)

	print "samtools markdups -r "  + sample_and_ref +"_mit_F4_sort.bam " + sample_and_ref + "_mit_F4_markdup.bam 2>>" + sample_and_ref + "_mit_alignment.log"
	call("samtools markdups -r "  + sample_and_ref +"_mit_F4_sort.bam " + sample_and_ref + "_mit_F4_markdup.bam 2>> " + sample_and_ref + "_mit_alignment.log",shell=True)

	print "samtools markdups -r "  + sample_and_ref +"_pe_mit_F4_sort_fixmate_resort.bam " + sample_and_ref + "_pe_mit_F4_markdup.bam 2>>" + sample_and_ref + "_pe_mit_alignment.log"
	call("samtools markdups -r "  + sample_and_ref +"_pe_mit_F4_sort_fixmate_resort.bam " + sample_and_ref + "_pe_mit_F4_markdup.bam 2>>" + sample_and_ref + "_pe_mit_alignment.log",shell=True)

	call("rm " + sample_and_ref + "*_mit_F4_sort.bam",shell=True)

	call("samtools flagstat " + sample_and_ref + "_mit_F4_markdup.bam > " + sample_and_ref + "_mit_F4_markdup.flagstat",shell=True)
	call("samtools flagstat " + sample_and_ref + "_pe_mit_F4_markdup.bam > " + sample_and_ref + "_pe_mit_F4_markdup.flagstat",shell=True)


def merge_and_process_mit(RG_file, reference, trim):

	reference_sequence = reference[0]

	reference_path = reference[1]

	#merge each lane then each sample
	#account for the fact that we are aligning to different mitochondrial refereneces

	merged_mit_bam_list = merge_lanes_and_sample(RG_file, trim,reference_sequence,"yes",reference_sequence)

	for bam in merged_mit_bam_list:

		bam_root = bam.split(".")[0]

		print "samtools flagstat " + bam + "  > " + bam_root + ".flagstat"
		call("samtools flagstat " + bam + "  > " + bam_root + ".flagstat",shell=True)

		print "samtools markdups -r " + bam_root + ".bam " + bam_root + "_markdup.bam "
		call("samtools markdups -r " + bam_root + ".bam " + bam_root + "_markdup.bam ",shell=True)

		print "samtools flagstat " + bam_root + "_markdup.bam > " + bam_root + "_markdup.flagstat"
		call("samtools flagstat " + bam_root + "_markdup.bam > " + bam_root + "_markdup.flagstat",shell=True)

		#filter for just q30
		for QC in ["30"]:

			print "samtools view -b -q" + QC + " " +  bam_root + "_markdup.bam > " + bam_root + "_markdup_q" + QC + ".bam"
			call("samtools view -b -q" + QC + " " +  bam_root + "_markdup.bam > " + bam_root + "_markdup_q" + QC + ".bam",shell=True)

			print "samtools index -@ 24 " + bam_root + "_markdup_q" + QC + ".bam"
			call("samtools index -@ 24 " + bam_root + "_markdup_q" + QC + ".bam",shell=True)

			cmd="java -Xmx10g -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -nt 24 -R " + reference_path + " -o DoC_" + bam_root + "_q" + QC + " -I " + bam_root + "_markdup_q" + QC + ".bam --omitDepthOutputAtEachBase --omitIntervalStatistics"
			print cmd
			call(cmd, shell=True)

			print "samtools idxstats " +  bam_root + "_markdup_q" + QC + ".bam >" + bam_root + "_markdup_q" + QC + ".idx"
			call("samtools idxstats " +  bam_root + "_markdup_q" + QC + ".bam >" + bam_root + "_markdup_q" + QC + ".idx",shell=True)

			for minD in ["1", "2", "3"]:

				print "angsd -doFasta 2 -i " + bam_root + "_markdup_q" + QC + ".bam  -doCounts 1 -out " + bam_root + "_angsd-consensus-min" + minD + "_q" + QC + " -setMinDepth " + minD + " -minQ 20 minMapQ " + QC + " -trim 4"
				call("angsd -doFasta 2 -i " + bam_root + "_markdup_q" + QC + ".bam  -doCounts 1 -out " + bam_root + "_angsd-consensus-min" + minD + "_q" + QC + " -setMinDepth " + minD + " -minQ 20 -minMapQ " + QC + " -trim 4 ",shell=True)

				call("gunzip " + bam_root + "_angsd-consensus-min" + minD + "_q" + QC + ".fa.gz; python ~/programs/scripts_for_goat_project/decircularize.py "  + bam_root + "_angsd-consensus-min" + minD + "_q" + QC + ".fa > " + bam_root + "_angsd-consensus-min" + minD + "_q" + QC + "_decirc.fa",shell=True)

		call("mkdir " + bam_root + "_angsd-consensus ; mv *angsd-conse*fa *angsd-conse*arg *angsd*.fa* -t " + bam_root + "_angsd-consensus",shell=True)


def align_bam(sample, RG_file, alignment_option, reference, trim, species):

	print "bam alignment"

	if (trim == "yes"):

		sample = sample + "_dummy"

	#used to do mate-discard, but I think it's best if it's removed and not included
	#for sample_split in ["_".join(sample.split("_")[:-1]) + "_mate-discard", "_".join(sample.split("_")[:-1])]:

	for sample_split in ["_".join(sample.split("_")[:-1])]:

		trimmed_fastq = sample_split + "_trimmed.fastq.gz"

		sample_ref = sample_split + "_" + species

		print(alignment_option + reference + " " + trimmed_fastq + " > " + sample_ref + ".sai")
		call(alignment_option + reference + " " + trimmed_fastq + " > " + sample_ref + ".sai 2>" + sample_ref + "_alignment.log",shell=True)

		with open(RG_file) as file:

			for line in file:

				split_line = line.split("\t")

				if (sample_split.replace("_mate-discard","")  == split_line[0].split(".")[0].replace("_r1","")):

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

		#we want to run the samse step in parallel, also flagstat
		#big question if the -@ 2 is actually efficient with parallel
		#need to use the print command due to needing the ' symbol
		f=open("samse.sh","a+")
		f.write("bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference + " " + sample_ref + ".sai " + trimmed_fastq + " | samtools view -@ 2 -Sb - > " + sample_ref + ".bam 2> "+ trimmed_fastq + "_" + species + "_alignment.log; samtools flagstat -@ 2 " + sample_ref + ".bam > " + sample_ref + ".flagstat\n")
		f.close()

		call("echo bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference + " " + sample_ref + ".sai " + trimmed_fastq + " \"|\" samtools view -@ 2 -Sb - \">\" " + sample_ref + ".bam \"2>\" " + sample_ref + "_alignment.log\";\" samtools flagstat -@ 2 " + sample_ref + ".bam \">\" " + sample_ref + ".flagstat",shell=True)
		#call("echo bwa samse -r \\'" + RG.rstrip("\n").replace("+", "\\t") + "\\' " + reference + " " + sample + ".sai " + trimmed_fastq + " \"|\" samtools view -@ 2 -Sb - \">\" " + sample + ".bam \"2>\" "+ trimmed_fastq + "_alignment.log\";\" samtools flagstat -@ 2 " + sample + ".bam \">\" " + sample + ".flagstat >> samse.sh",shell=True)


def align_bam_pe(sample, RG_file, alignment_option, reference, trim, species):

        print "bam alignment"

	sample = sample.replace("_r1","")

	for trimmed_fastq in [ sample + "_r1_trimmed.fastq.gz", sample + "_r2_trimmed.fastq.gz"]:

		sample_ref =  "_".join(trimmed_fastq.split("_")[:-1]) + "_" + species

		print(alignment_option + reference + " " + trimmed_fastq + " > " + sample_ref + ".sai")
		call(alignment_option + reference + " " + trimmed_fastq + " > " + sample_ref + ".sai 2>" + sample_ref + "_alignment.log",shell=True)

	with open(RG_file) as file:

		for line in file:

			split_line = line.split("\t")

			print sample

			print split_line[0].split(".")[0]

			if (sample  == split_line[0].split(".")[0].replace("_r1","")):

				RG = split_line[1].rstrip("\n")

				if not RG:

					print "No RGs were detected for this sample - please check sample names in fastq files and in RG file agree"

					break

				else:

					print "Reads groups being used are:"

					print RG

		file.seek(0)

		f=open("sampe.sh","a+")

		f.write("bwa sampe -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference + " " + sample + "_r1_" + species + ".sai " +  sample + "_r2_" + species + ".sai " +  sample + "_r1_trimmed.fastq.gz " +  sample + "_r2_trimmed.fastq.gz  | samtools view -@ 2 -Sb - > " + sample + "_pe_" + species + ".bam 2> "+ trimmed_fastq + "_" + species + "_pe_alignment.log; samtools flagstat -@ 2 " + sample + "_pe_" + species + ".bam > " + sample + "_pe_" + species + ".flagstat\n")



def process_bam(sample_name,species):

	#the input list will be different depending on whether the fastq files have been trimmed prior to the script being run
	#ie if the script was interupted and had to be restarted
	#including a fix for that

	if sample_name.endswith("_trimmed"):

		sample_name = "_".join(sample_name.split("_")[:-1])

	sample_name = sample_name.replace("_r1","")

	print "Processing step for: " + sample_name

	#for sample_name in ["_".join(sample_name.split("_")[0]) + "_" + species, "_".join(sample_name.split("_")[0]) + "_pe_" + species]:

	for sample_name in [sample_name + "_" + species, sample_name + "_pe_" +  species]:

		print "With species: " + sample_name

		sample_name.rstrip("_")

		#if the bam is gzipped, gunzip
		if os.path.isfile(sample_name + ".bam.gz"):

			call("gunzip " + sample_name + ".bam.gz",shell=True)

		#flagstat the bam
		print "samtools flagstat -@ 24 " + sample_name + ".bam > " + sample_name + ".flagstat"
		call("samtools flagstat -@ 24 " + sample_name + ".bam > " + sample_name + ".flagstat",shell=True)

		#sort this bam
		if not "_pe_" in sample_name:
			print "samtools sort -@ 24 " + sample_name + ".bam -O BAM -o " + sample_name + "_sort.bam"
			call("samtools sort -@ 24 " + sample_name + ".bam -O BAM -o " + sample_name + "_sort.bam",shell=True)

		else:

			print "samtools sort -@ 24 -n " + sample_name + ".bam -O BAM -o " + sample_name + "_sort.bam"
			call("samtools sort -@ 24 -n " + sample_name + ".bam -O BAM -o " + sample_name + "_sort.bam",shell=T)

		print "samtools view -@ 24 -q20 -b "  + sample_name + "_sort.bam > "  + sample_name + "_sort_q20.bam"
		call("samtools view -@ 24  -q20 -b "  + sample_name + "_sort.bam > "  + sample_name + "_sort_q20.bam" ,shell=True)

		print "rm " + sample_name + "_sort.ba*"
		call("rm " + sample_name + "_sort.ba*" ,shell=True)

		print "samtools flagstat -@ 24 " + sample_name + "_sort_q20.bam > " + sample_name + "_sort_q20.flagstat"
		call("samtools flagstat -@ 24 " + sample_name + "_sort_q20.bam > " + sample_name + "_sort_q20.flagstat",shell=True)

		if "_pe_" in sample_name:

			print "samtools fixmate -m -@ 24 " + sample_name + "_sort_q20.bam " +  sample_name + "_sort_q20_fix.bam"
			call("samtools fixmate -m -@ 24 " + sample_name + "_sort_q20.bam " +  sample_name + "_sort_q20_fix.bam",shell=T)

			print "samtools sort -@ 24 " +  sample_name + "_sort_q20_fix.bam -O BAM -o " +  sample_name + "_sort_q20_fix_resort.bam"
			call("samtools sort -@ 24 " +  sample_name + "_sort_q20_fix.bam -O BAM -o " +  sample_name + "_sort_q20_fix_resort.bam",shell=T)

		#gzip the original bam
		call("gzip " + sample_name + ".bam",shell=True)

def merge_lanes_and_sample(RG_file, trim, species,mit="no", mit_reference="no"):

	#get sample list from the RG file

	sample_list = []

	with open(RG_file) as r:

        	for line in r:

            		sample = line.split("\t")[2].rstrip("\n")

			if [sample] not in sample_list:

				sample_list.append([sample])

	print "List of samples for merging"
	print sample_list

	#cycle through the RG file and associate each lane with the correct sample
	for sample in sample_list:

		lane_list = []

		with open(RG_file) as r:

			for line in r:

				if sample[0] == line.split("\t")[2].rstrip("\n"):

					lane = line.split("\t")[2]

					if lane not in lane_list:

						lane_list.append(lane)

		sample.append(lane_list)

	#create a list of final merged,markdup bams that will be returned
	merged_bam_list = []

	#Merge each bam for each sample
	for sample in sample_list:

		merged_sample_list = []

		for lane in sample[1]:

			sample_files = []

			merge_cmd = "java -Xmx100g -jar /home/kdaly/programs/picard/picard.jar MergeSamFiles VALIDATION_STRINGENCY=SILENT "

			with open(RG_file) as r:

				for line in r:

					if (sample[0] == line.split("\t")[2].rstrip("\n")):

						#at this stage I have an issue with naming the sample
						#need to straighten out the name of the sample depending on if I have already trimmed prior to running the script
						if (mit == "yes"):

							if (trim=="no"):

								sample_files.append(line.split("\t")[0].split(".")[0].replace("_r1","")  + mit_reference +  "_mit_F4_markdup.bam")

								sample_files.append(line.split("\t")[0].split(".")[0].replace("_r1","") +  "_pe_mit_F4_markdup.bam")

							else:

								sample_files.append(line.split("\t")[0].split(".")[0].replace("_r1","") +  "_mit_F4_markdup.bam")

								sample_files.append(line.split("\t")[0].split(".")[0].replace("_r1","") +  "_pe_mit_F4_markdup.bam")


						#may have to deal with trimmed files here
						else:

							sample_files.append(line.split("\t")[0].split(".")[0] + "_" + species + "_q20_markdup.bam")


			#create a "sample name" variable to apply to final bams
			for bam in sample_files:

				print "Checking for bam: " + bam

				if os.path.isfile(bam):

					merged_sample_list.append(bam)
		#now, merge each lane for a given sample

		merge_cmd = "java -Xmx100g -jar /home/kdaly/programs/picard/picard.jar  MergeSamFiles VALIDATION_STRINGENCY=SILENT "

		for bam in merged_sample_list:

			merge_cmd = merge_cmd + "INPUT=" + bam + " "

		sample_name =  sample[0] + "_" + species

		if (mit == "yes"):

			sample_name = sample_name + "_" + mit_reference + "_mit"

		merge_cmd = merge_cmd + "OUTPUT=" + sample_name + "_q20_merged.bam 2>" + sample_name + "_q20_merged.log"

		print "Current merge command is: "
		print merge_cmd

		call(merge_cmd,shell=True)

		#flagstat the merged bam
		call("samtools flagstat -@ 20 " + sample_name + "_q20_merged.bam > " + sample_name+ "_q20_merged.flagstat",shell=True)

		call("samtools view -b -q 30 -@ 20 " + sample_name + "_q20_merged.bam > " + sample_name + "_merged_q30.bam",shell=True)

		#call("samtools markdups -r " + sample_name + "_merged_q30.bam " + sample_name + "_merged_q30_markdup.bam 2> " + sample_name + "_merged_q30_markdup.log",shell=True)

		call("samtools flagstat -@ 20 " + sample_name + "_q20_merged.bam > " + sample_name + "_q20_merged.flagstat",shell=True)

		call("samtools flagstat -@ 20 " + sample_name + "_merged_q30.bam > " + sample_name + "_merged_q30.flagstat",shell=True)

		#call("samtools flagstat -@ 20 " + sample_name + "_merged_q30_markdup.bam > " + sample_name + "_merged_q30_markdup.flagstat",shell=True)

		merged_bam_list.append(sample_name + "_merged_q30.bam")

	return merged_bam_list

def indel_realignment(markdup_bam, reference_genome):

	print "starting realignment for sample "+markdup_bam

	call("samtools index -@ 24 " + markdup_bam,shell=True)

	cmd = "java -Xmx20g -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 24 -R " + reference_genome + " -I " + markdup_bam + " -o forIndelRealigner_" + markdup_bam.split(".")[0] + ".intervals 2> " + markdup_bam.split(".")[0] + "_intervals.log"

	print cmd
	call(cmd, shell=True)

	cmd = "java -Xmx100g -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R " + reference_genome + " -I " + markdup_bam + " -targetIntervals forIndelRealigner_" + markdup_bam.split(".")[0] + ".intervals -o " +  markdup_bam.split(".")[0] + "_realigned.bam 2> " + markdup_bam.split(".")[0] + "_realignment.log"

	print cmd

	call(cmd,shell=True)

	return markdup_bam.split(".")[0] + "_realigned.bam"


def softclip_bam(bam,reference_genome, out_dir, to_clip = "4"):

	call("samtools view -h " + bam + " | python ~/programs/scripts_for_goat_project/softclip_mod.py - " + to_clip + " | samtools view -Sb - > " + bam.split(".")[0] + "_softclipped.bam",shell=True)

	call("echo samtools view -h " + bam + " \"|\" python ~/programs/scripts_for_goat_project/softclip_mod.py - " + to_clip + " \"|\" samtools view -Sb - \">\" " + bam.split(".")[0] + "_softclipped.bam",shell=True)

	return (bam.split(".")[0] + "_softclipped.bam")


def process_realigned_bams(realigned_bam, reference_genome, clip, output_dir, species):

	print "Realigned bam files is: "
	print realigned_bam

	if (clip == "yes"):

		#rescale at this point
		softclip_bam(realigned_bam,reference_genome,output_dir)

		realigned_bam = realigned_bam.split(".")[0] + "_softclipped.bam"


	#call("samtools view -b -F4 " + realigned_bam.split(".")[0] + ".bam > " + realigned_bam.split(".")[0] + "_F4.bam && samtools view -q30 -b " + realigned_bam.split(".")[0] + "_F4.bam > " + realigned_bam.split(".")[0] + "_F4_q30.bam",shell=True)

	#gzip the markdup bam
	#call("gzip " + realigned_bam ,shell=True)
	call("samtools flagstat -@ 24 "  + realigned_bam + " > " + realigned_bam.split(".")[0] + ".flagstat",shell=True)
	#call("samtools flagstat " + realigned_bam.split(".")[0] + "_F4_q30.bam > " + realigned_bam.split(".")[0] + "_F4_q30.flagstat",shell=True)

	call("samtools index -@ 24 "  + realigned_bam,shell=True)
	#call("samtools index -@ 24 " + realigned_bam.split(".")[0] + "_F4_q30.bam",shell=True)

	#call("samtools idxstats " + realigned_bam.split(".")[0] + "_F4_q30.bam > " + realigned_bam.split(".")[0].split("_")[0] + ".idx",shell=True)

	cmd="java -Xmx10g -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -nt 24 --omitIntervalStatistics -R " + reference_genome + " -o DoC_" + realigned_bam.split(".")[0] + " -I " + realigned_bam + "  --omitDepthOutputAtEachBase 2>  DoC_" + realigned_bam.split(".")[0] + ".log"
	#cmd="java -Xmx10g -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -nt 24 --omitIntervalStatistics -R " + reference_genome + " -o DoC_" + realigned_bam.split(".")[0] + " -I " + realigned_bam.split(".")[0] + "_F4_q30.bam --omitDepthOutputAtEachBase"

	call("sh /home/kdaly/programs/scripts_for_goat_project/run_DoC_autosomes.sh " + reference_genome + " " + realigned_bam,shell=True)


def clean_up(out_dir):

	#clean up files
	call("gunzip *flagstat.gz",shell=True)

	call("mkdir flagstat_files; mv *flagstat -t flagstat_files",shell=True)

	call("mkdir DoC; mv DoC_* DoC",shell=True)

	call("mkdir log_files; mv *log log_files", shell=True)

	call("mkdir trimmed_fastq_files_and_logs",shell=True)

	call("mv *trimmed*  -t trimmed_fastq_files_and_logs/",shell=True)

	call("mkdir idx_files; mv *idx* -t idx_files; mkdir auxillary_files; mv *txt *interval* RG* *md5sum* -t auxillary_files",shell=True)

	call("gzip *bam",shell=True)

	call("mkdir final_bams ; mv *rescaled* *softcli* *realigned.b* -t final_bams/ ; mv final_mit_bams final_bams -t " + out_dir + "; mkdir intermediate_bams; mv *bam* *bai -t intermediate_bams",shell=True)

	call("gzip trimmed_fastq_files_and_logs/*",shell=True)

	call("mkdir mapDamage; mv results_* -t mapDamage/",shell=True)

	call("mv mit_idx_files mit_logs flagstat_files log_files angsd_consensus_sequences trimmed_fastq_files_and_logs idx_files auxillary_files intermediate_bams mapDamage DoC fastq_files -t " + out_dir,shell=True)

	call("gzip intermediate_bams/*bam",shell=True)

def clean_up_mit(mitochondrial_references_file,out_dir):

	#make output directories and dump files
	call("mkdir auxillary_files; mv " + mitochondrial_references_file + " auxillary_files" ,shell=True)

	call("gzip *.bam",shell=True)

	call("mkdir mit_DoC; mv *DoC* mit_DoC",shell=True)

	call("mkdir angsd_consensus; mv *angsd* -t angsd_consensus",shell=True)

	call("mkdir mit_logs; mv *mit*.log -t mit_logs; mv *flagstat* -t flagstat_files; mkdir mit_idx_files; mv *mit*idx -t mit_idx_files", shell=True)

	call("bgzip *mit*bam.gz; mkdir final_mit_bams; mv *mit*q30.bam* -t final_mit_bams; mkdir intermediate_mit_bam_files; mv *_mit*.bam.gz -t intermediate_mit_bam_files ",shell=True)

	call("mv trimmed_fastq_files_and_logs mit_DoC mit_logs flagstat_files mit_idx_files final_mit_bams intermediate_mit_bam_files angsd_consensus -t " + out_dir,shell=True)

try:
	date_of_hiseq  = sys.argv[1]
	meyer = sys.argv[2]
	threads=sys.argv[3]
	species = sys.argv[4]
	mit = sys.argv[5]
	skip_mit_align = sys.argv[6]
	trim = sys.argv[7]
	align = sys.argv[8]
	process = sys.argv[9]
	merge = sys.argv[10]
	clip = sys.argv[11]
	RG_file  = sys.argv[12]
	output_dir = sys.argv[13]

except IndexError:
	print "Incorrect number of variables have been provided"
	print "Input variables are date_of_hiseq, meyer, number of threads, reference genome, mit, skip_mit_align, trim, align, process, merge, clip, RG_file, and the directory to put output directories/files"
	print "Exiting program..."
	sys.exit()

if not output_dir[-1] == "/":

	output_dir = output_dir + "/"

main(date_of_hiseq, meyer, threads, species, mit, skip_mit_align, trim, align, process, merge, clip, RG_file, output_dir)
