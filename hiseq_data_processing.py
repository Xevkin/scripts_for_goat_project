#!/usr/bin/env python
'''
Script is run in a directory with bam files to be aligned to either goat, wild goat or sheep genome
Reads will be trimmed, aligned, read groups added, duplicates removed, then 

Also need to supply a date for the Hiseq run - will automatically make results file
python script <date_of_hiseq> <meyer> <species> <mit_file> <trim> <align> <process_individual_bams> <merge+process> <rescale> <read group file> <directory in which to place output dir/>

mit_file should be tab deliminated, col1 with name of reference, col2 with path
read group file needs to be in the follow format (\t are actual tab characters)
FASTQ_FILE.GZ\t@RG+ID:X+SM:X+PL:X+LB:X\tLANE\tSAMPLE_NAME
ID should be in the format <sample_name>-<macrogen_index_number>-<lane_number>-<hiseq_number>
LB should refer to PCR: <sample>-<lab index>-<macrogen-index>-<PCR_number>

fastq files should be in the format <sample>-<PCR number and letter, if any>_<underscores and numbers>

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

	"goat1" : "/kendrick/reference_genomes/goat_CHIR1_0/goat_CHIR1_0.fasta",

	"sheep" : "/kendrick/reference_genomes/sheep_oviAri3/oviAri3.fa",

	"bezoar" : "/kendrick/reference_genomes/bezoar_CapAeg_1_0/CapAeg_renamed.fa",

	"goat1_with_mit" : "/kendrick/reference_genomes/goat_CHIR1_0/goat_CHIR1_with_mit.fa"
}


def main(date_of_hiseq, meyer, species, mit,trim, align, process, merge, rescale, RG_file, output_dir):
	
	#run the set up function.#set up will create some output directories
	#and return variables that will be used in the rest of the script
	
	files, reference, mit_references, out_dir, cut_adapt, alignment_option, fastq_list = set_up(date_of_hiseq, meyer, species, mit, RG_file, output_dir, trim) 
	
	#make a folder for flagstat files	
	call("mkdir flagstat_files",shell=True)

	#sample is the file root
	#trim fastq files and produce fastqc files - if we want to
	#the masterlist will change each time so it needs be equated to the function
	
	if (trim == "yes" or trim == "trim"):
 
		for fastq in fastq_list:
		
			trim_fastq(fastq, cut_adapt, out_dir)
	
	#at this stage we have our fastq files with adaptors trimmed
	#We can now move on to the next step: alignment
	
	#going to align to CHIR1.0, as that what was used for AdaptMap
	#this step will align each fastq and produce raw bams with RGs
	#they should have the same stem of the initial bam


	if (mit_references != "no" ):

		#align to every mitochondrial reference provided
		for mitochondrial_reference in mit_references:

			print "Doing mit alignment to " + mitochondrial_reference[0]

			map (lambda fastq : align_process_mit(fastq, RG_file, alignment_option, mitochondrial_reference, trim), fastq_list)

			merge_and_process_mit(RG_file, mitochondrial_reference)

		clean_up_mit(mit,out_dir)

	if (align == "yes" or align == "align"):

		map (lambda fastq : align_bam(fastq, RG_file, alignment_option, reference, trim), fastq_list)

	#Remove the sai files that we don't care about
	call("rm *sai",shell=True)

	#sort and remove duplicates from each bam	
	if (process == "yes" or process == "process"):

		map(process_bam, fastq_list)
	
	
	#add an option here to kill the script if you do not want merging to occur
	if (merge == "no"):

		sys.exit("Script is terminated as no merging was desired")	
 
    	#now merge the lanes for each sample
	#NOTE: if all the options up to process are "no", then this is where the script will start
	#expects rmdup bams for each index in each lane, ungziped
	merged_bam_list = merge_lanes_and_sample(RG_file)
	
	realigned_bam_list = []

	#indel realignment
	for i in merged_bam_list:

		print "Indel realignment on sample " + i

		realigned_bam_list.append(indel_realignment(i,reference))

	#now run the script to process the realigned bams
	#this function will invoke a function that rescales the bams
	#duplicates will also be removed at this point
	for i in realigned_bam_list:

		process_realigned_bams(i,reference,rescale,output_dir)
	
	clean_up(out_dir)
	

def set_up(date_of_hiseq, meyer, species, mit, RG_file, output_dir, trim):
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
	
	print "Species selected is " + species

	print "Path to reference genome is " + reference
	
	#if mit isn't no, pick a mitochondrial reference to use
	if not (mit == "no"):

		mit_references = []

		with open(mit) as file:

			for line in file:

				reference_name = line.rstrip("\n").split("\t")[0]

				reference_path = line.rstrip("\n").split("\t")[1]

				mit_references.append([reference_name, reference_path])

				print "Path to " +  reference_name +" reference is " + reference_path
	
	else:

		mit_reference = "no"
	
	#define default cut_adapt

	cut_adapt = "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 "

	#Prepare output directory

	out_dir = output_dir + date_of_hiseq +  "/"

	call("mkdir " + out_dir, shell=True)

	#allow meyer option to be used
	meyer_input = meyer.rstrip("\n").lower()

	alignment_option = "bwa aln -t 12 -l 1024 "  

	if (meyer_input == "meyer"):
		
		print "Meyer option selected."
		
		alignment_option = "bwa aln -t 8 -l 1000 -n 0.01 -o 2 " 

	#variable for RG file
	RG_file = RG_file.rstrip("\n")

	fastq_list = []
	
	#create a list of all fastq files
	for file in files:

		current_file = file.split(".")[0] 

		fastq_list.append(current_file.rstrip("\n"))
	
	return files, reference, mit_references, out_dir, cut_adapt, alignment_option, fastq_list


def trim_fastq(current_sample, cut_adapt, out_dir):
	
	print "Trimming: current sample is: " + current_sample
	
	zipped_fastq = current_sample + ".fastq.gz"
	
	#Get number of lines (and from that reads - divide by four) from raw fastq
	trimmed_fastq = current_sample + "_trimmed" + ".fastq" 
	
	#cut raw fastq files
	call(cut_adapt + zipped_fastq + " > " + trimmed_fastq + " 2> " + trimmed_fastq + ".log", shell=True)
	

def align_process_mit(fastq, RG_file, alignment_option, reference, trim):

    reference_sequence = reference[0]

    reference_path = reference[1]

    sample =  fastq.split(".")[0]

    sample_and_ref = fastq.split(".")[0] + "_" + reference_sequence
    
    trimmed_fastq = sample + "_trimmed.fastq"

    if (trim != "yes"):

        trimmed_fastq = sample + ".fastq"

        sample = "_".join(sample.split("_")[:-1])

    print alignment_option + reference_path + " " + trimmed_fastq + " > " + sample_and_ref + "_mit.sai 2>>"+ sample_and_ref + "_mit_alignment.log"
    call(alignment_option + reference_path + " " + trimmed_fastq + " > " + sample_and_ref + "_mit.sai 2>>"+ sample_and_ref + "_mit_alignment.log",shell=True)

    with open(RG_file) as file:

        print "Looking for RG. Current sample is " + sample

        for line in file:
		
		split_line = line.split("\t")
                
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
        print sample + " aligning to " + reference_path

        print RG
        print "bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference_path + " " + sample_and_ref + "_mit.sai " + trimmed_fastq + " | samtools view -Sb -F 4 - > " + sample_and_ref + "_mit_F4.bam + 2> " + trimmed_fastq + "_mit_alignment.log"
        call("bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference_path + " " + sample_and_ref + "_mit.sai " + trimmed_fastq + " | samtools view -Sb -F 4 - > " + sample_and_ref +"_mit_F4.bam", shell=True)

	call("samtools flagstat " + sample_and_ref +"_mit_F4.bam > " +sample_and_ref + "_mit_F4.flagstat 2>> " + sample_and_ref + "_mit_alignment.log",shell=True)

	call ("rm "+ sample_and_ref + "_mit.sai ",shell=True)

	print "samtools sort "  + sample_and_ref +"_mit_F4.bam " + sample_and_ref + "_mit_F4_sort 2>>" + sample_and_ref + "_mit_alignment.log"
	call("samtools sort "  + sample_and_ref +"_mit_F4.bam " + sample_and_ref + "_mit_F4_sort 2>> " + sample_and_ref + "_mit_alignment.log",shell=True)

	print "samtools rmdup -s "  + sample_and_ref +"_mit_F4_sort.bam " + sample_and_ref + "_mit_F4_rmdup.bam 2>>" + sample_and_ref + "_mit_alignment.log"
	call("samtools rmdup -s "  + sample_and_ref +"_mit_F4_sort.bam " + sample_and_ref + "_mit_F4_rmdup.bam 2>> " + sample_and_ref + "_mit_alignment.log",shell=True)
	
	call("rm " + sample_and_ref + "_mit_F4_sort.bam",shell=True)
	
	call("samtools flagstat " + sample_and_ref + "_mit_F4_rmdup.bam > " + sample_and_ref + "_mit_F4_rmdup.flagstat",shell=True)	

def merge_and_process_mit(RG_file, reference):

	reference_sequence = reference[0]

	reference_path = reference[1]

	#merge each lane then each sample
	#account for the fact that we are aligning to different mitochondrial refereneces
	
	merged_mit_bam_list = merge_lanes_and_sample(RG_file,"yes",reference_sequence)

	for bam in merged_mit_bam_list:

		bam_root = bam.split(".")[0]

		call("samtools flagstat " + bam + "  > " + bam_root + ".flagstat",shell=True)

		call("samtools rmdup -s " + bam_root + ".bam " + bam_root + "_rmdup.bam ",shell=True)

		call("samtools flagstat " + bam_root + "_rmdup.bam > " + bam_root + "_rmdup.flagstat",shell=True)

		#filter for both q25 and q30
		for QC in ["25","30"]:

			call("samtools view -b -q" + QC + " " +  bam_root + "_rmdup.bam > " + bam_root + "_rmdup_q" + QC + ".bam",shell=True)

			call("samtools index " + bam_root + "_rmdup_q" + QC + ".bam",shell=True)

			cmd="java -Xmx20g -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T DepthOfCoverage -R " + reference_path + " -o DoC_" + bam_root + "_q" + QC + " -I " + bam_root + "_rmdup_q" + QC + ".bam --omitDepthOutputAtEachBase"
		
			call(cmd, shell=True)
			
			call("samtools idxstats " +  bam_root + "_rmdup_q" + QC + ".bam >" + bam_root + "_rmdup_q" + QC + ".idx",shell=True)
			
			for minD in ["2", "3"]:

				call("angsd -doFasta 2 -i " + bam_root + "_rmdup_q" + QC + ".bam  -doCounts 1 -out " + bam_root + "_angsd-consensus-min" + minD + "_q" + QC + " -setMinDepth " + minD + " -minQ " + QC,shell=True)

				call("gunzip " + bam_root + "_angsd-consensus-" + minD + "_q" + QC + ".fa.gz; decircularize.py "  + bam_root + "_angsd-consensus-min" + minD + "_q" + QC + ".fa > " + bam_root + "_angsd-consensus-min" + minD + "_decirc_q" + QC + ".fa",shell=True)

		call("mkdir " + bam_root + "; mv *angsd-conse* " + bam_root,shell=True)

		
def align_bam(sample, RG_file, alignment_option, reference, trim):

    trimmed_fastq = sample + "_trimmed.fastq"

    if (trim != "yes"):

    	trimmed_fastq = sample + ".fastq"

	sample = "_".join(sample.split("_")[:-1]) 

    print(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai")
    call(alignment_option + reference + " " + trimmed_fastq + " > " + sample + ".sai 2>" + trimmed_fastq + "_alignment.log",shell=True)
    
    with open(RG_file) as file:

	for line in file:
		
		split_line = line.split("\t")	

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
                    		
	print "bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam 2>" + trimmed_fastq + "_alignment.log"
        call("bwa samse -r \'" + RG.rstrip("\n").replace("+", "\\t") + "\' " + reference + " " + sample + ".sai " + trimmed_fastq + " | samtools view -Sb - > " + sample + ".bam 2> "+ trimmed_fastq + "_alignment.log",shell=True)
	
        #flagstat the bam
        call("samtools flagstat " + sample + ".bam > " + sample + ".flagstat",shell=True)

        #sort this bam
        print "samtools sort " + sample + ".bam " + sample + "_sort"
        call("samtools sort " + sample + ".bam " + sample + "_sort",shell=True)

        #gzip the original bam
        call("gzip " + sample + ".bam",shell=True)

	#remove duplicates from the sorted bam
        print "samtools rmdup -s " + sample + "_sort.bam " + sample + "_rmdup.bam"

        call("samtools rmdup -s " + sample + "_sort.bam " + sample + "_rmdup.bam 2> " + trimmed_fastq + "_alignment.log", shell=True)

        #remove the "sorted with duplicates" bam
        call("rm " + sample + "_sort.bam",shell=True)

        #flagstat the rmdup bam
        call("samtools flagstat " + sample + "_rmdup.bam > " + sample + "_rmdup.flagstat",shell=True)


def process_bam(sample_name):

	#the input list will be different depending on whether the fastq files have been trimmed prior to the script being run
	#ie if the script was interupted and had to be restarted
	#including a fix for that

	if sample_name.endswith("_trimmed"):	
	
		sample_name = "_".join(sample_name.split("_")[:-1])		
			
	print "Processing step for: " + sample_name
	
	#flagstat the bam
	call("samtools flagstat " + sample_name + ".bam > " + sample_name + ".flagstat",shell=True)

	#sort this bam
	print "samtools sort " + sample_name + ".bam " + sample_name + "_sort"
	call("samtools sort " + sample_name + ".bam " + sample_name + "_sort",shell=True)
	
	#gzip the original bam
	call("gzip " + sample_name + ".bam",shell=True)
	
	print "samtools rmdup -s " + sample_name + "_sort.bam " + sample_name + "_rmdup.bam 2>" + sample_name + "_alignment.log"

	#remove duplicates from the sorted bam
	call("samtools rmdup -s " + sample_name + "_sort.bam " + sample_name + "_rmdup.bam 2> "  + sample_name + "_alignment.log", shell=True)
	
	#remove the "sorted with duplicates" bam
	call("rm " + sample_name + "_sort.bam",shell=True)

	#flagstat the rmdup bam
	call("samtools flagstat " + sample_name + "_rmdup.bam > " + sample_name + "_rmdup.flagstat",shell=True)

	
def merge_lanes_and_sample(RG_file, mit="no", mit_reference="no"):

	#get sample list from the RG file
	
	sample_list = []
	
	with open(RG_file) as r:

        	for line in r:
			
            		sample = line.split("\t")[3].rstrip("\n")
			
			if [sample] not in sample_list:
                        
				sample_list.append([sample])

	print sample_list
	
	#cycle through the RG file and associate each lane with the correct sample
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
	
	#for each sample, go through each lane for that sample and merge each bam for that sample
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
						if (mit == "yes"):
					
							files_in_lane.append(line.split("\t")[0].split(".")[0] + "_trimmed_" + mit_reference +  "_mit_F4_rmdup.bam")

						else:
						
							files_in_lane.append(line.split("\t")[0].split(".")[0] + "_rmdup.bam")

			
			#create a "sample name" variable to apply to final bams
			for bam in files_in_lane:
					
				if os.path.isfile(bam):
					
					merge_cmd = merge_cmd + "INPUT=" + bam + " "
											
					if not "-" in bam:
				
						sample_name = bam.split("_")[0]
 
					else:
			
						sample_name = bam.split("-")[0]

				else:
		
					print bam  + " is not in the current directory"

			if (mit == "yes"):

				sample_name = sample_name + "_"  + mit_reference +  "_mit"

			sample_lane = sample_name + "_"	+ lane + "_merged"		 

			merge_cmd = merge_cmd + "OUTPUT=" + sample[0]+ "_" + sample_lane + ".bam 2>" + sample[0] + "_" + sample_lane + ".log"

			print merge_cmd
		
			#now merge each bam file associated with a given lane
			call(merge_cmd,shell=True)

			#flagstat the merged bam
			call("samtools flagstat "+ sample[0]+ "_" + sample_lane + ".bam >"+ sample[0]+ "_" + sample_lane + ".flagstat" ,shell=True)
	
			#remove duplicates from the merged lane bam
			cmd = "samtools rmdup -s " + sample[0]+ "_" + sample_lane + ".bam " + sample[0]+ "_" + sample_lane + "_rmdup.bam 2> " + sample[0]+ "_" + sample_lane + "_rmdup.log"  

			call(cmd,shell=True)

			#gzip the original merged bam
			call("gzip " + sample[0]+ "_" + sample_lane + ".bam",shell=True)
		
			#flagstat the rmdup_merged file
			call("samtools flagstat "+ sample[0]+ "_" + sample_lane + "_rmdup.bam > "+ sample[0]+ "_" + sample_lane + "_rmdup.flagstat" ,shell=True)

			merged_lane_list.append(sample[0]+ "_" + sample_lane + "_rmdup.bam")
		
		#each lane has been merged
		#now, merge each lane for a given sample

		merge_cmd = "java -Xmx20g -jar /research/picard-tools-1.119/MergeSamFiles.jar VALIDATION_STRINGENCY=SILENT "
	
		for bam in merged_lane_list:

			merge_cmd = merge_cmd + "INPUT=" + bam + " "
		
		sample_name =  sample[0]

		if (mit == "yes"):

			sample_name = sample_name + "_" + mit_reference + "_mit"
	
		merge_cmd = merge_cmd + "OUTPUT=" + sample_name + "_merged.bam 2>" + sample_name + "_merged.log"

		print merge_cmd

		call(merge_cmd,shell=True)

		#flagstat the merged bam
		call("samtools flagstat " + sample_name + "_merged.bam > " + sample_name+ "_merged.flagstat",shell=True)
	
		merged_bam_list.append(sample_name + "_merged.bam")

	#note: output is not rmdup. rmdup occurs after indel realignment
	return merged_bam_list

def indel_realignment(rmdup_bam, reference_genome):
	
	print "starting realignment for sample "+rmdup_bam
	
	call("samtools index " + rmdup_bam,shell=True)
	
	cmd = "java -Xmx20g -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T RealignerTargetCreator -R " + reference_genome + " -I " + rmdup_bam + " -o forIndelRealigner_" + rmdup_bam.split(".")[0] + ".intervals 2> " + rmdup_bam.split(".")[0] + "_intervals.log"

	call(cmd, shell=True)
	
	call("java -Xmx20g -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T IndelRealigner -R " + reference_genome + " -I " + rmdup_bam + " -targetIntervals forIndelRealigner_" + rmdup_bam.split(".")[0] + ".intervals -o " +  rmdup_bam.split(".")[0] + "_realigned.bam 2> " + rmdup_bam.split(".")[0] + "_realignment.log",shell=True)

	return rmdup_bam.split(".")[0] + "_realigned.bam"


def mapDamage_rescale(bam,reference_genome, out_dir):

	call("mapDamage -i " + bam + " -r " + reference_genome + " --rescale --verbose",shell=True)
	
	#gzip the input bam
	call("gzip " + bam,shell=True)

	print "mv results_" + bam.split(".")[0] + "/" +bam.split(".")[0] + ".rescaled.bam ./" + bam.split(".")[0] + "_rescaled.bam"
	call("mv results_" + bam.split(".")[0] + "/" +bam.split(".")[0] + ".rescaled.bam ./" + bam.split(".")[0] + "_rescaled.bam",shell=True)

	print bam.split(".")[0] + "_rescaled.bam"

	return (bam.split(".")[0] + "_rescaled.bam")


def process_realigned_bams(realigned_bam, reference_genome, rescale, output_dir):
	
	print "Realigned bam files is: "
	print realigned_bam

	if (rescale == "yes"):

		#rescale at this point
		mapDamage_rescale(realigned_bam,reference_genome,output_dir)
	
		realigned_bam = realigned_bam.split(".")[0] + "_rescaled.bam"


	call("samtools view -b -F4 " + realigned_bam.split(".")[0] + ".bam > " + realigned_bam.split(".")[0] + "_F4.bam && samtools view -q25 -b " + realigned_bam.split(".")[0] + "_F4.bam > " + realigned_bam.split(".")[0] + "_F4_q25.bam",shell=True)

	#gzip the rmdup bam
	call("gzip " + realigned_bam ,shell=True)
	
	call("samtools flagstat " + realigned_bam.split(".")[0] + "_F4_q25.bam > " + realigned_bam.split(".")[0] + "_F4_q25.flagstat",shell=True)

	call("samtools index " + realigned_bam.split(".")[0] + "_F4_q25.bam",shell=True)

	call("samtools idxstats " + realigned_bam.split(".")[0] + "_F4_q25.bam > " + realigned_bam.split(".")[0].split("_")[0] + ".idx",shell=True)

	cmd="java -Xmx20g -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T DepthOfCoverage -R " + reference_genome + " -o DoC_" + realigned_bam.split(".")[0] + " -I " + realigned_bam.split(".")[0] + "_F4_q25.bam --omitDepthOutputAtEachBase"

	call(cmd,shell=True)
	

def clean_up(out_dir):

	#clean up files
	call("gunzip *flagstat.gz",shell=True)
	
	call("mkdir flagstat_files; mv *flagstat flagstat_files",shell=True)
	
	call("mkdir DoC; mv DoC_* DoC",shell=True)
	
	call("mkdir log_files; mv *log log_files", shell=True)
	
	call("mkdir angsd_consensus_sequences; mv *angsd* angsd_consensus_sequences",shell=True)
	
	call("mkdir trimmed_fastq_files_and_logs",shell=True)
	
	call("mv *trimmed* trimmed_fastq_files_and_logs/",shell=True)

	call("mkdir idx_files; mv *idxmv --help *idx.gz idx_files; mkdir auxillary_files; mv *txt *interval* RG.tsv* *md5sum* auxillary_files",shell=True)	
	
	call("bgzip *bam",shell=True)

	call("mkdir final_bams; mv *F4* finals_bams/; mv final_mit_bams final_bams " + out_dir + "; mkdir intermediate_bams; mv *bam* *bai intermediate_bams",shell=True)
	
	call("gzip trimmed_fastq_files_and_logs/*",shell=True)
	
	call("mkdir mapDamage; mv results_* mapDamage/",shell=True)
	
	call("mkdir fastqc; mv *fastqc* fastqc", shell=True)

	call("mv mit_idx_files mit_logs flagstat_files log_files angsd_consensus_sequences trimmed_fastq_files_and_logs idx_files fastqc auxillary_files intermediate_bams mapDamage DoC fastq_files " + out_dir,shell=True)

def clean_up_mit(mitochondrial_references_file,out_dir):

	#make output directories and dump files
	call("mkdir auxillary_files; mv " + mitochondrial_references_file + " auxillary_files",shell=True)
	
	call("gzip *.bam",shell=True) 

	call("mkdir mit_DoC; mv *DoC* mit_DoC",shell=True)
	
	call("mkdir mit_logs; mv *mit*.log mit_logs; mv *flagstat* flagstat_files; mkdir mit_idx_files; mv *mit*idx mit_idx_files", shell=True)
	
	call("bgzip *mit*bam.gz; mkdir final_mit_bams; mv *mit*q25*bam* *mit_merged_rmdup* final_mit_bams; mkdir intermediate_mit_bam_files",shell=True)
	
	call("mkdir angsd_consensus; mv *angsd* angsd_consensus; mv *_mit* intermediate_mit_bam_files ; mv intermediate_mit_bam_files/final_mit_bams ./",shell=True)

	call("mv mit_DoC mit_logs flagstat_files mit_idx_files final_mit_bams intermediate_mit_bam_files angsd_consensus " + out_dir,shell=True)

try:
	date_of_hiseq  = sys.argv[1]
	meyer = sys.argv[2]
	species = sys.argv[3]
	mit = sys.argv[4]
	trim = sys.argv[5]
	align = sys.argv[6]
	process = sys.argv[7]
	merge = sys.argv[8]
	rescale = sys.argv[9]
	RG_file  = sys.argv[10]
	output_dir = sys.argv[11]

except IndexError:
	print "Incorrect number of variables have been provided"
	print "Input variables are date_of_hiseq, meyer, species, mit, trim, align, process, merge, rescale, RG_file, and the directory to put output directories/files"
	print "Exiting program..."
	sys.exit()

if not output_dir[-1] == "/":

	output_dir = output_dir + "/"

main(date_of_hiseq, meyer, species, mit, trim, align, process, merge, rescale, RG_file, output_dir)
