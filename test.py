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

