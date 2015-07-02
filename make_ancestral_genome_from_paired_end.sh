#!/bin/sh

#Script assumes alignment to the deer_no_U.fq

#Generates a fasta file that can be used as an ancestral genome


REF="/deer_wg/deer_HS_data/deer_ref/deer_ref_no_U.fa"

for i in $(ls *_1.fq.gz | cut -f1 -d'_');

	do gunzip $i"_1.fq.gz" $i"_2.fq.gz";

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 $i"_1.fq" > $i"_1_trimmed.fq";

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 $i"_2.fq" > $i"_2_trimmed.fq";

	gzip $i"_1.fq" $i"_2.fq";

	bwa mem -M -t 16 $REF $i"_1_trimmed.fq" $i"_2_trimmed.fq" | samtools view -Sb - > $i"_to_deer.bam";

	gzip $i"_1_trimmed.fq" $i"_2_trimmed.fq";

	samtools sort  $i"_to_deer.bam" $i"_sort_to_deer";

	samtools rmdup $i"_sort_to_deer.bam" $i"_rmdup_to_deer.bam";

	rm $i"_sort_to_deer.bam"; gzip $i"_to_deer.bam"; 
	
	angsd -i $i"_rmdup_to_deer.bam" -doFasta 1;

	gzip $i"_rmdup_to_deer.bam"; 

done
