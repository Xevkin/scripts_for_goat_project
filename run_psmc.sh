#!/bin/sh
#input must be in the form of an input file, then path to the folder with the bams of interest, the reference genome, and the gen time



INPUT_FILE=$1
INPUT_FOLD=$2
REF=$3

#must use the correct/most recent version of samtools - bcftools will complain if otherwise
SAMTOOLS="/home/kdaly/bin/samtools-1.2/samtools"

PSMC="/home/kdaly/bin/psmc/"

#using the time interval settings, generation time and mutation rate as suggested by Bosse et al, 2015 
PSMC_SETTINGS="-N25 -t15 -r5 -p '4+50*1+4+6'"
GEN_TIME=$4
MUTATION_RATE="1e-8"

#input file should have one line per sample, space seperated
#first col is sample name without bam handle
#second col is minimum depth (one third of mean coverage)
#third col is max depth (double mean coverage)
#eg KD100_rescaled_rmdup 8 40 0.1

while read sample min_depth max_depth FNR
do
	#run samtools pileup, then bcftools to create correct .vcf for psmc
	$SAMTOOLS mpileup -C50 -uf $REF $INPUT_FOLD$sample".bam" | bcftools call -c -O v -o $sample".vcf" -
	
	#run vcfutils.pl, a script that converts vcf file to  .fq format
	#this is whole diploid concensus sequence for the individual
	#we filter variants passed on coverage - minimum of one third mean depth, max of twice mean depth

	vcfutils.pl vcf2fq -d $min_depth -D $max_depth $sample".vcf" | gzip > $sample".fq.gz"
	
	#Now we create the bin file  that is the input to psmc
	#this is a .psmcfa file, which actually looks like a fasta file. Each character represents a bin of size 100bp
        #N indicates too many N's in the bin; T indicates a homozygous bin; K is a heterozygous bin

	SAMPLE_NAME=`echo $sample | cut -f1 -d'_'`

	$PSMC"utils/fq2psmcfa" -q20 $sample".fq.gz" > $SAMPLE_NAME".psmcfa"
	
	#We now run psmc based on the settings above
	#note this does not bootstrap - do seperately
	#using a generic FNR for now
	echo $PSMC"psmc" $PSMC_SETTINGS -o $SAMPLE_NAME".psmc" $SAMPLE_NAME".psmcfa" > tmp.sh; sh tmp.sh
	#echo $PSMC"psmc" -N25 -t15 -r5 -p \"4+50*1+4+6\" -o $SAMPLE_NAME".psmc" $SAMPLE_NAME".psmcfa" > tmp.sh; sh tmp.sh
	$PSMC"utils/psmc2history.pl" $SAMPLE_NAME".psmc" | $PSMC"utils/history2ms.pl" > ms-cmd.sh
	echo $PSMC"utils/psmc_plot.pl" -u $MUTATION_RATE -g $GEN_TIME -M '"'$SAMPLE_NAME"="$FNR'"' $SAMPLE_NAME"-"$FNR $SAMPLE_NAME".psmc" >tmp.sh
	sh tmp.sh
	mv ms-cmd.sh $SAMPLE_NAME"-ms-cmd.sh"
	rm tmp.sh
done<$INPUT_FILE


