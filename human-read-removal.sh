BAM_ROOT=`echo $1 | cut -f1 -d'.'`

SAMPLE=`echo $1 | cut -f1 -d'_'`

REF=/raid_md0/Reference_Genomes/For_Fastq_Screen/hg38.fa

ARS1=/raid_md0/Reference_Genomes/goat/ARS1.fa

samtools bam2fq  $1 | gzip -c - > ${BAM_ROOT}.fastq.gz

python2 ~/programs/scripts_for_goat_project/quick_alignment.py ${BAM_ROOT}.fastq.gz $REF

samtools view ${SAMPLE}_ARS1_merged_q30_hg38_F4.bam |cut -f1 > ${SAMPLE}_hg38_reads_to_remove.txt

java -jar /raid_md0/Software/picard.jar FilterSamReads I=${1} READ_LIST_FILE=${SAMPLE}_hg38_reads_to_remove.txt O=${BAM_ROOT}_hg38-removed.bam FILTER=excludeReadList
