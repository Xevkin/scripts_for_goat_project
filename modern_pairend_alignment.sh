#align all fastq pairs in this directory, using $1 threads and $2 mapQ
#$3 indicates the path of the reference genome, $4 the reference name to attach to the files


for i in $(ls *1.fastq.gz | cut -f1 -d'_');
do bwa mem -t $1 -R "@RG\tID:$i\tSM:$i\tLB:$i" $3 ${i}_1.fastq.gz ${i}_2.fastq.gz | samtools fixmate -r -@ $1 -O bam - ${i}_${4}_fixmate.bam &&
rm ${i}_1.fastq.gz ${i}_2.fastq.gz &&
samtools sort -@ $1 ${i}_${4}_fixmate.bam -o ${i}_${4}_sorted.bam -T temp &&
rm ${i}_${4}_fixmate.bam &&
samtools view -@ 24 -q ${2} -b  ${i}_${4}_sorted.bam > ${i}_${4}_q${2}.bam &&
rm ${i}_${4}_sorted.ba* &&
samtools index ${i}_${4}_q${2}.bam -@ $1 &&
java -Xmx60g -Djava.io.tmpdir=`pwd`/tmp -jar ~/programs/picard/picard.jar MarkDuplicates I=${i}_${4}_q${2}.bam O=${i}_${4}_q${2}_dups-removed.bam METRICS_FILE=metrics_"$i".txt VALIDATION_STRINGENCY=SILENT MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50 MAX_RECORDS_IN_RAM=100000 TMP_DIR=`pwd`/tmp REMOVE_DUPLICATES=true &&
rm  ${i}_${4}_q${2}.ba* &&
samtools index ${i}_${4}_q${2}_dups-removed.bam -@ $1 &&
java -Xmx5g -jar ~/programs/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/ARS1/ARS1.fa -I ${i}_${4}_q${2}_dups-removed.bam -o forIndelRealigner_"$i".intervals -nt $1 2> "$i"_intervals.log &&
java -Xmx60g -jar ~/programs/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ~/ARS1/ARS1.fa -I ${i}_${4}_q${2}_dups-removed.bam -targetIntervals forIndelRealigner_"$i".intervals -o ${i}_${4}_q${2}_dups-removed_realigned.bam &&
rm ${i}_${4}_q${2}_dups-removed.bam &&
samtools index ${i}_${4}_q${2}_dups-removed_realigned.bam -@ $1;
done
