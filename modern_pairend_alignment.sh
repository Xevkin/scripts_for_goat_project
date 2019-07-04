#align all fastq pairs in this directory, using $1 threads

for i in $(ls *1.fastq.gz | cut -f1 -d'_');
do bwa mem -t $1 -R "@RG\tID:$i\tSM:$i\tLB:$i" ~/ARS1/ARS1.fa "$i"_1.fastq.gz "$i"_2.fastq.gz | samtools fixmate -r -@ $1 -O bam - "$i"_fixmate.bam &&
rm "$i"_1.fastq.gz "$i"_2.fastq.gz &&
samtools sort -@ $1 "$i"_fixmate.bam -o "$i"_sorted.bam -T temp &&
rm "$i"_fixmate.bam &&
samtools view -@ 24 -q 30 -b  "$i"_sorted.bam > "$i"_q30.bam &&
rm "$i"_sorted.ba* &&
samtools index "$i"_q30.bam -@ $1 &&
java -Xmx60g -Djava.io.tmpdir=`pwd`/tmp -jar ~/programs/picard/picard.jar MarkDuplicates I="$i"_q30.bam O="$i"_q30_dups-removed.bam METRICS_FILE=metrics_"$i".txt VALIDATION_STRINGENCY=SILENT MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50 MAX_RECORDS_IN_RAM=100000 TMP_DIR=`pwd`/tmp REMOVE_DUPLICATES=true &&
rm  "$i"_q30.ba* &&
samtools index "$i"_q30_dups-removed.bam -@ $1 &&
java -Xmx5g -jar ~/programs/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/ARS1/ARS1.fa -I "$i"_q30_dups-removed.bam -o forIndelRealigner_"$i".intervals -nt $1 2> "$i"_intervals.log &&
java -Xmx60g -jar ~/programs/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ~/ARS1/ARS1.fa -I "$i"_q30_dups-removed.bam -targetIntervals forIndelRealigner_"$i".intervals -o "$i"_q30_dups-removed_realigned.bam &&
rm "$i"_realigned_dups-removed.ba* &&
samtools index "$i"_q30_dups-removed_realigned.bam -@ $1;
done
