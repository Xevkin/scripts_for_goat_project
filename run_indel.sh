#Genome should be $1, bam $2 , threads $3

REF=$1
BAMROOT=`echo $2 | cut -f1 -d'.'`
THREADS=$3

samtools index -@ $THREADS "$BAMROOT".bam

java -Xmx20g -jar /home/kdaly/programs/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $THREADS -R $REF -I "$BAMROOT".bam -o forIndelRealigner_"$BAMROOT".intervals 2> "$BAMROOT"_intervals.log

java -Xmx100g -jar /home/kdaly/programs/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I "$BAMROOT".bam -targetIntervals forIndelRealigner_"$BAMROOT".intervals -o "$BAMROOT"_realigned.bam 2> "$BAMROOT"_realignment.log
