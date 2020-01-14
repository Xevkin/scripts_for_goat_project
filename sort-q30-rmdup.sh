
#input should be the bam file

ROOT=`echo $1 | cut -f1 -d'.'`

samtools flagstat -@ 12 ${ROOT}.bam > ${ROOT}.flagstat

samtools sort -@ 12 ${ROOT}.bam > ${ROOT}_sorted.bam

samtools view -@ 12 -q 30 -b ${ROOT}_sorted.bam > ${ROOT}_sorted_q30.bam

samtools flagstat -@ 12 ${ROOT}_sorted_q30.bam > ${ROOT}_sorted_q30.flagstat

samtools rmdup -s ${ROOT}_sorted_q30.bam ${ROOT}_q30_rmdup.bam

samtools flagstat -@ 12 ${ROOT}_q30_rmdup.bam > ${ROOT}_q30_rmdup.flagstat
