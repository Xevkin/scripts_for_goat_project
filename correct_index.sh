
rm index_parallel.txt;

rm machine_prefixes.txt; for j in $(ls *fastq.gz *fq.gz | grep -v "exact") ;

	do zcat $j | head -n1 | cut -c2-6 >> machine_prefixes.txt; done;

	sort -u machine_prefixes.txt > tmp; mv tmp machine_prefixes.txt; PREFIX="";

	while read i; do PREFIX=`echo $PREFIX'\|'$i`; done < machine_prefixes.txt; PREFIX=`echo $PREFIX | cut -c3-`;


for j in $(ls *fastq.gz *fq.gz | grep -v "exact") ; do i=`echo $j | cut -f1 -d'.'`; rm ${i}_exact.txt;

	INDEX=`zcat ${j} | head -n 40000| grep $PREFIX |rev | sed -e "s/1\///g" |  cut -c1-8 | rev | sort | uniq -c  | awk -v OFS="\t" '$1=$1' | sort -n -k1,1 | tail -n1 | cut -f2 | sed -e "s/://g"`;

	echo "zcat ${j} | grep \"@*$INDEX\" | cut -f1 -d' ' | cut -f2 -d'@' >> ${i}_exact.txt ; ~/programs/seqtk/seqtk subseq ${j} ${i}_exact.txt | gzip -c - > ${i}_exact.fastq.gz  " >> index_parallel.txt; done

parallel -a index_parallel.txt -j $1

#rm index_parallel.txt
