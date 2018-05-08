rm trimlist.sh

for i in $(ls *1.fastq.gz | cut -f1 -d'.'); do

echo cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -m 30 "$i".fastq.gz -o "$i"_trimmed.fastq.gz "2>" "$i"_trimmed.log >> trimlist.sh

done

echo parallel -x $1 -a trimlist.sh

parallel -j $1 -a trimlist.sh

rm trimlist.sh
