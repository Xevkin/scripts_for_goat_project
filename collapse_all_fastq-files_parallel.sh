rm collapselist.sh

for R2 in $(ls *R2.fastq.gz *_2.fastq.gz);

do

R1=`echo $R2 | sed -e "s/_R2.fastq/_R1.fastq/g" | sed -e "s/-R2.fastq/-R1.fastq/g" | sed -e "s/_2.fastq/_1.fastq/g"`

ROOT=`echo $R1 | sed -e "s/_R1.fastq.gz//g" | sed -e "s/-R1.fastq.gz//g"`

echo AdapterRemoval --threads 2 --collapse --minadapteroverlap 1 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minlength 30 --gzip --trimns --trimqualities  \
--file1 $R1 --file2 $R2 --basename $ROOT >> collapselist.sh

done

parallel -a collapselist.sh -j 10
