rm collapselist.sh

for R1 in $(ls *R1.fastq.gz);

do

R2=`echo $R1 | sed -e "s/_R1.fastq/_R2.fastq/g" | sed -e "s/-R1.fastq/-R2.fastq/g"`

ROOT=`echo $R1 | sed -e "s/_R1.fastq.gz//g" | sed -e "s/-R1.fastq.gz//g"`

echo AdapterRemoval --threads 2 --collapse --minadapteroverlap 1 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minlength 30 --gzip --trimns --trimqualities  \
--file1 $R1 --file2 $R2 --basename $ROOT >> collapselist.sh

done

parallel -a collapselist.sh -j 10
