#expects samples_SXX_XX[_-]R[12].fastq.gz

rm counts.sh counts.out

rm fastqc.sh

rm adaptor.sh

rm fastqc_trim.sh

rm fastqscreen.sh

mkdir fastq_screen

mkdir fastqc

for R1 in $(ls *R1.fastq.gz)

do

R2=`echo $R1 |sed -e "s/R1.fastq/R2.fastq/g"`

ROOT=`echo $R1 | cut -f1 -d'_'`

echo J='`'zcat ${R1} "|" wc -l "|" awk \'{print '$0 / 4' }\' '`' "&&" echo $R1 $J ">>" counts.out >> counts.sh

echo fastqc $R1 $R2 >> fastqc.sh

echo AdapterRemoval  --threads 2 --collapse --minadapteroverlap 1 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minlength 30 --gzip --trimns --trimqualities  --file1 $R1 --file2 $R2 -- basename $ROOT >> adaptor.sh

echo J='`'zcat ${ROOT}.collapsed.gz "|" wc -l "|" awk \'{print '$0 / 4' }\' '`' "&&" echo ${ROOT}.collapsed.gz $J ">>" counts.out >> counts.sh

echo fastqc ${ROOT}.collapsed.gz >> fastqc_trim.sh

echo fastqscreen --conf ~/fastq_screen.conf --aligner bowtie --force --outdir ./fastq_screen/ ${ROOT}.collapsed.gz >> fastqscreen.sh

done

parallel -a fastqc.sh -j ${1} && mv *html fastqc && parallel -a adaptor.sh -j ${1} && parallel -a fastqc_trim.sh -j $1 &&  parallel -a fastqscreen.sh -j $1 && parallel -a counts.sh -j $1

