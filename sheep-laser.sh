BAMFILE=$1
OUTNAME=`echo $BAMFILE | cut -f1 -d'.'`
DATASET_PATH=/st1/hdd/pg/kdaly/sheep/datasets/
REF=/Reference_Genomes/sheep/oviAri3_mod.fa

RAND=`echo $RANDOM`; rm parallel${RAND}.sh
while read i
do
SAMPLE=`echo $i | rev | cut -f1 -d'/' | rev | cut -f1 -d'_'`
for DATASET in andrew andrew_trans
do
echo samtools mpileup -q 30 -Q 20 -B -f $REF -l ${DATASET_PATH}/${DATASET}.bed ${i} ">" ${SAMPLE}_${DATASET}.pileup >> parallel${RAND}.sh
done
done < $BAMFILE

parallel -a  parallel${RAND}.sh -j 12

rm  parallel${RAND}.sh

for DATASET in andrew andrew_trans
do
python /home/kdaly/programs/LASER-2.04/pileup2seq/pileup2seq.py -f ${REF} -m ${DATASET_PATH}/${DATASET}.site -b ${DATASET_PATH}/${DATASET}.bed -o ${DATASET}_${OUTNAME} *_${DATASET}.pileup
done
