BAMFILE=${1}
#paths must end in /
EXTRA_PILEUP_PATH=${2}
EXTRA_PILEUP_PATH2=${3}
EXTRA_PILEUP_PATH3=${4}

OUTNAME=`echo $BAMFILE | cut -f1 -d'.'`
DATASET_PATH=/raid_md0/goat/vargoats/phased_june2022/
REF=/raid_md0/Reference_Genomes/goat/ARS1.fa
LASER=/home/kdaly/.local/bin/laser
PILEUPS_PATH=/raid_md0/goat/vargoats/phased_june2022/pileups/published/
LASER_PILEUP_PY=/home/kdaly/programs/LASER/pileup2seq/pileup2seq.py
RAND=`echo $RANDOM`; rm parallel${RAND}.sh

while read i

do

SAMPLE=`echo $i | rev | cut -f1 -d'/' | rev | cut -f1 -d'_'`

for DATASET in vargoats_snps_filtered_1372_beagle5-3_phased_combined_trans_maf05

do

DATASETNAME="vargoats-trans-maf05"

echo samtools mpileup -q 30 -Q 20 -B -f $REF -l ${DATASET_PATH}${DATASET}.bed ${i} ">" ${SAMPLE}_${DATASETNAME}.pileup >> parallel${RAND}.sh

done

done < $BAMFILE

parallel -a  parallel${RAND}.sh -j 12

rm  parallel${RAND}.sh

for DATASET in vargoats_snps_filtered_1372_beagle5-3_phased_combined_trans_maf05

do

python2 $LASER_PILEUP_PY -f ${REF} -m ${DATASET_PATH}${DATASET}.site -b ${DATASET_PATH}${DATASET}.bed -o ${DATASETNAME}_${OUTNAME} *_${DATASETNAME}.pileup ${PILEUPS_PATH}*_${DATASETNAME}.pileup ${EXTRA_PILEUP_PATH}*_${DATASETNAME}.pileup ${EXTRA_PILEUP_PATH2}*_${DATASETNAME}.pileup ${EXTRA_PILEUP_PATH3}*_${DATASETNAME}.pileup

echo -e GENO_FILE ${DATASET_PATH}${DATASET}.geno\\nNUMTHREADS 8\\nSEQ_FILE ${DATASETNAME}_${OUTNAME}.seq\\nOUT_PREFIX ${DATASETNAME}_${OUTNAME}\\nDIM 6\\nMIN_LOCI 2500 > ${DATASETNAME}_${OUTNAME}.par

${LASER} -p ${DATASETNAME}_${OUTNAME}.par -r 10

sed -f /home/kdaly/programs/scripts_for_goat_project/goat-id_to_name.sed ${DATASETNAME}_${OUTNAME}.SeqPC.coord | sed -e "s/_${DATASETNAME}//g" > tmp1 && mv tmp1  ${DATASETNAME}_${OUTNAME}.SeqPC.coord

#Rscript /home/kdaly/raid/kiel/programs/sheep_laser_pca.r ${DATASETNAME}_${OUTNAME}.RefPC.coord  ${DATASETNAME}_${OUTNAME}.SeqPC.coord

done
