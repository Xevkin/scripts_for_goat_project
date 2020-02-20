#!/bin/bash

#input must not be gzip

#goes from VCF to several filtered ped files

MODERN=`echo $1`

RANDOMA=`echo $RANDOM`

#for VCF in $(ls ${MODERN}.vcf); do gzip $VCF; done

rm parallel${RANDOMA}.txt ;

#for MISS in 0 05; do

#	INVMISS=`echo "1 - 0.${MISS}"|bc`;

#	vcftools --gzvcf ${MODERN}.vcf.gz --max-missing ${INVMISS} --recode --recode-INFO-all --out ${MODERN}_missing${MISS} &&

#	mv ${MODERN}_missing${MISS}.recode.vcf ${MODERN}_missing${MISS}.vcf &&

#	gzip -f ${MODERN}_missing${MISS}.vcf &&

#	i=`echo ${MODERN}_missing${MISS}` &&

#	for FILT in 01 02 025 03 05; do

#		INPUT="" &&

#		INFILT=$(bc <<< "1 - 0.${FILT}") &&

#		INPUT=`echo $INPUT "&&" vcftools  --gzvcf ${i}.vcf.gz --recode --recode-INFO-all --maf 0.${FILT}  --max-maf ${INFILT} --out ${i}_MAF${FILT} "&&" mv ${i}_MAF${FILT}.recode.vcf ${i}_MAF${FILT}.vcf` &&

#		echo $INPUT | sed -e "s/^&& //g" >> parallel${RANDOMA}.txt;

#	done

#done

parallel -j 20 -a parallel${RANDOMA}.txt && \

rm parallel${RANDOMA}.txt; rm  CpG${RANDOMA}.txt

for i in $(ls ${MODERN}*MAF*vcf | grep -v "CpG\|Trans" |cut -f1 -d'.'); do

	echo $i

	grep -v "#" ${i}.vcf | awk '{print $1"\t"$2-1"\t"$2}' > ${i}.bed && \

	bgzip -f ${i}.bed && tabix -p bed ${i}.bed.gz && \

	gunzip -c ${i}.bed.gz > ${i}.bed && \

	vcftools --vcf ${i}.vcf --plink --out $i && \

	rm CpG${RANDOMA}.txt; \

	while read A B C D; do E=`expr ${D} - 1` && F=`expr ${D} + 1` && \

	echo G=\`samtools faidx ~/ARS1/ARS1.fa ${A}:${E}-${F}\` "&&" echo '$G' "|" grep \"GC\\\|CG\" "|" cut -f1 -d \' \' "|" sed -e \"s/\>//g\" "|" sed -e \"s/-/ /g\" "|" sed -e \"s/:/ /g\" "|" awk \'{print '$1'\":\"'$2'+1}\' ">>" CpG${RANDOMA}.txt >> parallel${RANDOMA}.txt;  done < ${i}.map && \

	parallel -j 40 -a parallel${RANDOMA}.txt && rm parallel${RANDOMA}.txt && \

	plink --noweb --cow --file ${i} --extract CpG${RANDOMA}.txt --recode --out ${i}_CpG && rm CpG${RANDOMA}.txt && \

	plink --noweb --cow --file ${i}_CpG --missing --out ${i}_CpG && plink --file ${i}_CpG --indep-pairwise 50 5 0.5 --noweb --cow --out ${i}_CpG && \

	plink --noweb --cow --file ${i}_CpG --extract ${i}_CpG.prune.in --recode --out ${i}_CpG_LD50 && \

	plink --noweb --cow --file ${i}_CpG_LD50 --missing --out ${i}_CpG_LD50 && plink --file $i --indep-pairwise 50 5 0.5 --noweb --cow --out $i && \

	plink --noweb --cow --file $i --extract ${i}.prune.in --recode --out ${i}_LD50 && plink --noweb --cow --file ${i}_LD50 --missing --out ${i}_LD50 && j=`echo ${i}_LD50` && \

	grep -v "#" ${i}.vcf | awk '{if ( !( $4 == "C" && $5 == "T" ) ) print }' | awk '{if ( !( $4 == "T" && $5 == "C" ) ) print }' | awk '{if ( !( $4 == "G" && $5 == "A" ) ) print }' | awk '{if ( !( $4 == "A" && $5 == "G" ) ) print }' > tmp.txt && cat ${i}.vcf | grep "#" > header2_filt.txt && cat header2_filt.txt tmp.txt >"$i"_Trans.vcf && \

	vcftools --vcf ${i}_Trans.vcf --plink --out ${i}_Trans &&

	plink --noweb --cow --file ${i}_Trans --missing --out ${i}_Trans && plink --file ${i}_Trans --indep-pairwise 50 5 0.5 --noweb --cow --out ${i}_Trans && \

	plink --noweb --cow --file ${i}_Trans --extract ${i}_Trans.prune.in --recode --out ${i}_Trans_LD50 && \

	plink --noweb --cow --file ${i}_Trans_LD50 --missing --out ${i}_Trans_LD50 && \

	for j in $(ls ${i}*map | cut -f1 -d'.'); do

		cat ${j}.map | awk '{print $1"\t"$4-1"\t"$4}' > ${j}.bed && bgzip -f ${j}.bed && tabix -p bed ${j}.bed.gz && zcat ${j}.bed.gz > ${j}.bed ;

	done

done
