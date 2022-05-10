rm ~/raid/papers/direkli4_paper2021/analyses/damage.txt; 

for i in $(ls -d *mapDamage); do I=`echo $i | cut -f1 -d'_'`; 

for BASES in 1 2 3 4 5 6 7 8 9 10; 

do for DIRECT in 5p 3p; 

do grep $DIRECT ${i}/misincorporation.txt | cut -f3,5,10-12  |  awk -v var=${BASES} '{ if ($2 == var) print}' | awk -v varDIR=${DIRECT} -v varSAM=${I} -v varBASE=${BASES}  '{sum1+=$3; sum2+=$4; sum3+=$5;} END{print varSAM" "varDIR" "varBASE" "sum1" "sum2" "sum3;}' >> ~/raid/papers/direkli4_paper2021/analyses/damage.txt; 

done; 

done; 

done
