OUT=`echo $1 | cut -f1 -d'.'`

j=0; while read i; do if echo $i | grep -q "^>" ; then j=`expr $j + 1`; echo $i | sed -e "s/>//g" | cut -f1 -d' ' >> ${OUT}_${j}.fa ; continue; fi; echo $j >> ${OUT}_${j}.fa; done < $1
