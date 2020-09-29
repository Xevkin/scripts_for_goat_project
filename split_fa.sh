OUT=`echo $j | cut -f1 -d'.'`


j=0; while read i; do if echo $i | grep -q "^>" ; then j=`expr $j + 1`; fi; echo $j >> ${OUT}_${j}; done < $1
