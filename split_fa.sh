OUT=`echo $1 | cut -f1 -d'.'`

j=0; while read i; do if echo $i | grep -q "^>" ; then j=`expr $j + 1`; NAME=`echo $i | sed -e "s/>//g" |cut -f1 -d' '`; echo $i | cut -f1 -d' ' >> ${OUT}_${NAME}.fa;  continue; fi; echo $i >> ${OUT}_${NAME}.fa; done < $1
