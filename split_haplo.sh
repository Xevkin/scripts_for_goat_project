WINDOW=`echo $2`

WINDOW_MINUS=`expr $WINDOW - 1`

START=1

LENGTH=`tail -n1 ${1}.haplo | cut -f2`

rm ${1}_${START}-${WINDOW}.haplo

rm split_parallel.sh

while [ $START -le $LENGTH ]

do

	END=`expr $START + $WINDOW_MINUS`

	rm ${1}_${START}-${END}.haplo

	head -n1 ${1}.haplo > ${1}_${START}-${END}.haplo

	echo awk  \'{if '($2 >=' $START '&& $2 <' $END')' print}\' ${1}.haplo ">>" ${1}_${START}-${END}.haplo >> split_parallel.sh

	END=`expr $END + 1`

	START=`echo $END`

done

parallel -a split_parallel.sh -j 10
