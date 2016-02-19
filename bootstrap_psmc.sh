

PSMC_FILE=$1
PSMCFA_FILE=$2
GEN_TIME=$3
MUTATION_RATE=$4
PSMC="/home/kdaly/bin/psmc/"

$PSMC"/bin/utils/splitfa" $PSMCFA_FILE > split.psmcfa

seq 100 | xargs -P 4 -i echo $PSMC"psmc" -N25 -t15 -r5 -b -p "4+50*1+4+6" -o round-{}.psmc split.psmcfa | sh 

cat $PSMC_FILE round-*.psmc > combined.psmc

echo $PSMC"utils/psmc_plot.pl" -u $MUTATION_RATE -g $GEN_TIME -M '"'$SAMPLE_NAME"="$FNR'"' $SAMPLE_NAME"-"$FNR $SAMPLE_NAME".psmc" >tmp.sh
	
