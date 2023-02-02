import sys
import subprocess

FASTQ=sys.argv[1]

SAMPLE = FASTQ.split(".")[0]

REFPATH=sys.argv[2]

REF=REFPATH.split("/")[-1].split(".")[0]

subprocess.call("bwa aln -l 1000 -n 0.01 -o 2 -t 10 " + REFPATH + " " + FASTQ + " > " + SAMPLE + "_" + REF + ".sai" ,shell=True)

subprocess.call("bwa samse " + REFPATH + " " + SAMPLE + "_" + REF + ".sai " + FASTQ + " | samtools view -F4 -q30 -Sb - > " + SAMPLE + "_" + REF + "_F4-q30.bam"  ,shell=True)

subprocess.call("samtools sort "  + SAMPLE + "_" + REF + "_F4-q30.bam -o "  + SAMPLE + "_" + REF + "_F4-q30_sort.bam",shell=True)

subprocess.call("samtools rmdup -s " + SAMPLE + "_" + REF + "_F4-q30_sort.bam " + SAMPLE + "_" + REF + "_F4-q30_rmdup.bam && rm " + SAMPLE + "_" + REF + "_F4-q30_sort.bam " + SAMPLE + "_" + REF + ".sai", shell=True)

subprocess.call("for i in $(ls " + SAMPLE + "*bam | cut -f1 -d'.'); do samtools flagstat -@ 24 ${i}.bam > ${i}.flagstat; done ",shell=True)
