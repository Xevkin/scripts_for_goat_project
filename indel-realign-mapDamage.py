import sys

from subprocess import call

bamroot = sys.argv[1]

ref = sys.argv[2]

call("java -Xmx100g -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 24 -R " + ref + " -I " + bamroot + ".bam -o forIndelRealigner_" + bamroot + ".intervals 2> " + bamroot + "_intervals.log",shell=True)

call("java -Xmx100g -jar /home/kdaly/programs/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R " + ref + " -I " + bamroot + ".bam -targetIntervals forIndelRealigner_"  + bamroot + ".intervals -o " + bamroot + "_realigned.bam  2> "+ bamroot + "_realigned.log",shell=True)

call("samtools index -@ 24 " + bamroot + "_realigned.bam",shell=True)

call("mapDamage -r " + ref + " -i " + bamroot + "_realigned.bam " + "--rescale --merge-reference-sequences",shell=True)
