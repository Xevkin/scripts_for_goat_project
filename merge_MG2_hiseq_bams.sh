#this script is to merge the MG2 sample bams following meyer alignment
#bam is then indexed and DepthOfCoverage is run - large file produced is then removed
java -jar /research/picard-tools-1.119/MergeSamFiles.jar \
INPUT=../../MG2_2hs_q25_RG.bam \
INPUT=../../first_hiseqs/MG2/MG2_hs1_q25.bam \
INPUT=/d3/giant_deer_data/MG2_deer_q25_HC_meyer_F4_fixed_RGs.bam \
OUTPUT=MG2_merged_q25.bam
samtools index MG2_merged_q25.bam 
java -jar /research/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T DepthOfCoverage -R /deer_wg/deer_HS_data/deer_ref/deer_ref_no_U.fa -I=MG2_merged_q25.bam -o=MG2_merged_DoC
rm MG2_merged_DoC
