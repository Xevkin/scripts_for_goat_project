#script adds read groups to MG2 bams with different indexes, then merges the bams, sorts and rmdups, and generates F4 and q25 bams
java -jar /research/picard-tools-1.119/AddOrReplaceReadGroups.jar \
INPUT=MG2_index15.bam \
OUTPUT=MG2_index15_RG.bam \
RGID=MG2_1st_hiseq_index15 \
RGLB=MG2 \
RGPL=illumina_hiseq \
RGSM=MG2 \
RGPU=unknown_platform_unit
java -jar /research/picard-tools-1.119/AddOrReplaceReadGroups.jar \
INPUT=MG2_index16.bam \
OUTPUT=MG2_index16_RG.bam \
RGID=MG2_1st_hiseq_index16 \
RGLB=MG2 \
RGPL=illumina_hiseq \
RGSM=MG2 \
RGPU=unknown_platform_unit
java -jar /research/picard-tools-1.119/AddOrReplaceReadGroups.jar \
INPUT=MG2_index17.bam \
OUTPUT=MG2_index17_RG.bam \
RGID=MG2_1st_hiseq_index17 \
RGLB=MG2 \
RGPL=illumina_hiseq \
RGSM=MG2 \
RGPU=unknown_platform_unit
java -jar /research/picard-tools-1.119/AddOrReplaceReadGroups.jar \
INPUT=MG2_index18.bam \
OUTPUT=MG2_index18_RG.bam \
RGID=MG2_1st_hiseq_index18 \
RGLB=MG2 \
RGPL=illumina_hiseq \
RGSM=MG2 \
RGPU=unknown_platform_unit
java -jar /research/picard-tools-1.119/MergeSamFiles.jar \
INPUT=MG2_index15_RG.bam \
INPUT=MG2_index16_RG.bam \
INPUT=MG2_index17_RG.bam \
INPUT=MG2_index18_RG.bam \
OUTPUT=MG2_hs1.bam
gzip MG2_index15.bam MG2_index16.bam MG2_index17.bam MG2_index18.bam
samtools sort MG2_hs1.bam MG2_hs1_sort
samtools rmdup -s MG2_hs1_sort.bam MG2_hs1_rmdup.bam
rm MG2_hs1_sort.bam
samtools view -b -F 4 MG2_hs1_rmdup.bam > MG2_hs1_F4.bam
samtools view -b -q25 MG2_hs1_F4.bam > MG2_hs1_q25.bam
gzip MG2_hs1_rmdup.bam MG2_hs1_F4.bam MG2_hs1_q25.bam
