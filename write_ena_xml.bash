SAMPLE=`echo ${1^} | cut -f1 -d'_'`;

ERS=`echo $SAMPLE | sed -f /raid_md0/kdaly/papers/hare2022/name_to_ERS.sed -`

echo '<'?xml version=\"1.0\" encoding=\"US-ASCII\"?'>'
echo '<'ANALYSIS_SET'>'
echo -e \\t'<'ANALYSIS alias=\"${SAMPLE}_oviAri3_alignment\"'>'
echo -e \\t\\t'<'TITLE'>'Alignment of ancient sheep sample to oviAri3 assembly'</'TITLE'>'
echo -e \\t\\t'<'DESCRIPTION'>'Sequence alignment of sample ${SAMPLE} to oviAri3, presented in bam format after duplicate removal and mapQ30 filtering.'</'DESCRIPTION'>'
echo -e \\t\\t'<'STUDY_REF accession=\"ERP146890\"'/>'
echo -e \\t\\t'<'SAMPLE_REF accession=\"$ERS\"'/>'
grep ${ERS} /raid_md0/kdaly/papers/hare2022/runs-2023-07-06.csv > tmp.csv
while read i; do I=`echo $i | cut -f1 -d',' `; echo -e \\t\\t'<RUN_REF accession='\"$I\"'/>'; done < tmp.csv
echo -e \\t\\t'<'ANALYSIS_TYPE'>'
echo -e \\t\\t\\t'<'REFERENCE_ALIGNMENT'>'
echo -e \\t\\t\\t\\t'<'ASSEMBLY'>'
echo -e \\t\\t\\t\\t\\t'<'STANDARD accession=\"GCA_000298735.1\"'/>'
echo -e \\t\\t\\t\\t'</'ASSEMBLY'>'
echo -e \\t\\t\\t\\t'<'SEQUENCE accession=\"CM001582.1\"'/>'
echo -e \\t\\t\\t\\t'<'SEQUENCE accession=\"CM001583.1\"'/>'
echo -e \\t\\t\\t'</'REFERENCE_ALIGNMENT'>'
echo -e \\t\\t'</'ANALYSIS_TYPE'>'
echo -e \\t\\t'<'FILES'>'
echo -e \\t\\t\\t'<'FILE checksum=\"`grep $1 md5sum.upload | cut -f1 -d' '`\" checksum_method=\"MD5\" filename=\"${1}\" filetype=\"bam\"'/>'
echo -e \\t\\t'</'FILES'>'
echo -e \\t'<'/ANALYSIS'>'
echo '</'ANALYSIS_SET'>'
