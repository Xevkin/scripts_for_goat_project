#provide study accession, sample accession, sample name, coverage

echo STUDY $1\\nSAMPLE $2\\nASSEMBLYNAME ${3}_mtDNA\\nASSEMBLY_TYPE isolate\\nCOVERAGE ${4}\\nPROGRAM bwa_aln'/'ANGSD_doFasta\\nPLATFORM ILLUMINA
