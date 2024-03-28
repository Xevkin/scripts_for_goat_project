"""
April 2019 - modified Victoria's script
To mask base qualities of bam files.
How to run:
samtoool view -h samfile | python script - [bases to clip]| samtools view -Sb - > file.bam
To do different amounts of masking change the number of #### and then [?:-?]
"""

import sys
import fileinput

bases_to_softclip = int(sys.argv[2])

mask = ""

for i in range(0,bases_to_softclip):

	mask = mask + "#"

ct = 0

for line in fileinput.input(sys.argv[1]):
   if line[0:1] == '@':
       print line.strip()
   else:
       test = line.strip().split('\t')
       if(len(mask+test[10][bases_to_softclip:-bases_to_softclip]+mask) != len(test[10])):

            raise SystemExit("new qual line not equal to old " + test[10] + " " +  mask+test[10][bases_to_softclip:-bases_to_softclip]+mask)
       else:
           print '\t'.join(test[0:10])+'\t'+mask+test[10][bases_to_softclip:-bases_to_softclip]+mask+'\t'+'\t'.join(test[11:])

