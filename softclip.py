"""
Feb 2017
To mask base qualities of bam files.
How to run:
samtoool view -h samfile | script | samtools view -Sb - > file.bam
To do different amounts of masking change the number of #### and then [?:-?]
"""

import fileinput

ct = 0

for line in fileinput.input():
   if line[0:1] == '@':
       print line.strip()
   else:
       test = line.strip().split('\t')
       if(len('####'+test[10][4:-4]+'####') != len(test[10])):

            raise SystemExit("new qual line not equal to old")
       else:
           print '\t'.join(test[0:10])+'\t'+'####'+test[10][4:-4]+'####'+'\t'+'\t'.join(test[11:])

