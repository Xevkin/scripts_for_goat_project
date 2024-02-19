args = commandArgs(trailingOnly=TRUE)

sample<-read.table(args[1],header=T,sep=",")

lineages<-read.table(args[2],header=T,sep=" ")

colnames(lineages)<-c("CHROM", "POS", "REF.lin", "ALT.lin","LIN")

merged<-merge(sample, lineages)

merged
