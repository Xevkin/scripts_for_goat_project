#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(ape)

#1=matrix
#2=sample list
#3=root
#4=output

ibs<-read.table(args[1])
rownames(ibs)<-read.table(args[2])$V1
colnames(ibs)<-read.table(args[2])$V1
boot.d<-as.dist(ibs)
boot.tree<-root(nj(boot.d),args[3])
write.tree(boot.tree, file=args[4])
