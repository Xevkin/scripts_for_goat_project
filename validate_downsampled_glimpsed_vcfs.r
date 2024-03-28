#!/usr/bin/env Rscript

require(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

vcf.true<-read.table(args[1])

vcf.sample<-read.table(args[2])

vcf.true$V11<-str_split_fixed(vcf.true$V10, pattern = ":", n=2)[,1]

vcf.sample$V11<-str_split_fixed(vcf.sample$V10, pattern = ":", n=2)[,1]

vcf.sample$V12<-as.numeric(str_split_fixed(str_split_fixed(vcf.sample$V8, pattern = ";", n=2)[,1],"=",n=2)[,2])

vcf.sample<-data.frame(vcf.sample)

vcf.true<-data.frame(vcf.true)

vcf.merged<-merge(vcf.sample, vcf.true, by=c("V1","V2"))[c("V1","V2","V4.x","V5.x","V12","V11.x","V11.y")]

vcf.merged$sample<-args[3]

vcf.merged$coverage<-as.numeric(args[4])

write.table(x = vcf.merged,file=paste0(args[5]),row.names = F, col.names = F, sep = ",")

