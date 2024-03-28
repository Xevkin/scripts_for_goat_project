library(dplyr)

args = commandArgs(trailingOnly=TRUE)

output.file<-as.character(args[3])

data.to.filter <- read.table(args[1],header=T)

data.to.filter.col<-colnames(data.to.filter)

data.to.filter.col[1] <- "chr"

data.to.filter.col[2] <- "position"

colnames(data.to.filter)<-data.to.filter.col

filter.by <- read.table(args[2],header=F)

colnames(filter.by) <- c("chr", "position")

filtered.data<-anti_join(data.to.filter, filter.by, by=c("chr","position"))

write.table(x=filtered.data, file=output.file, sep = "\t", row.names = F, append = F, quote = F)
