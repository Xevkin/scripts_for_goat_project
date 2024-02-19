args = commandArgs(trailingOnly=TRUE)

sample<-read.table(args[1])

colnames(sample)<-c("CHROM","POS","REF","OBSERVED")

lineages<-read.table(args[2],header=T,sep=" ")

colnames(lineages)<-c("CHROM", "POS", "REF.lin", "ALT.lin","LIN")

merged<-merge(sample, lineages)

merged.not.match<-merged[merged$OBSERVED != merged$LIN,]

#print(merged.match)

for (x in 1:nrow(merged.not.match)) {

	a<-merged.not.match[x,]

	if (a$LIN == a$REF.lin) {

	print(paste(a$CHROM, a$POS, a$POS, paste0(a$REF.lin,"/",a$ALT.lin),"+",sep=" "))
	} else {

	print(paste(a$CHROM, a$POS, a$POS, paste0(a$ALT.lin,"/",a$REF.lin),"+",sep=" "))
}

}
