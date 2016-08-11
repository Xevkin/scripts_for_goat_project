#!/usr/bin/R

###### Companion script for plotting the results of PMDtools (Skoglund, Northoff, Shunkov, Dervianko, Paabo, Krause, Jakobsson, 2014, PNAS)
###### Contact: pontus.skoglund@gmail.com 
###### 
###### Usage example: 
######        samtools view mybam.bam | python pmdtools.py --deamination > PMD_temp.txt
######        R CMD BATCH plotPMD.R
######        cp PMD_plot.pdf mybam_PMD_plot.pdf



pdf("PMD_plot.pdf",width=12,height=5)


par(mfrow=c(1,2))


data<-read.table("PMD_temp.txt",header= TRUE)


mismatch_CC<-data$CC
mismatch_CT<-data$CT
mismatch_CA<-data$CA
mismatch_CG<-data$CG

mismatch_GC<-data$GC
mismatch_GT<-data$GT
mismatch_GA<-data$GA
mismatch_GG<-data$GG

read_position <-data$z

plot(read_position,mismatch_CC,xlab="Distance from 5' end of sequence read",ylab="Base frequency in read given C in reference",cex=3,type="n",ylim=c(0,0.4),xlim=c(0,max(read_position)))
lines(read_position,mismatch_CC,col="black",lwd=2)
lines(read_position,mismatch_CT,col="blue",lwd=2)
lines(read_position,mismatch_CG,col="green",lwd=2)
lines(read_position,mismatch_CA,col="red",lwd=2)
legend("topright",c("T","A","G","C"),lty=c(1,1,1,1),lwd=c(2,2,2,2),col=c("blue","red","green","black"),bty="b",cex=1,bg="white")

plot(read_position,mismatch_GG,xlab="Distance from 3' end of sequence read",ylab="Base frequency in read given G in reference",cex=3,type="n",ylim=c(0,0.4),xlim=c(max(read_position),0))
lines(read_position,mismatch_GC,col="black",lwd=2)
lines(read_position,mismatch_GT,col="blue",lwd=2)
lines(read_position,mismatch_GG,col="green",lwd=2)
lines(read_position,mismatch_GA,col="red",lwd=2)
legend("topleft",c("T","A","G","C"),lty=c(1,1,1,1),lwd=c(2,2,2,2),col=c("blue","red","green","black"),bty="b",cex=1,bg="white")

dev.off()
