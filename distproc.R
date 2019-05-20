#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)
n<-as.numeric(args[1])
nset<-as.numeric(args[2])
distmat <- read.table("distances", header=FALSE, row.names=1 ,sep=" ", stringsAsFactors=FALSE, skip=1, colClasses = c("character"))
distmat<-distmat[,-c(1,2,3,4)]

#index begin with 0
distance_significance <- function(i,j){
	if(i==0){N1<-"HEP_CR"}
	if(i==1){N1<-"HEP_VR"}
	if(i==2){N1<-"LEP_CR"}
	if(i==3){N1<-"LEP_VR"}
	if(i==4){N1<-"PHE"}
	if(j==0){N2<-"HEP_CR"}
	if(j==1){N2<-"HEP_VR"}
	if(j==2){N2<-"LEP_CR"}
	if(j==3){N2<-"LEP_VR"}
	if(j==4){N2<-"PHE"}
	submat <- distmat[((i*n)+1):((i+1)*n),((j*n)+1):((j+1)*n)]
	real.dist <- as.numeric(submat[1,1])
	submat <- as.numeric(submat[upper.tri(submat,diag = F)])
	submat.mean <- mean(submat)
	submat.SEM <- sd(submat)/sqrt(length(submat))
	if(i==j){
		submat1 <- submat[1:(length(submat)/2)]
		submat2 <- submat[((length(submat)/2)+1):length(submat)]
		submat.stats <- wilcox.test(submat1,submat2, alternative = "two.sided")
		stable<-data.frame(Seq.SETs=paste(N1,'.vs.', N2,sep=''), Dist=real.dist, Mean=submat.mean, SEM=submat.SEM, p.value.AB.AB=submat.stats$p.value, p.value.A.AB="NA", p.value.B.AB="NA", stringsAsFactors=FALSE)
		return(stable)
	} else {
		submat1 <- submat[1:(length(submat)/2)]
		submat2 <- submat[((length(submat)/2)+1):length(submat)]
		submat3 <- distmat[((i*n)+1):((i+1)*n),((i*n)+1):((i+1)*n)]
		submat3 <- as.numeric(submat3[upper.tri(submat3,diag = F)])
		submat4 <- distmat[((j*n)+1):((j+1)*n),((j*n)+1):((j+1)*n)]
		submat4 <- as.numeric(submat4[upper.tri(submat4,diag = F)])
		submat.stats.AB.AB <- wilcox.test(submat1,submat2, alternative = "two.sided")
		submat.stats.A.AB <- wilcox.test(submat1,submat3, alternative = "two.sided")
		submat.stats.B.AB <- wilcox.test(submat1,submat4, alternative = "two.sided")
		stable<-data.frame(Seq.SETs=paste(N1,'.vs.', N2,sep=''), Dist=real.dist, Mean=submat.mean, SEM=submat.SEM, p.value.AB.AB=submat.stats.AB.AB$p.value, p.value.A.AB=submat.stats.A.AB$p.value, p.value.B.AB=submat.stats.B.AB$p.value,stringsAsFactors=FALSE)
		return(stable)
	}	
}

Final.table <- NULL
k<-0
for (i in 0:(nset-1)){
	for (j in k:(nset-1)){
		stable<-distance_significance(i,j)
		Final.table <- rbind(Final.table,stable)
	}
	k<-k+1
}

write.table(Final.table, "dist.table", row.names = F, sep=' ', col.names = T,quote=F)


