#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)
n<-as.numeric(args[5])-1

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 


#libraries
library(doParallel)
library(doRNG)
ncores <- detectCores()

merge_data <- function(a,b){
    rbind(a, b)
}

#Function for bootstrap sequence generation
bootstrap_codons <- function(seq){
	bootseq<-sample(seq, length(seq), replace = T)
	colnames(bootseq)<-1:length(bootseq)
	return(bootseq)
}

#Open seq, Parallel calculation of n bootstrapped seqs and save
calculate_bootseqs <- function(fn){
	fn1<-paste(gsub('.concat.fas', '', fn),".boot.fas",sep='')
	seq <- read.table(fn, header=FALSE, sep=" ", stringsAsFactors=FALSE, skip=1, colClasses = c("character"))
	seq <- read.fwf(textConnection(seq[1,1]), widths = rep(3, ceiling(max(nchar(seq[1,1]) / 3))), stringsAsFactors = FALSE)
	bootseq<-NULL
	bootseq<-rbind(bootseq,seq)
	colnames(bootseq)<-1:length(bootseq)
	registerDoParallel(cores=ncores)
	cl <- makeCluster(ncores, type = "FORK")
	registerDoParallel(cl)
	semilla<-runif(1,1,1000)
	registerDoRNG(seed = semilla)
	bootseqs = foreach(i=1:n, .combine=merge_data) %dopar% {bootstrap_codons(seq)}
	stopImplicitCluster()
	bootseq<-rbind(bootseq,bootseqs)
	file.create(fn1)
	for(i in 1:dim(bootseq)[1]){
		write.table(paste('>Seq_',i,sep=''), fn1, row.names = F, col.names = F,append=T,quote=F)
		write.table(cbind("ATG",bootseq[i,],"TAA"), fn1, row.names = F, sep='', col.names = F,append=T,quote=F)
	}
	#return(bootseq)
}


fn <- args[1]
calculate_bootseqs(fn)
fn <- args[2]
calculate_bootseqs(fn)
fn  <- args[3]
calculate_bootseqs(fn)
fn  <- args[4]
calculate_bootseqs(fn)
fn  <- args[6]
calculate_bootseqs(fn)







