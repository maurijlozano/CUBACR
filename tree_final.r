#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


destfile<-paste(args[1],"modal_seqs_NC_GC3s.tab",sep="") 
if (!file.exists(destfile)) {
stop("Nc data missing!", call.=FALSE)
}

#libraries
library('ape')

fn<-paste(args[1],"freqs2.tree",sep='')
tree <- read.tree(file = fn)
fn<-paste(args[1],"CCindex2.txt",sep='')
names <- read.table(fn, sep=" ",stringsAsFactors=F)
fn<-paste(args[1],"modal_seqs_NC_GC3s.tab",sep='')
NCs <- read.table(fn, header=T, sep="\t",stringsAsFactors=F)
NCs <- NCs[,-4]

fn<-paste(args[1],"freqs2.dist",sep='')
dists <- read.table(fn, sep=" ",stringsAsFactors=F, skip=1)
dists <- dists[,6:dim(dists)[2]]
colnames(dists) <- names[,1]
rownames(dists) <- names[,1]
fn<-paste(args[1],"tree.dist.table",sep='')
write.table(dists, fn)

otus <- NULL
for(i in 1:dim(names)[1]){
	if(i < 10){otun <-paste("otu0000",i,sep='')}
	if(i > 9){ 
		if (i < 100){otun <-paste("otu000",i,sep='')}
	} else if (i > 99){otun <-paste("otu00",i,sep='')}
	otus <- c(otus, otun)
}

names <- cbind(otus,names)
names[1:4,2]<-paste(names[1:4,2],"-Nc: ",NCs[1:4,2],sep="")
names[grep('[Ss]ingle',names[,2]),2]<-paste("Single-Nc: ",NCs[grep('[Ss]ingle',names[,2]),2],sep="")

tree$tip.label<-names[[2]][match(tree$tip.label, names[[1]])]

fn<-paste(args[1],"tree2.svg",sep='')
svg(fn,width=5,height=4)
	plot(tree,type="unrooted",no.margin=TRUE,lab4ut="axial", edge.width=2,cex=0.7,label.offset=0.02)
dev.off()


