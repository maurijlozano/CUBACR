#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

#libraries
library('pheatmap')
library('ape')
library('phytools')

#load data files
fn<-paste(args[1],"fmdata.txt",sep='')
hdata <- read.table("freq_head.txt", header=FALSE, sep=" ")
data = read.table(fn, header=F, sep=" ",row.names=1)
data <- data[,-1]
data <- cbind(data, rep(1, length(data[,1])), rep(1,length(data[,1])))
colnames(data) <- c(as.character(unlist(hdata[,2:60])), "AUG", "UGG")

#filename = "pheheatmap.pdf", cellwidth = 25, cellheight = 15, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = FALSE, color=colorRampPalette(c("white", "navy"))(50))

#rename tree tips
fn<-paste(args[1],"freqs.tree",sep='')
tree <- read.tree(file = fn)
fn<-paste(args[1],"CCindex.txt",sep='')
names <- read.table(fn, sep=" ",stringsAsFactors=F)

otus <- NULL
for(i in 1:dim(names)[1]){
	if(i < 10){otun <-paste("otu0000",i,sep='')}
	if(i > 9){ 
		if (i < 100){otun <-paste("otu000",i,sep='')}
	} else if (i > 99){otun <-paste("otu00",i,sep='')}
	otus <- c(otus, otun)
}

names <- cbind(otus,names)
tree$tip.label<-names[[2]][match(tree$tip.label, names[[1]])]

fn<-paste(args[1],"renamed.tree",sep='')
write.tree(tree,file = fn)

colors<-colorRampPalette(colors=c("red","white","blue"))(100)

#plot(tree)
#identify(tree, nodes = TRUE, tips = FALSE,labels = FALSE, quiet = FALSE)
#nodelabels() -> adds node labels to plot
#rotate(tree, node)
fn<-paste(args[1],"tree_heetmap.svg",sep='')

svg(fn,width=14,height=10)
	phylo.heatmap(tree,data[,-c(60,61)],colors=colors,fsize=c(1,0.7,1),standardize=T)
dev.off()


