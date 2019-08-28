#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

#Calculate Nc vs GC3s


args = commandArgs(trailingOnly=TRUE)
name <- gsub('[/]*$','',args)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
## program...
library('data.table')
library('naturalsort')
#library('ggpubr')


#functions

setwd(args[1])

#data: FUC vs Cores
#test if file exist

destfile="NC_GC3s.tab"    
if (!file.exists(destfile)) {
    stop("Nc vs GC3s table missing", call.=FALSE)
}

#data 
Ncgenes <- read.table("NC_GC3s.tab", header=F, sep="\t", stringsAsFactors=F)
Ncgenes <- Ncgenes[,1:2]
colnames(Ncgenes) <- c("title","Nc")
Ncgenes[grep('[Gg]enom[ea]',Ncgenes[,1]),1] <- "GNM"
Ncgenes[grep('[Ss]ingle',Ncgenes[,1]),1] <- "Single"

Ncglevels <-c("GNM",as.character(unique(naturalsort(factor(Ncgenes$title[grep('GNM',Ncgenes$title,invert=T)])))))
Ncgenes$title <- factor(Ncgenes$title, levels=Ncglevels)
LC <- length(unique(Ncgenes[grep('^C[0-9]{1,2}',Ncgenes[,1]),1])) 


library(agricolae)
ncgroups <- kruskal(Ncgenes$Nc, Ncgenes$title, group=TRUE, p.adj="holm")
Nc.sum <- data.frame(title=rownames(ncgroups$mean),Nc=ncgroups$mean[,1],SD=ncgroups$mean[,3],SEM=ncgroups$mean[,3]/sqrt(ncgroups$mean[,4]))
Nc.sum <- cbind(SP=rep(name,dim(Nc.sum)[1]),Nc.sum)


write.table(Nc.sum,file="Nc_means.txt")



