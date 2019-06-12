#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

#Calculate Nc vs GC3s


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
## program...
library('ggplot2')
library('ggrepel')
library('ggthemes')
library('data.table')
library('tAI')
library('naturalsort')
#library('ggpubr')


#functions
fgc3s <- function (gc3) 
{
    a = -6
    b = 34
    c = 1.025
    x = a + gc3 + (b/(gc3^2 + (c - gc3)^2))
    return(x)
}




setwd(args[1])

#data: FUC vs Cores
#test if file exist
 
destfile="modal_seqs_NC_GC3s.tab"    
if (!file.exists(destfile)) {
stop("Nc vs GC3s modal table missing", call.=FALSE)
}

destfile="NC_GC3s.tab"    
if (!file.exists(destfile)) {
    stop("Nc vs GC3s table missing", call.=FALSE)
}


#load data
hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")
corder = read.table('../codOrder.txt', header=F, sep=" ")

#data 
Ncmodal <- read.table("modal_seqs_NC_GC3s.tab", header=T, sep="\t")
Ncmodal <- Ncmodal[,-length(Ncmodal)]
Ncmodal[,1] <- gsub('_[^_]*_.*','',Ncmodal[,1])
Ncmodal[grep('[Gg]enom[ea]',Ncmodal[,1]),1] <- "GNM"
Ncmodal[grep('[Ss]ingle',Ncmodal[,1]),1] <- "Single"

Ncmodal$title <- factor(Ncmodal$title, levels=c("GNM",as.character(naturalsort(factor(Ncmodal$title[grep('GNM',Ncmodal$title,invert=T)])))))
LC <- length(Ncmodal[grep('^C[0-9]{1,2}',Ncmodal[,1]),1]) 


svg(filename="Nc_modal.svg", width=11, height=6, pointsize=10)

ggplot(data=Ncmodal, aes(x=title, y=Nc, fill=title)) +
geom_bar(stat="identity", show.legend=F) +
scale_fill_manual(values = c("dodgerblue4",rep("dodgerblue3",LC),'firebrick4','firebrick3','firebrick2', 'maroon4','maroon3','maroon2','gold3','aquamarine3')) +
#coord_cartesian(ylim = c(yminimo,ymaximo)) + #ylim(yminimo, ymaximo) hace desaparecer las barras
ylab("Nc") + 
xlab("Set") +
theme_classic() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_blank(),
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()

#nc.adj(nc, gc3)


ncplot <- data.frame(x=seq(0,1,0.01), y=fgc3s(seq(0,1,0.01)))



svg(filename="Nc_modal_vs_GC3s.svg", width=8, height=6, pointsize=10)

ggplot(data=Ncmodal, aes(x=GC3s, y=Nc)) +
geom_line(data=ncplot, aes(x=x,y=y))+
geom_point(aes(colour=title),show.legend=T) +
labs(color = "Gene Set") +
scale_color_manual(values = c("dodgerblue4",rep("dodgerblue3",LC),'firebrick4','firebrick3','firebrick2', 'maroon4','maroon3','maroon2','gold3','aquamarine3')) +
#coord_cartesian(ylim = c(yminimo,ymaximo)) + #ylim(yminimo, ymaximo) hace desaparecer las barras
ylab("Nc") + 
xlab("GC3s") +
theme_classic() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"), #element_blank()
	  axis.title.y = element_text(size=18, face="bold"),
) 

dev.off()







Ncgenes <- read.table("NC_GC3s.tab", header=T, sep="\t", stringsAsFactors=F)
Ncgenes <- Ncgenes[,-length(Ncgenes)]
colnames(Ncgenes) <- colnames(Ncmodal)
Ncgenes[grep('[Gg]enom[ea]',Ncgenes[,1]),1] <- "GNM"
Ncgenes[grep('[Ss]ingle',Ncgenes[,1]),1] <- "Single"

Ncglevels <-c("GNM",as.character(unique(naturalsort(factor(Ncgenes$title[grep('GNM',Ncgenes$title,invert=T)])))))
Ncgenes$title <- factor(Ncgenes$title, levels=Ncglevels)
LC <- length(unique(Ncgenes[grep('^C[0-9]{1,2}',Ncgenes[,1]),1])) 


#complist <- function(table){
#	mycomp <- NULL
#	k<-1
#	tlevels <-c("GNM",as.character(unique(naturalsort(factor(table[grep('GNM',table[,1],invert=T),1])))))
#	mycomb <- combinations(length(tlevels), 2,tlevels)
#	npc <- pairwise.wilcox.test(table[,2],table[,1], p.adjust.method = "holm")
#	for (i in 1:(length(tlevels)-1)){
#	for (j in 1:(length(tlevels)-1)){
#		pvalue <- npc$p.value[i,j]
#		if(is.na(pvalue)){} else if (pvalue < 0.05){ 
#		mycomp[[k]] <- c(rownames(npc$p.value)[i],colnames(npc$p.value)[j])
#		k<-k+1 
#		}
#	}
#	}
#	return(mycomp)
#}
#
#my_comparisons <- complist(Ncgenes)

library(agricolae)
ncgroups <- kruskal(Ncgenes$Nc, Ncgenes$title, group=TRUE, p.adj="holm")
Nc.sum <- data.frame(title=rownames(ncgroups$mean),Nc=ncgroups$mean[,6])
Nc.sum <- naturalsort(Nc.sum )
ncgroups <- ncgroups$groups[naturalsort(rownames(ncgroups$groups)),]
ncstats <- kruskal(Ncgenes$Nc, Ncgenes$title, group=F, p.adj="holm")


svg(filename="Nc_genes.svg", width=11, height=6, pointsize=10)

ggplot(data=Ncgenes, aes(x=title, y=Nc, fill=title)) +
geom_boxplot(outlier.colour="gray55", outlier.shape=8, outlier.size=1, show.legend=F) +
scale_fill_manual(values = c("dodgerblue4",rep("dodgerblue3",LC),'firebrick4','firebrick3','firebrick2', 'maroon4','maroon3','maroon2','gold3','aquamarine3')) +
#stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels for ggpubr
geom_text(data=Nc.sum ,aes(x=Nc.sum$title,y=(Nc.sum$Nc+(Nc.sum$Nc*.03)),label=ncgroups$groups),vjust=0)+
ylab("Nc") + 
xlab("Set") +
theme_classic() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()





ncplot <- data.frame(x=seq(0,1,0.01), y=fgc3s(seq(0,1,0.01)))
svg(filename="Nc_genes_vs_GC3s.svg", width=8, height=6, pointsize=10)

ggplot(data=Ncgenes, aes(x=GC3s, y=Nc)) +
geom_line(data=ncplot, aes(x=x,y=y))+
geom_point(aes(colour=title)) +
labs(color = "Gene Set") +
scale_color_manual(values = c("dodgerblue4",rep("dodgerblue3",LC),'firebrick4','firebrick3','firebrick2', 'maroon4','maroon3','maroon2','gold3','aquamarine3')) +
#coord_cartesian(ylim = c(yminimo,ymaximo)) + #ylim(yminimo, ymaximo) hace desaparecer las barras
ylab("Nc") + 
xlab("GC3s") +
theme_classic() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"), #element_blank()
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()



