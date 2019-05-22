#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

setwd(args[1])
name <- args[2]
name <- strtrim(name,5)

library("tidyr")
#load data files

hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")

data = read.table("fmdata.txt", header=F, sep=" ") 
data <- cbind(data, rep(1, length(data[,1])), rep(1,length(data[,1])))
data <- data[,-2]
colnames(data) <- c(as.character(unlist(hdata[,1:60])), "AUG", "UGG")

corder2 <- c("GCU","GCC","GCA","GCG","CGU","CGC","CGA","CGG","AGA","AGG","AAU","AAC","GAU","GAC","UGU","UGC","CAA","CAG","GAA","GAG","GGU","GGC","GGA","GGG","CAU","CAC","AUU","AUC","AUA","UUA","UUG","CUU","CUC","CUA","CUG","AAA","AAG","AUG","UUU","UUC","CCU","CCC","CCA","CCG","UCU","UCC","UCA","UCG","AGU","AGC","ACU","ACC","ACA","ACG","UGG","UAU","UAC","GUU","GUC","GUA","GUG")
t<-t(data)
t <- t[match(corder2,row.names(t)),]
fdata <- t(t)
fdata <- cbind(fdata[,1:37],fdata[,39:54],fdata[,56:61])
fdata <- data.frame(fdata)
rownames(fdata) <- data[,1]


fdata1 <- fdata[ grep("HEP_CR$|HEP$|PHE$",rownames(fdata)), ]
for (i in 1:ncol(fdata1)){fdata1[,i]<-as.numeric(as.character(fdata1[,i]))}


deltaCR_HEP <- fdata1[1,]-fdata1[2,]
deltaCR_PHE <- fdata1[1,]-fdata1[3,] 
deltaHEP_PHE <- fdata1[2,]-fdata1[3,]


rownames(deltaCR_HEP)="deltaCR-HEP"
rownames(deltaCR_PHE)="deltaCR-PHE"
rownames(deltaHEP_PHE)="deltaHEP-PHE"

aminoacids <- c("C","D","E","F","H","K","N","Q","Y","I","A","G","P","T","V","L","R","S")
aminoacids3letter <- c("Cys","Asp","Glu","Phe","His","Lys","Asn","Gln","Tyr","Ile","Ala","Gly","Pro","Trp","Val","Leu","Arg","Ser")
aancod <- c(2,2,2,2,2,2,2,2,2,3,4,4,4,4,4,6,6,6)
aa.table <- as.data.frame(cbind(aminoacids,aminoacids3letter,aancod))
aa.table$aancod <- as.numeric(as.character(aa.table$aancod))

aatable <- gather(deltaCR_HEP, Codon, Freq)
aatable2 <- gather(deltaCR_PHE, Codon, Freq)
aatable3 <- gather(deltaHEP_PHE, Codon, Freq)
aatable <- cbind(aatable, DCR_PHE=aatable2[,2],DHEP_PHE=aatable3[,2])

chead <- read.table("../c_head.txt", header=FALSE, sep=" ", stringsAsFactors = FALSE)
chead <- rbind(chead,c("W","UGG"),c("M","AUG"),c("TER","UAA"),c("TER2","UAG"),c("TER3","UGA"))


custom.sort3 <- function(x){
	data.frame(x)
	thdbase <- substring(x, 3)
	scndbase <- substring(x, 2,2)
	fstbase <- substring(x, 1,1)
	x <- cbind(x,fstbase,scndbase,thdbase)
	ord <- c("U","C","A","G")
	factor(x[,2],levels=ord)
	factor(x[,3],levels=ord)
	x[,4] <- factor(x[,4],levels=ord)
	x[order(x[,2],x[,4]),]
}

mydata <- NULL
for (i in aminoacids){
	codons <- chead[which(chead[,1] == i),2]
	codons <- custom.sort3(codons)[,1]
	for (j in codons){
		mycodons <- aatable[which(aatable$Codon == codons[which(codons == j)]),]
		mycodons <- cbind(mycodons, AA = i, AA3letter=aa.table[which(aa.table$aminoacids == i),2])
		mydata <- rbind(mydata,mycodons)	
	}
}
rownames(mydata)=NULL


heatmap.table <- mydata[,c(6,1,2,3,4)]
rownames(heatmap.table)=paste(heatmap.table[,1],"|",heatmap.table[,2], sep="")
colnames(heatmap.table)[1]<-"AA"
colnames(heatmap.table)[2]<-"Codon"
colnames(heatmap.table)[3]<-"CR-HEP"
colnames(heatmap.table)[4]<-"CR-PHE"
colnames(heatmap.table)[5]<-"HEP-PHE"
heatmap.tablebk <- heatmap.table

library('pheatmap')

heatmap.table <- heatmap.table[,3:(ncol(heatmap.table))]



myColor <- c(colorRampPalette(c("firebrick3", "#D13838"))(2),             colorRampPalette(c("#D54A4A", "#DD6E6E"))(2), colorRampPalette(c("#E18080", "#F2C8C8"))(8), colorRampPalette(c("white"))(1), colorRampPalette(c("#D0D0E7", "#7373B9"))(8), colorRampPalette(c("#5C5CAE", "#2E2E97"))(2), colorRampPalette(c("#17178B", "#000080"))(2))

myBreaks <- c(seq(-1, -0.7, length.out=2), seq(-0.65, -0.4, length.out=2), seq(-0.35, -0.05, length.out=8), seq(-0.025, 0.025, length.out=2), seq(0.05, 0.35, length.out=8), seq(0.4, 0.65, length.out=2), seq(0.7, 1, length.out=2))

#annotation<-data.frame(Strain=rep(name,length(heatmap.table)))
#rownames(annotation) <- colnames(heatmap.table)

svg(filename="delta-heatmap.svg", width=0.3, height=10, pointsize=10)
p<-pheatmap(heatmap.table, cellwidth = 5, cellheight = 10, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = F,color=myColor, breaks=myBreaks,fontsize = 4, fontsize_col = 6, gaps_row =c(2,4,6,8,10,12,14,16,18,21,25,29,33,37,41,47,53), show_rownames=F,main=name,)
print(p)
dev.off()
# annotation_names_col = F, annotation_legend = F, annotation_col=annotation, 




       
       
       
       
       
       
