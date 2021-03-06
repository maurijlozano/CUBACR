#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

setwd(args[1])
name <- args[2]
name <- paste(strtrim(name,6),'.',sep='')

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

fdata1 <- fdata[ grep("HEP_CR$|HEP_VR$|LEP_CR$|LEP_VR$|HEP$|LEP$",rownames(fdata)), ]
for (i in 1:ncol(fdata1)){fdata1[,i]<-as.numeric(as.character(fdata1[,i]))}

deltaHCR_HVR <- fdata1[1,]-fdata1[2,]
deltaHCR_LCR <- fdata1[1,]-fdata1[3,] 
deltaHVR_LVR <- fdata1[2,]-fdata1[4,]
deltaLCR_LVR <- fdata1[3,]-fdata1[4,]
deltaHEP_LEP <- fdata1[5,]-fdata1[6,]

rownames(deltaHCR_HVR)="deltaHEP.CR-HEP.VR"
rownames(deltaHCR_LCR)="deltaHEP.CR-LEP.VR"
rownames(deltaHVR_LVR)="deltaHEP.VR-LEP.VR"
rownames(deltaLCR_LVR)="deltaLEP.CR-LEP.VR"
rownames(deltaHEP_LEP)="deltaHEP-LEP"


aminoacids <- c("C","D","E","F","H","K","N","Q","Y","I","A","G","P","T","V","L","R","S")
aminoacids3letter <- c("Cys","Asp","Glu","Phe","His","Lys","Asn","Gln","Tyr","Ile","Ala","Gly","Pro","Thr","Val","Leu","Arg","Ser")
aancod <- c(2,2,2,2,2,2,2,2,2,3,4,4,4,4,4,6,6,6)
aa.table <- as.data.frame(cbind(aminoacids,aminoacids3letter,aancod))
aa.table$aancod <- as.numeric(as.character(aa.table$aancod))

aatable <- gather(deltaHCR_HVR, Codon, Freq)
aatable2 <- gather(deltaHCR_LCR, Codon, Freq)
aatable3 <- gather(deltaHVR_LVR, Codon, Freq)
aatable4 <- gather(deltaLCR_LVR, Codon, Freq)
aatable5 <- gather(deltaHEP_LEP, Codon, Freq)
aatable <- cbind(aatable, DHCR_LCR=aatable2[,2],DHVR_LVR=aatable3[,2],DLCR_LVR=aatable4[,2], DHEP_LEP=aatable5[,2])

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


#Function to calculate olsen distance, x and y are the column numbers of mydata for distance calculation
olsen_dist <- function (x,y){
	odist <- 0
	odist2 <- 0
	for (i in aminoacids){
		codons <- chead[which(chead[,1] == i),2]
		codons <- custom.sort3(codons)[,1]
		for (j in codons){
			mycodonx <- aatable[which(aatable$Codon == codons[which(codons == j)]),x+1]
			mycodony <- aatable[which(aatable$Codon == codons[which(codons == j)]),y+1]
			odist <- odist + (abs(mycodonx-mycodony))
		}
		odist2 <- odist2 + (odist)**2
		odist <- 0
	}
	res <- sqrt(odist2/4)
	return(res)
}


#delta index lopez

i1 <- sum(abs (deltaHCR_HVR - deltaHEP_LEP))
i2 <- sum(abs (deltaHCR_LCR - deltaHEP_LEP))
i3 <- sum(abs (deltaHVR_LVR - deltaHEP_LEP))
i4 <- sum(abs (deltaLCR_LVR - deltaHEP_LEP))
i <- format(c(i1,i2,i3,i4), digits=3)

#euclidian dist
j1 <- dist(t(aatable[,c(2,6)]), method = "euclidean")
j2 <- dist(t(aatable[,c(3,6)]), method = "euclidean")
j3 <- dist(t(aatable[,c(4,6)]), method = "euclidean")
j4 <- dist(t(aatable[,c(5,6)]), method = "euclidean")
j <- format(c(j1,j2,j3,j4), digits=3)

#olsen distance
k1 <- olsen_dist(1,5) 
k2 <- olsen_dist(2,5) 
k3 <- olsen_dist(3,5) 
k4 <- olsen_dist(4,5) 
k <- format(c(k1,k2,k3,k4), digits=3)


heatmap.table <- mydata[,c(8,1,2,3,4,5,6)]
rownames(heatmap.table)=paste(heatmap.table[,1],"|",heatmap.table[,2], sep="")
colnames(heatmap.table)[1]<-"AA"
colnames(heatmap.table)[2]<-"Codon"
colnames(heatmap.table)[3]<-format(paste("HEP_cr-HEP_vr [",k[1],"]",sep=''),width=15)
colnames(heatmap.table)[4]<-format(paste("HEP_cr-LEP_cr [",k[2],"]",sep=''),width=15)
colnames(heatmap.table)[5]<-format(paste("HEP_vr-LEP_vr [",k[3],"]",sep=''),width=15)
colnames(heatmap.table)[6]<-format(paste("LEP_cr-LEP_vr [",k[4],"]",sep=''),width=15)
colnames(heatmap.table)[7]<-"HEP-LEP"
heatmap.tablebk <- heatmap.table

library('pheatmap')

heatmap.table <- heatmap.table[,c(3,6,4,5,7)]



myColor <- c(colorRampPalette(c("firebrick3", "#D13838"))(2),             colorRampPalette(c("#D54A4A", "#DD6E6E"))(2), colorRampPalette(c("#E18080", "#F2C8C8"))(8), colorRampPalette(c("white"))(1), colorRampPalette(c("#D0D0E7", "#7373B9"))(8), colorRampPalette(c("#5C5CAE", "#2E2E97"))(2), colorRampPalette(c("#17178B", "#000080"))(2))

myBreaks <- c(seq(-1, -0.7, length.out=2), seq(-0.65, -0.4, length.out=2), seq(-0.35, -0.05, length.out=8), seq(-0.025, 0.025, length.out=2), seq(0.05, 0.35, length.out=8), seq(0.4, 0.65, length.out=2), seq(0.7, 1, length.out=2))

myColor<-colorRampPalette(c("darkred","firebrick4","firebrick","firebrick3", "white","dodgerblue3","dodgerblue4","navyblue","midnightblue"))(100)
myBreaks <- seq(-1, 1, length.out=101)


#annotation<-data.frame(Strain=rep(name,length(heatmap.table)))
#rownames(annotation) <- colnames(heatmap.table)

svg(filename="delta-heatmap3.svg", width=0.5, height=11, pointsize=10)
p<-pheatmap(heatmap.table, cellwidth = 7, cellheight = 10, cluster_cols = FALSE,cluster_rows = FALSE, border_color=NA,legend = F,color=myColor, breaks=myBreaks,fontsize = 6, fontsize_col = 6, gaps_row =c(2,4,6,8,10,12,14,16,18,21,25,29,33,37,41,47,53), show_rownames=F,main=name)
print(p)
dev.off()
# annotation_names_col = F, annotation_legend = F, annotation_col=annotation, 





       
       
       
       
       
       
