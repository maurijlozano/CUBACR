#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

#Calculate tAI vs dist on R using complete genome


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

setwd(args[1])

#data: FUC vs Cores
#test if file exist
destfile="modal.counts"    
if (!file.exists(destfile)) {
stop("Counts data file does not exist", call.=FALSE)
}

destfile="s_opts_DCBS_GNM.txt"    
if (!file.exists(destfile)) {
    stop("S_BCBS data file does not exist", call.=FALSE)
}


#load data
hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")
sopt_DCBS <- read.table("s_opts_DCBS_GNM.txt", header=TRUE, sep="\t", dec=".")
sopt_DCBS <- sopt_DCBS[order(sopt_DCBS[,10], decreasing=T),]
corder = read.table('../codOrder.txt', header=F, sep=" ")

#data for tai calcualtion - modals


f.counts <- read.table("modal.counts", header=FALSE, row.names = 1, sep=" ")
colnames(f.counts) <- c(as.character(unlist(hdata[,2:60])), "AUG", "UGG")
corder2 <- corder[-c(11,12,15),]
t<-t(f.counts)
t <- t[match(corder2$V2,row.names(t)),]
fdata <- t(t)


#trnas
trna = read.table('trna.txt', header=T, row.names = NULL, sep=" ")
trna <- trna[match(corder$V2,trna$codon),]
t<-t(fdata)
t <- t[match(corder2$V2,row.names(t)),]
f.o <- t(t)


library('tAI')

Sopt_DCBS <- as.numeric(sopt_DCBS[1,1:9])

if (trna$trna[35] == 0){sk <- 1 } else {sk <- 0}
ws_DCBS <- get.ws(tRNA=trna$trna,s=Sopt_DCBS, sking=sk)
tai_DCBS <- get.tai(f.o[,-33], ws_DCBS)


fdata <- as.data.frame(cbind(SET=as.character(rownames(f.counts)),tai_DCBS))
rownames(fdata)[rownames(fdata) %like% "^single"] = "Single"
rownames(fdata)[rownames(fdata) %like% "LEP_CR$"] = "LEP_CR"
rownames(fdata)[rownames(fdata) %like% "HEP_CR$"] = "HEP_CR"
rownames(fdata)[rownames(fdata) %like% "LEP_VR$"] = "LEP_VR"
rownames(fdata)[rownames(fdata) %like% "HEP_VR$"] = "HEP_VR"
rownames(fdata)[rownames(fdata) %like% "^PHE"] = "PHE"
fdata[,1] <- rownames(fdata)
fdata$SET <- factor(fdata$SET, levels = c("Single","LEP","LEP_CR","LEP_VR","PHE","HEP","HEP_CR","HEP_VR"))
fdata$tai_DCBS <- as.numeric(as.character(fdata$tai_DCBS))

#plots
yminimo <- min(fdata$tai_DCBS) - (max(fdata$tai_DCBS)/10)
ymaximo <- max(fdata$tai_DCBS) + (max(fdata$tai_DCBS)/10)


svg(filename="tAi_CR_VR.svg", width=8, height=6, pointsize=10)

ggplot(data=fdata, aes(x=SET, y=tai_DCBS, fill=SET)) +
geom_bar(stat="identity", show.legend=F)+
scale_fill_manual(values=c('aquamarine3','dodgerblue3','dodgerblue1','dodgerblue2','gold3','firebrick1','firebrick3','firebrick2'))+
coord_cartesian(ylim = c(yminimo,ymaximo)) + #ylim(yminimo, ymaximo) hace desaparecer las barras
ylab("m-tAI") + 
theme_calc() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_text(size=18, face="bold"),
	  axis.title.y = element_text(size=18, face="bold")
) 

dev.off()



