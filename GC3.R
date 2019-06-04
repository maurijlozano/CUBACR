#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

#Calcular el indice en R


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

hdata <- read.table("freq_head.txt", header=FALSE, sep=" ")

setwd(args[1])

#data: FUC vs Cores
#test if file exist
destfile="fcdata.txt"    
if (!file.exists(destfile)) {
stop("Frequence data file does not exist", call.=FALSE)
}


#Calculate %GC on the third condon base: NN@

fcdata <- read.table("fcdata.txt", header=FALSE, sep=" ")
colnames(fcdata) <- as.character(unlist(hdata[1,]))
fcdata[,1] <- as.character(fcdata[,1])
fcdata[fcdata[,1] %like% "[Ss]ingle",1] <- rep("Single",length(fcdata[fcdata[,1] %like% "[Ss]ingle",1]))
fcdata[fcdata[,1] %like% "[Pp][Hh][Ee]",1] <- rep("PHE",length(fcdata[fcdata[,1] %like% "[Pp][Hh][Ee]",1]))

fcdata[fcdata[,1] %like% "HEP_CR",1] <- rep("HEP_CR",length(fcdata[fcdata[,1] %like% "HEP_CR",1]))
fcdata[fcdata[,1] %like% "LEP_CR",1] <- rep("LEP_CR",length(fcdata[fcdata[,1] %like% "LEP_CR",1]))
fcdata[fcdata[,1] %like% "HEP_VR",1] <- rep("HEP_VR",length(fcdata[fcdata[,1] %like% "HEP_VR",1]))
fcdata[fcdata[,1] %like% "LEP_VR",1] <- rep("LEP_VR",length(fcdata[fcdata[,1] %like% "LEP_VR",1]))



GC3 <- fcdata[grep("..G|..C",names(fcdata))]
fcdata.sum <- cbind.data.frame(SET=fcdata[,"SET"], total = apply(as.matrix(fcdata[,-1]),1, sum))
GC3.sum <- cbind.data.frame(SET=fcdata[,"SET"], GC3 = apply(as.matrix(GC3),1,sum))

GC3.res <- cbind.data.frame(SET=fcdata[,"SET"], GC3 = ((GC3.sum[,-1])/(fcdata.sum[,-1]))*100 )
GC3.res$SET <- factor(GC3.res$SET, levels = c("Single","LEP","LEP_CR","LEP_VR","PHE","HEP","HEP_CR","HEP_VR"))
GC3.res$GC3 <- as.numeric(as.character(GC3.res$GC3))


#plots
yminimo <- min(GC3.res$GC3) - (max(GC3.res$GC3)/10)
ymaximo <- max(GC3.res$GC3) + (max(GC3.res$GC3)/10)



svg(filename="GC3_CR_VR.svg", width=8, height=6, pointsize=10)

ggplot(data=GC3.res, aes(x=SET, y=GC3, fill=SET)) +
geom_bar(stat="identity", show.legend=F)+
scale_fill_manual(values=c('aquamarine3','dodgerblue3','dodgerblue1','dodgerblue2','gold3','firebrick1','firebrick3','firebrick2'))+
coord_cartesian(ylim = c(yminimo,ymaximo)) + #ylim(yminimo, ymaximo) hace desaparecer las barras
ylab("GC3") + xlab("")+
theme_classic() +
theme(axis.text.x = element_text(size=14),
	  axis.text.y = element_text(size=14),
	  axis.title.x = element_blank(),
	  axis.title.y = element_text(size=18, face="bold")
) 


dev.off()

