#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

## program...
library('ggplot2')
library('ggrepel')
library('ggthemes')

setwd(args[1])
data = read.table(args[2], header=TRUE, sep=" ") #check field separator

genes <- data[ which(data$label=='Gene'), ]
modal <- data[ grep("C[0-9]{1,2}$|^[Ss]ingle|^[gG]enom|^[Gg]NM",data$label), ]
HEP <- data[ grep("^HEP",data$label), ]
LEP <- data[ grep("^LEP",data$label), ]
PHE <- data[ grep("^PHE",data$label), ]
labels <- data[ grep("^Gene|^[gG]enom|^[Gg]NM",data$label, invert=T), ]
labels[,1]<-as.character(labels[,1])
labels[grep("[Ss]ingle",labels$label),1]<-"Single"

varpc = read.table("coa_var_axis1_2.txt", header=TRUE, sep=" ") #check field separator
varpc <- varpc*100
xlabtext <- paste("C1 (",varpc[1,1]," %)", sep="")
ylabtext <- paste("C2 (",varpc[1,2]," %)", sep="")


png(filename = "CA_genes.png", width = 1000, height = 800, units = "px", pointsize = 9,
    bg = "white", res = 230, type = "cairo-png")

ggplot (data, aes( y = Axis2, x = Axis1)) +
geom_point(data = genes, color='gray', size = 1, alpha = I(0.3)) +
geom_point(data = modal, fill = "dodgerblue3", size=3, colour="black", shape=21, stroke = 0.5) +
geom_point(data = LEP, fill = "darkslategray3", size=3, colour="black", shape=21, stroke = 0.5) +
geom_point(data = HEP, fill = "firebrick3", size=3, colour="black", shape=21, stroke = 0.5) +
geom_point(data = PHE, fill = "dodgerblue3", size=3, colour="black", shape=21, stroke = 0.5) +
xlab(xlabtext) + ylab(ylabtext) + labs(fill = "Legend")  +
geom_text_repel(data = labels, aes(label=label, fontface=2), size = 3, point.padding = 0.1 , arrow = arrow(length = unit(0.01, 'npc'))) +
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 

dev.off()

