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

dataU <- data[c(2,13,14,16,30,44,45,47),1:3]
dataC <- data[c(5,7,8,19,33,35,50),1:3]

varpc = read.table("coa_var_axis1_2.txt", header=TRUE, sep=" ") #check field separator
varpc <- varpc*100
xlabtext <- paste("C1 (",varpc[1,1]," %)", sep="")
ylabtext <- paste("C2 (",varpc[1,2]," %)", sep="")

png(filename = "CA_codons.png", width = 1000, height = 800, units = "px", pointsize = 12,
    bg = "white", res = 230, type = "cairo-png")

ggplot (data, aes( y = Axis2, x = Axis1)) +
geom_point(fill = "grey65", size=3, colour="black", shape=21, stroke = 0.5) +
geom_point(data=dataU, fill = "#61b8ffff", size=3, colour="black", shape=21, stroke = 0.5) +
geom_point(data=dataC, fill = "#8a0000ff", size=3, colour="black", shape=21, stroke = 0.5) +
xlab(xlabtext) + ylab(ylabtext) + labs(fill = "Legend")  +
geom_text_repel(aes(label=label, fontface=2), size = 2, point.padding = .1 , arrow = arrow(length = unit(0.01, 'npc'))) +
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 

dev.off()

