#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)
#args <- c("","")
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 


table <- args[1]
dir <- sub('[^/]*$','',args[1])


table <- read.table(table, header=T, sep=" ", stringsAsFactors=FALSE, colClasses = c("character"))
table[,2] <- as.numeric(table[,2])
table[,3] <- as.numeric(table[,3])
table[,4] <- as.numeric(table[,4])
table[,5] <- as.numeric(table[,5])

library('ggplot2')
library('ggrepel')
library('ggthemes')

#gsub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

# "Expression.Protein.file" "AA.ID" "CODON.ID" "AA.VAR" "AA.Length"      
new.col<-dim(table)[2]+1
for (i in 1:dim(table)[1]){
	table[i,(new.col)] <- gsub('(^[^/]*).*','\\1',table[i,1])

}
colnames(table)[6]<-"SET"

table[,6]<-factor(table[,6])
new.col<-new.col+1
for (i in 1:dim(table)[1]){
	table[i,(new.col)] <- table[i,2]*table[i,5]/100

}
colnames(table)[7]<-"Length.AA.100"

new.col<-new.col+1
for (i in 1:dim(table)[1]){
	table[i,(new.col)] <- table[i,3]*table[i,5]/100

}
colnames(table)[8]<-"Length.NT.100"

new.col<-new.col+1
for (i in 1:dim(table)[1]){
	table[i,(new.col)] <- table[i,4]*table[i,5]/100

}
colnames(table)[9]<-"Length.AA.Var"


#more libraries
library(plyr)
library(tidyr)



long_table <- gather(table, SEQ.ID, ID_Length, c(2,3,4,5,7,8,9))
attach(long_table)


long_table$SEQ.ID = factor(long_table$SEQ.ID, c("AA.ID","CODON.ID", "AA.VAR", "AA.Length", "Length.AA.100", "Length.NT.100","Length.AA.Var"))

mu <- ddply(long_table, SET ~ SEQ.ID, summarise, grp.mean=mean(ID_Length))
  
  
fn <- paste(dir,"coverage_plot.svg",sep='')



svg(filename=fn, width=9, height=4, pointsize=10)

	ggplot (long_table, aes( x = ID_Length, fill=SET)) +
	geom_density(alpha=0.4) +
	geom_vline(data=mu, aes(xintercept=grp.mean, color=SET), linetype="dashed") +
	geom_text(data = mu, aes(x=grp.mean, y = median(density(long_table$ID_Length)$y)), color="black", label = round(mu$grp.mean,2), hjust = 0, size = 2, angle = 90, fontface=2) +
	theme_calc() +
	xlab("% ID / Sequence Length") +
	ylab("Density") +
	theme(axis.text.x = element_text(size=8),
		  axis.text.y = element_text(size=8),
		  axis.title.x = element_text(size=10, face="bold"),
		  axis.title.y = element_text(size=10, face="bold")
	) +
	facet_wrap(SEQ.ID ~ ., ncol=4, scales="free")
	#changed facet.grid for facet_wrap

dev.off()



















