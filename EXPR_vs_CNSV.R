#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R


destfile="Expr.table"    
if (!file.exists(destfile)) {
stop("Expression table named Expr.table required. The table must contain the following tab separated fields: Gene PA SET[C:core, IC:Intermediate core, AC:Ancestral core] Species", call.=FALSE)
}

exprtable <- read.table("Expr.table", header=T, sep="\t", stringsAsFactors=FALSE)
constable <- read.table("CNSV.table", header=T, sep=" ", stringsAsFactors=FALSE)

figtable <- cbind(exprtable[match(constable[,2],exprtable[,1]),],constable[,c(1,3:5)])
figtable <- figtable[,-3]
colnames(figtable) <- c("Gene","PA","SP","EL","AA100","C100","AAV")

library('ggplot2')
library('ggrepel')
library('ggthemes')



#gsub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

# "Expression.Protein.file" "AA.ID" "CODON.ID" "AA.VAR" "AA.Length"      


svg(filename="EXPRvsCNSV.svg", width=9, height=7, pointsize=10)

	ggplot (figtable, aes( x = log(PA), y = AA100, color=SP, shape=EL)) +
	geom_point(size=2) +
	theme_classic() +
	xlab("log(PA)") +
	ylab("% Conservation") +
	theme(axis.text.x = element_text(size=8),
		  axis.text.y = element_text(size=8),
		  axis.title.x = element_text(size=10, face="bold"),
		  axis.title.y = element_text(size=10, face="bold")
	)
	#facet_wrap(SP ~ ., ncol=4, scales="free")
	#changed facet.grid for facet_wrap

dev.off()



















