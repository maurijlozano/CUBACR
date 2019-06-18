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
figtable <- figtable[which(! is.na(figtable[,1])),]


library('ggpmisc')
library('ggplot2')
library('ggrepel')
library('ggthemes')



#gsub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)

# "Expression.Protein.file" "AA.ID" "CODON.ID" "AA.VAR" "AA.Length"      
my.formula <- y ~ x


svg(filename="EXPRvsCNSV.svg", width=18, height=7, pointsize=10)

	ggplot (figtable, aes( x = log(PA), y = AA100, color=SP)) +
	geom_point(size=2, show.legend=F) +
	geom_smooth(method = "lm", se = T, colour="black",formula = my.formula)+
	stat_poly_eq(formula = my.formula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE, colour="black") +
	theme_classic() +
	xlab("log(PA)") +
	ylab("% Conservation") +
	theme(axis.text.x = element_text(size=8),
		  axis.text.y = element_text(size=8),
		  axis.title.x = element_text(size=10, face="bold"),
		  axis.title.y = element_text(size=10, face="bold")
	)+
	facet_wrap(SP ~ ., ncol=4, scales="free")
	#changed facet.grid for facet_wrap

dev.off()



















