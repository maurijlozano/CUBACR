#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

setwd(args[1])


if (!file.exists("modelFreqs.txt")) {
stop("model Freqs missing!", call.=FALSE)
}

table <- read.table("modelFreqs.txt", sep=" ", dec=".", header=F, row.names=1)
table <- rbind(table,TAA=0,TGA=0,TAG=0)
codonOrder <- c("TTT","TTC","TTA","TTG",
				"TCT","TCC","TCA","TCG",
				"TAT","TAC","TAA","TAG",
				"TGT","TGC","TGA","TGG",
				"CTT","CTC","CTA","CTG",
				"CCT","CCC","CCA","CCG",
				"CAT","CAC","CAA","CAG",
				"CGT","CGC","CGA","CGG",
				"ATT","ATC","ATA","ATG",
				"ACT","ACC","ACA","ACG",
				"AAT","AAC","AAA","GAG",
				"AGT","AGC","AGA","GGG",
				"GTT","GTC","GTA","GTG",
				"GCT","GCC","GCA","GCG",
				"GAT","GAC","GAA","GAG",
				"GGT","GGC","GGA","GGG"
)
table2 <- table[match(codonOrder, row.names(table)),]
table3 <- NULL

j <- 1
k <- 1
vec <- NULL
table3 <- NULL

while (j < 5){
	if (k == 65){break}
	vec <- c(vec, table2[k] )
	j<-j+1
	k<-k+1
	if (j == 5){ 
		j <- 1
		table3 <- rbind(table3,vec)
		vec <- NULL
	}
}
#table normalization
table3 <- table3/(sum(table3))

write.table(table3, "modelFreqs.txt", col.names=F, row.names=F,sep=" ",dec=".")
