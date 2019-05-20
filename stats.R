#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)
#args <- c("","")
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

aa.ali <- args[1]
nt.ali <- args[2]

#load data files

aa.table <- read.table(aa.ali, header=FALSE, sep=" ", stringsAsFactors=FALSE, colClasses = c("character"))

nt.table <- read.table(nt.ali, header=FALSE, sep=" ", stringsAsFactors=FALSE, colClasses = c("character"))

aminoacids <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X","-")
codons <-c("TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG","ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG","TCT","TCC", "TCA","TCG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG","TAT","TAC","TAA", "TAG","CAT","CAC","CAA","CAG","AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG","TGT","TGC","TGA","TGG", "CGT","CGC","CGA","CGG","AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG","---")
aa <- c("F","F","L","L","L","L","L","L","I","I","I","M","V","V","V","V","S","S","S","S","P","P","P","P", "T","T","T","T","A","A","A","A","Y","Y","X","X","H","H","Q","Q","N","N","K","K","D","D","E","E","C","C", "X","W","R","R","R","R","S","S","R","R","G","G","G","G","-")
aa_cod <- cbind(AA=aa,Codons=codons)



#AA 100%C
aa.freq <- data.frame(cbind(AA=aminoacids, Freq=rep(0,length(aminoacids))))
aa.freq[,2]<-as.numeric(as.character(aa.freq[,2]))
aa.freqs.pos <- cbind(AA=aminoacids)

for (j in 2:dim(aa.table)[2]){
	for (i in 1:dim(aa.table)[1]){
		aa.freq[which(aa.freq[,1]==aa.table[i,j]),2] <- aa.freq[which(aa.freq[,1]==aa.table[i,j]),2] +1 	
	}
	aa.freqs.pos <- cbind(aa.freqs.pos, aa.freq[,2]/sum(aa.freq[,2]))
	aa.freq[,2]<-0
}
colnames(aa.freqs.pos)[2:dim(aa.freqs.pos)[2]]<-1:(dim(aa.freqs.pos)[2]-1)

idx<-NULL
for (j in 2:dim(aa.freqs.pos)[2]){
	for (i in 1:dim(aa.freqs.pos)[1]){
		if(aa.freqs.pos[i,j]==1){idx<-c(idx, j)}
	}
}

nt.aa100C <- nt.table[1,idx]
nt.aa100C <- nt.aa100C[which(nt.aa100C != "---")]
nt.aa100C.p <- length(idx)/dim(nt.table)[2]*100


#Codon 100%C
nt.freq <- data.frame(cbind(aa_cod[,2], Freq=rep(0,length(aa_cod[,1]))))
nt.freq[,2]<-as.numeric(as.character(nt.freq[,2]))
nt.freqs.pos <- cbind(Codon=aa_cod[,2])

for (j in 2:dim(nt.table)[2]){
	for (i in 1:dim(nt.table)[1]){
		nt.freq[which(nt.freq[,1]==nt.table[i,j]),2] <- nt.freq[which(nt.freq[,1]==nt.table[i,j]),2] +1 	
	}
	nt.freqs.pos <- cbind(nt.freqs.pos, nt.freq[,2]/sum(nt.freq[,2]))
	nt.freq[,2]<-0
}
colnames(nt.freqs.pos)[2:dim(nt.freqs.pos)[2]]<-1:(dim(nt.freqs.pos)[2]-1)

idx<-NULL
for (j in 2:dim(nt.freqs.pos)[2]){
	for (i in 1:dim(nt.freqs.pos)[1]){
		if(nt.freqs.pos[i,j]==1){idx<-c(idx, j)}
	}
}

nt.nt100C <- nt.table[1,idx]
nt.nt100C <- nt.nt100C[which(nt.nt100C != "---")]
nt.nt100C.p <- length(idx)/dim(nt.table)[2]*100


#Codon Variable
idx<-NULL
for (j in 2:dim(aa.freqs.pos)[2]){
	if(mean(as.numeric(as.character(aa.freqs.pos[which(aa.freqs.pos[,j]!=0),j])))<0.2){idx<-c(idx, j)}
}
nt.aaV <- nt.table[1,idx]
nt.aaV <- nt.aaV[which(nt.aaV != "---")]
nt.aaV.p <- length(idx)/dim(nt.table)[2]*100



#write fasta files
fn <- gsub('.aa.table', '', args[1])
fn1<-paste(fn,'.aa100C.fas',sep='')

file.create(fn1)
if(length(nt.aa100C)==0){
	write.table(paste('>',nt.table[1,1],sep=''), fn1, , row.names = F, col.names = F,append=T,quote=F)
	write.table(paste('ATG','TAA',sep=''), fn1, , row.names = F, sep='', col.names = F,append=T,quote=F)
} else {
	if(nt.aa100C[1]=="ATG"){nt.aa100C<-nt.aa100C[-1]}
	if(length(nt.aa100C)!=0){
		if(nt.aa100C[length(nt.aa100C)]=="TAA" || nt.aa100C[length(nt.aa100C)]=="TGA" || nt.aa100C[length(nt.aa100C)]=="TGA"){nt.aa100C <- nt.aa100C[-length(nt.aa100C)]}
	}
	write.table(paste('>',nt.table[1,1],sep=''), fn1, , row.names = F, col.names = F,append=T,quote=F)
	write.table(cbind("ATG",nt.aa100C,"TAA"), fn1, , row.names = F, sep='', col.names = F,append=T,quote=F)
}



fn1<-paste(fn,'.nt100C.fas',sep='')

file.create(fn1)
if(length(nt.nt100C)==0){
	write.table(paste('>',nt.table[1,1],sep=''), fn1, row.names = F, col.names = F,append=T,quote=F)
	write.table(paste('ATG','TAA',sep=''), fn1, row.names = F, sep='', col.names = F,append=T,quote=F)
} else {
	if(nt.nt100C[1]=="ATG"){nt.nt100C<-nt.nt100C[-1]}
	if(length(nt.nt100C)!=0){
		if(nt.nt100C[length(nt.nt100C)]=="TAA" || nt.nt100C[length(nt.nt100C)]=="TGA" || nt.nt100C[length(nt.nt100C)]=="TGA"){nt.nt100C <- nt.nt100C[-length(nt.nt100C)]}
	}
	write.table(paste('>',nt.table[1,1],sep=''), fn1, row.names = F, col.names = F,append=T,quote=F)
	write.table(cbind("ATG",nt.nt100C,"TAA"), fn1, row.names = F, sep='', col.names = F,append=T,quote=F)
}



fn1<-paste(fn,'.aaV.fas',sep='')

file.create(fn1)
if(length(nt.aaV)==0){
	write.table(paste('>',nt.table[1,1],sep=''), fn1, row.names = F, col.names = F,append=T,quote=F)
	write.table(paste('ATG','TAA',sep=''), fn1, row.names = F, sep='', col.names = F,append=T,quote=F)
} else {
	if(nt.aaV[1]=="ATG"){nt.aaV<-nt.aaV[-1]}
	if(length(nt.aaV)!=0){
		if(nt.aaV[length(nt.aaV)]=="TAA" || nt.aaV[length(nt.aaV)]=="TGA" || nt.aaV[length(nt.aaV)]=="TGA"){nt.aaV <- nt.aaV[-length(nt.aaV)]}
	}
	write.table(paste('>',nt.table[1,1],sep=''), fn1, row.names = F, col.names = F,append=T,quote=F)
	write.table(cbind("ATG",nt.aaV,"TAA"), fn1, row.names = F, sep='', col.names = F,append=T,quote=F)
}

fp <- args[3]
fn1<-paste(fp,'.stats',sep='')
conservation.table<-cbind(AA_100C=nt.aa100C.p, Codon_100C=nt.nt100C.p, AA_V=nt.aaV.p)
write.table(conservation.table,fn1,  row.names = F, col.names = F, quote=F)












