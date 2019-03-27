chi.bc <- read.table("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/SegDist/chisq_stats_DP8_oxg.txt",header=TRUE,sep="\t")
chi.bc <- read.table("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/SegDist/chisq_stats_DP8_bc1a.txt",header=TRUE,sep="\t")
chi.bc <- read.table("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/SegDist/chisq_stats_DP8_bc2a.txt",header=TRUE,sep='\t')

chi.bc <- chi.bc[grep("EG9",as.character(chi.bc$CHROM)),]
chi.bc$CHROM <- factor(chi.bc$CHROM)
chi.bc$CHR_MOD <- factor(chi.bc$CHR_MOD)

chi.bc$Pval[chi.bc$Pval==0] <- 2.2e-16
chi.bc.na <- chi.bc[!is.na(chi.bc[,9]),]
bonf.bc <- data.frame(chi.bc.na,bonferroni=p.adjust(chi.bc.na[,9], method = "bonferroni"),fdr=p.adjust(chi.bc.na[,9],method="fdr"))

sig.bc <- bonf.bc[bonf.bc[,15]<0.05,]

library("qqman")
bc.df <- data.frame(SNP=as.character(paste(chi.bc.na[,3],chi.bc.na[,4],sep="_")),CHR=as.numeric(as.character(chi.bc.na[,2])),BP=chi.bc.na[,4]/1000000, P=as.numeric(chi.bc.na[,9])) 

#How many parental heterozygot snps on on each chromosome? or how many snps did we test for seg. distortion on each chromosome?
as.data.frame(table(bc.df$CHR))

manhattan(bc.df)
manhattan(subset(bc.df,CHR ==10.5))