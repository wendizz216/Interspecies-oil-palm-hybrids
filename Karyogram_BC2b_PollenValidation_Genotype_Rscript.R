#################################################
#Karyotype figures annotate genetic regions corresponding to wild and cultivated oil palm genes as well as regions identified to harbor 
#Karyotype file was generated using a python script karyotype_file.py in the seg_dist directory
#################################################

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

library(GenomicRanges)
library(ggbio)

data(darned_hg19_subset500, package = "biovizBase")
dn <- darned_hg19_subset500
library(GenomicRanges)


###Karyogram by chromosomes
kary.snp <- read.csv("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/SegDist/Karyogram/bc2a_specSNPs.csv",header=T)

kary.snp <- kary.snp[,colnames(kary.snp) %in% c('chr_idx','chr','start','end','type','SD_Pval')]
kary.snp <- kary.snp[!(is.na(kary.snp[,5])),]
kary.snp$type <- as.factor(kary.snp$type)


colnames(kary.snp)[1] <- "chr"
colnames(kary.snp)[2] <- "chr1"

#reorder factor levels for chromosomes
sort <- unique(kary.snp$chr)
kary.snp$chr <- factor(kary.snp$chr, levels=as.vector(sort))


gr.chrom <- makeGRangesFromDataFrame(kary.snp,keep.extra.columns=T)
autoplot(gr.chrom, layout="karyogram",aes(color=type,fill=type)) 


#MANUALLY CHANGING COLORS
https://support.bioconductor.org/p/77482/
#data frame must have a separate column of colors that correspond to the factors of interest before generating a genome range meta file.
?scale_fill_identity()

autoplot(gr.kary3, layout="karyogram",aes(color=color,fill=color)) + scale_fill_identity("Genotype", labels = gr.kary3$geno, breaks = gr.kary3$color, guide = "legend") + scale_color_identity()


####Identify palms to prioritize for breeders
#Filter genotypes for pollen functional validation
pollen_genos <- read.table("/Users/wvu/Dropbox/Oil_Palm/RAD_SequencingResults_Dec2017/SegDist/Karyogram/pollen_validation_genos.txt",header=F)

kary.snp <- kary.snp[!is.na(kary.snp[,1]),]
kary.snp <- kary.snp[as.character(kary.snp[,11]) %in% as.character(pollen_genos[,1]),]

kary.3 <- kary.snp[as.character(kary.snp[,2]) == 'EG9_Chr3',]
kary.10a <- kary.snp[as.character(kary.snp[,2]) == 'EG9_Chr10a',]
kary.10b <- kary.snp[as.character(kary.snp[,2]) == 'EG9_Chr10b',]

gr.kary3 <- makeGRangesFromDataFrame(kary.3,keep.extra.columns=T)
gr.kary10a <- makeGRangesFromDataFrame(kary.10a,keep.extra.columns=T)
gr.kary10b <- makeGRangesFromDataFrame(kary.10b,keep.extra.columns=T)

autoplot(gr.kary3, layout="karyogram",aes(color=geno,fill=geno)) 
autoplot(gr.kary10a, layout="karyogram",aes(color=geno,fill=geno)) 
autoplot(gr.kary10b, layout="karyogram",aes(color=geno,fill=geno))



