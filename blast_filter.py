# !/usr/bin/env python

#This script filters the best blast hits from BLASTP alignment output: mapped protein sequences have >90% sequence aligned. Then identify gene functions from Arabidopsis functional annotations.

import pandas as pd
import numpy as np
import sys, os

names = ["qseqid","qlen","qstart","qend","sseqid","slen","sstart","send","length","pident","bitscore","evalue"]

blast = pd.read_table("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/BLAST/BLASTP/EG9_cds.blastpout6",sep=None,delimiter=None,names=names)
#blast = pd.read_table("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/seg_distortion/BLAST/BLASTP/cds_Oleifera.blastout6",sep=None,delimiter=None,names=names)
#blast = pd.read_table("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/seg_distortion/BLAST/BLASTP/compar_EG8_EG9/EG9_3_AT.blastout6",sep=None,delimiter=None,names=names)

#reformat qseqid and add percent coverage of alignment for query and ref seq to blast results
qseqid = np.array(list(blast['qseqid'].str.split('|')))
print qseqid


print len(qseqid), type(qseqid[0]), len(qseqid[0])
blast_genes = [s.replace('ID=','') for s in list(qseqid[:,0])]
blast['qseqid'] = blast_genes
blast.insert(loc=1, column='position', value=qseqid[:,1])
blast['length'] =blast['length'].astype(float)
blast['prop_q'] = blast['length']/blast['qlen']
blast['prop_s'] = blast['length']/blast['slen']

print blast

#Filtering criteria
#1. length of prop_q and prop_s > 0.90
#2. sort each blast query block by pident,prop_q,prop_s,bitscore per group of alignments 

blast = blast[(blast['prop_q'] > 0.90) & (blast['prop_s'] > 0.90)]
uniq_genes = blast['qseqid'].unique()

#print blast['qseqid']
#print uniq_genes

sig_hits = []
for gene in uniq_genes:
        gene_blk = blast[blast['qseqid']==gene]
        gene_blk = gene_blk.sort(['pident','prop_q','prop_s','bitscore'],ascending=False)
        print gene
        sig_hits.append(list(gene_blk.iloc[0]))

print sig_hits

#reformat ref gene id:
ref_gene = [s.split('.')[0] for s in list(np.asarray(sig_hits)[:,5])]
sig_hits = np.hstack((np.asarray(sig_hits),np.reshape(np.asarray(ref_gene),(len(ref_gene),1))))

print sig_hits


#get annotations from Arabidopsis/Rice

gff = pd.DataFrame(np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/BLAST/BLASTP/A.thaliana_db/Arabidopsis_thaliana.TAIR10.38.gff3",dtype='str',delimiter= None,comments='#'))


gene_gff = gff[gff[2]=="gene"]

gene_anno = []
for hit in sig_hits:
	print hit
        anno = gene_gff[gene_gff[8].str.contains(hit[15])]
        gene_anno.append([list(anno.iloc[0])[8].split('[')[0]])
        

print np.asarray(gene_anno), len(gene_anno)


gene_anno = np.hstack((sig_hits[:,[0,1,2,5,6,9,10,11,12,13,14,15]],np.asarray(gene_anno)))
print gene_anno


s = "%s\t%s\t"*4 + "%s\t%s\t%s\t%s\t%s"

header = 'qseqid\tposition\tqlen\tsseqid\tslen\tlength\tpident\tbitscore\tevalue\tprop_q\tprop_s\tseqid_gene\tannotation'

np.savetxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/BLAST/BLASTP/anno_EG9_AllChrom.txt",gene_anno,fmt=s,header=header)







