# !/usr/bin/env python

'''
Identify all of the CDS that correspond to genes of interest (i.e. genes in 
segregation distortion regions, protein-coding genes).  

#Usage: python create_cds_fasta.py <output_name> <seq> <prefix> <gff> <genelist>
Usage: python create_cds_fasta.py <output_name> <seq> <gff> <genelist> <splice_delimiter>
- <output_name> The name you wish to give to the output file.
- <seq> The reference genome for the species in FASTA format. The sequence for 
  each chromosome must be in one line right beneath the header.
#- <prefix> String or regular expression that denotes the prefix used with the 
#  chromsome number in the <seq> file. Do not include the ">" that is at the 
#  beginning of every header line in a FASTA file. 
- <gff> A gff file that only contains CDS sequences. Must be sorted by 
  chromosome, then start position, then end position. The chromosome number 
  should not include the <prefix>.
- <genelist> File of genes of interest with the following columns:
  (Sample file)
  # chr type    start   end     gene    
  1     gene    82274   87582   ID=p8_sc0666_maker_00151691;Name=p8_sc0666_maker_00151691;Alias=maker-p8_sc0666-augustus-gene-1.18
  This file must also be sorted by chromosome, then start position, and 
  then end position. The chromosome number should not include the <prefix>.
- <splice_delimiter> The prefix common to all of the transcripts wihtin gffs 
  that denotes how the alternative splices are enumerated. Include the delimiter 
  used to separate the different parts of the name of the transcript.
'''

import sys, os
import numpy as np
import scipy as sp
import pandas as pd
from string import maketrans
import re

def group_by_heading(some_source):
    ''' Reading FASTA file line by line. In FASTA file, each chr is 
    already on 1 line.
    '''
    buffer = []
    for line in some_source:
        if re.match('>', line):
            if buffer: yield buffer
            buffer= [ line ]
        else:
            buffer.append(line)
    yield buffer


def main(args):
    #gff = np.loadtxt('/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/seg_distortion/BLAST/EG8_cds_sort.gff3',dtype='str',delimiter=None, skiprows=1)
    gff = np.loadtxt(os.path.expanduser(args[3]), dtype='str', delimiter=None, skiprows=1)
    #seg_genes = np.loadtxt('/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/seg_distortion/BLAST/GBS_segdist_genes_gff.txt',dtype='str',delimiter=None)
    seg_genes = np.loadtxt(os.path.expanduser(args[4]), dtype='str', delimiter=None)
    print len(seg_genes)
    #prefix = re.compile(args[3])

    cds=[]
    for line in gff:
        ID = line[8].split(':')[0]
        #ID = ID.split("_R")[0] #_R denotes CDS
        ID = re.split(args[5], ID)[0]
        cds_id = line[8].split(":")[0]
        cds.append([line[0],line[2],line[3],line[4],line[6],cds_id,ID])
    #print cds

    sd_gene = []
    # Retreive gene id from annotation field for genes in seg dist regions
    for gene in seg_genes:
        ID = gene[4].split(";")[0]
        sd_gene.append([gene[0],gene[1],gene[2],gene[3],ID])
    #print sd_gene

    # Only care about the CDS regions that overlap with the genes of interest
    sd_cds_idx = np.nonzero(np.in1d(np.asarray(cds)[:,6],np.asarray(sd_gene)[:,4]))
    #print sd_cds_idx
    seg_cds = np.asarray(cds)[list(sd_cds_idx)]
    cds_df = pd.DataFrame(seg_cds[:,0:6], columns=['chr','type','start','end','dir','cds'])
    print cds_df
    sd_gene = pd.DataFrame(np.asarray(sd_gene),columns=['chr','type','start','end','gene'])
    print sd_gene
    headings = []
    FASTA=[]


    with open(args[2], "r") as source:
        for heading_and_lines in group_by_heading(source):
            #chrom = re.sub(prefix, "", heading_and_lines[0]).replace(">","").strip("\n")
            chrom = heading_and_lines[0].replace(">","").strip("\n")
            print chrom
            gene_coord = cds_df[cds_df['chr'] == str(chrom)]
            #print gene_coord
            seg_gene = sd_gene[sd_gene['chr'] == str(chrom)]
            #print seg_gene
            chrom_seq = heading_and_lines[1].strip()

            for gene in np.asarray(seg_gene):
                #print gene
                gene_block = gene_coord[gene_coord['cds'].str.contains(gene[4])]
                #print gene_block
                GENE_START = gene[2]
                GENE_END = gene[3]

                cds_uniq = list(gene_block['cds'].unique())
                CDS_NUM = len(list(cds_uniq)) #number of differenct CDS within a gene
                #print cds_uniq, CDS_NUM

                for cds_type in cds_uniq: #iterate through the cds groups/types (denotes alternative splice variants for a gene)
                    #print cds_type
                    cds_group=gene_block[gene_block['cds'].str.contains(cds_type)]
                    #print cds_group
                    SEQ = []
                    plus = "TAGC"
                    complement = "ATCG"
                    translate_minus = maketrans(plus,complement) # (Reverse) complementation
                    CDS_CAT = len(np.asarray(cds_group)) #Number of CDS that were concatenated
                    #print CDS_CAT

                    for c in np.asarray(cds_group):
                        #print c
                        SEQ.append(chrom_seq[int(c[2])-1:int(c[3])])

                    if c[4] == '+':
                        cds_cat =  "".join(SEQ)
                        #print "plus:", cds_cat
                        FASTA.append(['>'+'%s' %cds_type+'|'+heading_and_lines[0].strip('\n').replace('>','')+'_'+'%s' %GENE_START+':'+'%s' %GENE_END+'|'+'CDS_NUM:'+'%s' %CDS_NUM+'|'+"CDS_CAT:"+'%s' %CDS_CAT+'|'+'CDS_LEN:'+'%s' %len(cds_cat)])
                        FASTA.append([cds_cat])
                    if c[4] == '-':
                        cds_cat = "".join(SEQ).translate(translate_minus)[::-1] #concatenate,complement, then reverse seq       
                        #print "minus:", cds_cat
                        FASTA.append(['>'+'%s' %cds_type+'|'+heading_and_lines[0].strip('\n').replace('>','')+'_'+'%s' %GENE_START+':'+'%s' %GENE_END+'|'+'CDS_NUM:'+'%s' %CDS_NUM+'|'+"CDS_CAT:"+'%s' %CDS_CAT+'|'+'CDS_LEN:'+'%s' %len(cds_cat)])
                        FASTA.append([cds_cat])


    np.savetxt(args[1], np.asarray(FASTA), fmt="%s")


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Missing command line arguments to script. Please read documentation:')
        print(__doc__)
        sys.exit(1)
    else:
        main(sys.argv)


