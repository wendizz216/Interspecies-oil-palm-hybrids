# Oil-Palm-Hybrids
Using genome seqeuncing technology to improve inter-species hybrid breeding

### PROJECT MOTIVATION/OBJECTIVE
Oil palm cultivars are highly susceptible to a number of diseases. The potential social and economic risk involved to the multibillion-dollar industry has prompted the Malaysian oil palm government to start breeding disease resistant oil palm cultivars. Breeders have identified disease resistant wild oil palm endemic to Latin America--a sister species of the African cultivated varieties (diverged 50 Mya).] To generate palms that express disease resistance as well as other agronomically valuable wild traits, the breeders cross-bred wild Latin American with African cultivated palms to produce several generations of hybrids. However, despite three generations of selective breeding, the hybrid palms did not express the desirable wild traits of interest. Here we use genomic sequencing analysis to identify the genetic source preventing the transfer of wild traits into hybrid populations.

### PARTNER
[Malaysian Oil Palm Board](http://www.mpob.gov.my/)

### DATA SETS
* Whole genome sequencing data (N=228)
* Gene sequences
  * Download files from [here](https://www.dropbox.com/sh/hmzssojdk0m4qi7/AACzlYJK3cD1m6K8IeOSzsJ2a?dl=0)


### METHODS
* Chi-squared test (custom Python script)
* Functional characterization of genes
  * Protein coding gene sequence extraction from genome (extract_gene_seq.py)
  * Translate gene sequences to protein sequences
    - Download EMBOSS command-line tool [here](http://emboss.sourceforge.net/download/)
    - command-line function: <input> transeq eg9_cds_genes.fasta eg9_proteins.pep -trim Y
  * Genomic sequence alignments ([BLAST sequence alignment tool](https://blast.ncbi.nlm.nih.gov/Blast.cgi))
    - To run BLAST, you can download command-line tool [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    - Filter BLAST results to extract reliable sequence results (custom Python script)
* Data visualization (R, ggplot)

### PROJECT DESCRIPTION
We hypothesize that genes associated with traits that under-represented in the hybrid populations are linked to unfavorable combinations of genes that result in hybrid inviability. We can identify unfavorable combinations of genes by comparing gene frequencies in the population of hybrids to that of the expected population frequency inferred from the parents. 

* I used a simple chi-square test to scan the 1.8Gb genome for regions that deviate from expected gene ratios (estimated from the parental genetic sequences) and see if these regions contain genes associated with wild traits of interest. 

* Because the gene functions were not well annotated in the oil palm genome, I wasn't able to infer potential genetic mechanisms. To get around this issue, I mapped oil palm protein sequences onto a functionally well annotated plant genome (e.g., Arabidopsis) and inferred gene functions based on sequence similarities. I used a customized python script to extract gene sequences that were then translated into protein sequences and implemented BLAST sequence alignment tools to infer oil palm gene functions. 

### RESULTS
* Identified several large genomic regions of unfavorable gene combinations in hybrid populations and found that genes significantly associated with wild traits of interests are located within these regions.

* Regions of unfavorable gene combinations harbor genes responsible for pollen viability/sterility and other reproductive traits, suggesting that genes associated with wild traits of interest may be linked to genes that cause reproductive abnormalities in hybrids thus are not represented in the hybrid populations.

* Results were experimentally validated and provided new innovative insights towards developing molecular gene editing techniques to improve inter-species hybrid breeding. 
