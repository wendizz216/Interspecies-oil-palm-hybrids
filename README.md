# Oil-Palm-Hybrids
Improving inter-species hybrid breeding

### PROJECT MOTIVATION/OBJECTIVE
Oil palm cultivars are highly susceptible to a number of diseases. The potential social and economic risk involved to the multibillion-dollar industry has prompted the Malaysian oil palm industry to start breeding disease resistant oil palm cultivars. Breeders have identified disease resistant wild oil palm originating North/South America, which are a sister species of the African cultivated varieties (diverged 50 Mya). To generate palms that express disease resistant traits as well as other agronomically valuable traits, the breeders cross-bred wild disease resistant oil palm with cultivated oil palm to produce several generations of hybrids. Unfortunately, the hybrid palms did not express the desirable wild traits of interest. Here we use genomic sequence analysis to identify the genetic source preventing the transfer of wild traits in oil palm hybrid populations.

### PARTNER
[Malaysian Oil Palm Board](http://www.mpob.gov.my/)

### DATA SETS
* Whole genome sequencing data (N=228)
  - 1st, 2nd, 3rd generation of hybrid populations
* Genome reference build

### METHODS
* Chi-squared test (custom Python script)
* Functional characterization of genes
  * Protein coding gene sequence extraction from genome (custom Python script)
  * Translate gene sequences to protein sequences
    - Download EMBOSS command-line tool [here](http://emboss.sourceforge.net/download/)
  * Genomic sequence alignments ([BLAST sequence alignment tool](https://blast.ncbi.nlm.nih.gov/Blast.cgi))
    - To run BLAST, you can download command-line tool [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    - Filter BLAST results to extract reliable sequence results (custom Python script)
* Data visualization (R, ggplot)

### PROJECT DESCRIPTION
The goal of this project was to identify the genetic source preventing the effective transfer of wild oil palm traits into hybrid populations using traditional breeding methods. Genetic theory suggests that inter-species hybrids often suffer from inviability and that traits linked to inviability will not be expressed in the hybrid populations. Here I assume that unfavorable gene combinations in the hybrid populations will be under-represented due to the negative effects on hybrid viability. We hypothesize that, within the hybrid population of the parental species crosses, the frequency of genes associated with wild traits linked to inviability genes will fall outside of the expected Mendelian segregation ratios ([Mendelian Law of Inheritance](https://en.wikipedia.org/wiki/Mendelian_inheritance)). 

* I used a simple chi-square test to scan the 1.8Gb genome for regions that deviate from expected Mendelian segregation ratios (those that have unfavorable gene combinations) and see if these regions contain genes associated with wild traits of interest. 

* To infer genetic mechanisms related to inviability, I needed to understand functions of the genes. I used a customized python script to extract gene sequences that were then translated into protein sequences to identify gene function. Because the gene functions are not well annotated in the oil palm genome, I implemented BLAST sequence alignment tool to infer gene function based on sequence similarities of a functionally well annotated plant genome (e.g., Arabidopsis). 

### Results
* Identified several large genomic regions of unfavorable gene combinations in hybrid populations (i.e., gene frequencies that deviate from mendalian segregation ratios) and found that genes significantly associated with wild traits of interests are located within these regions.
* Regions of unfavorable gene combinations harbor genes responsible for pollen viability and other reproductive traits, suggesting that genes associated with wild traits of interest may be linked to genes that cause reproductive abnormalities in hybrids thus preventing the expression of these traits.
* Results were experimentally validated and provided new innovative insights towards developing molecular gene editing techniques to improve inter-species hybrid breeding. 
