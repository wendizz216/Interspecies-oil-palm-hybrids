# !/usr/bin/env python


'''
Identify genes/genotypes with hybrid population frequencies deviating from the expected gene/genotype frequencies estimated from the parental genes/genotypes.

#Usage: python chisquared_test.py <input_file> <pop> <sample_ID> <parent_ID> <output_name>
- <input_file> The name of a tab file of gene sequence variations of all samples (palms). First 3 cols correspond to chrom/pos/ref-allele; remaining columns correspond to seqeunced samples and rows, sequence variation across positions in the genome.
- <pop> String that denotes the hybrid population to test. Enter one of the population ID options: OxG, BC1a, BC1b, BC2a, BC2b.
- <sample_ID> File with all column names of the gene sequence variation tab file
- <parent_ID> File with parental samples for each hybrid population. 
- <output_name> The name you wish to give to the output file.
'''

import os, sys, os.path
import numpy as np
import scipy as sp
import itertools
import collections
from collections import Counter
import pandas as pd
import scipy.stats as stats
import re

def main(args):

	#tab = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/github/test.tab",dtype='str',delimiter='\t')
	#samples = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/QC/sample_names.txt",dtype=str,delimiter='\t')
	#par = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/github/parental_samples.csv",delimiter=",",dtype=str)

	tab = np.loadtxt(os.path.expanduser(args[1]),dtype='str',delimiter='\t')
	samples = np.loadtxt(os.path.expanduser(args[3]),dtype='str',delimiter='\t')
	par = np.loadtxt(os.path.expanduser(args[4]),dtype='str',delimiter=',')
	pop_id = args[2]
	print pop_id
	
	###FIRST: PARENT SNP POS --> Extract parental genotype sequences from hybrid population of interest

	def parental_seq(pop_type):
		par_samples = par[par[:,0] == pop_type,1]
		par_idx = [i for i, item in enumerate(list(samples)) if item in list(par_samples)]
		par_seq = tab[:,par_idx]
		return par_seq
	
	
	PAR_SEQ = parental_seq(pop_id)

	###SECOND: HYBRID SNP POS --> Extract hybrid genotypes from hybrid population of interest

	#Remove parental genotypes from the columns of the tab file
	par_idx = [i for i, item in enumerate(list(samples)) if item in list(par[:,1])]
	tab_nopar = np.delete(tab,par_idx, axis=1)
	sample_nopar = np.delete(samples, par_idx, None)

	def hybrid_seq(pop_type):
		hyb_idx = [i for i, item in enumerate(sample_nopar) if pop_type in item]	
		hyb_seq = tab_nopar[:,hyb_idx]
		return hyb_seq

	HYB_SEQ = hybrid_seq(pop_id) 


	###THIRD: Filter/Skip SNPs in which both parents are homozygous and expected to only produce heterozygous hybrid offspring

	snp_skip_idx = [] #keeping track of snps where both parents are homozygous and is expected produce all hetero offspring
	SNP_TEST_IDX = [] #keeping track of snps that are used to test segregation distortion

	for idx, snp in enumerate(PAR_SEQ):
		#skip snp sites in which both parents are homozygous or a missing allele in genotype call
		if snp[0][0] == "." or snp[0][2]=="." or snp[1][0] == "." or snp[1][2] == ".":
			#print snp, "missing allele"
			snp_skip_idx.append(idx)
			continue

		if snp[0][0] == snp[0][2] and snp[1][0] == snp[1][2]:
			#print snp, "produces all homos"
			snp_skip_idx.append(idx)
			continue

		else:
			#print snp, "test snps"
			SNP_TEST_IDX.append(idx)

	print SNP_TEST_IDX

	###FOURTH: EXPECTED MENDELIAN GENOTYPES --> preparing unique hybrid genotype calls

	#Clean genotype calls: remove missing data
	HYB_GENO = [] #all genotype types that are found in the hybrid pop
	for idx, snp in enumerate(HYB_SEQ):
		snp1=list(snp)
		cleaned = [x for x in snp1 if "./." not in  x] #remove missing genotypes "./."
		#print list(np.unique(cleaned)), len(snp1)
		HYB_GENO.append(list(np.unique(cleaned)))

	print HYB_GENO

	EXP_HYB_GENO = [] #use this to rm unexpected genotype calls in hybrid pop based on parental genotype calls
	for idx, snp in enumerate(PAR_SEQ):
		#generate exp genotype counts for hybrid populaiton only for heterozygous sites.
		if idx in SNP_TEST_IDX:
			offspring = HYB_GENO[idx]

			exp_geno=[]

			geno1a = snp[0][0]+ "/" + snp[1][0]
			geno1b = snp[1][0]+ "/" + snp[0][0]

			geno2a = snp[0][0]+ "/" +snp[1][2]
			geno2b = snp[1][2]+ "/" +snp[0][0]

			geno3a = snp[0][2]+ "/" +snp[1][0]
			geno3b = snp[1][0]+ "/" +snp[0][2]

			geno4a = snp[0][2]+ "/" +snp[1][2]
			geno4b = snp[1][2]+ "/" +snp[0][2]
			#Must identify the correct order of heterozygote genotypes by using the genotypes called in offspring pop               
			if geno1a not in offspring:
				exp_geno.append(geno1b)
			else:
				exp_geno.append(geno1a)

			if geno2a not in offspring:
				exp_geno.append(geno2b)
			else:
				exp_geno.append(geno2a)

			if geno3a not in offspring:
				exp_geno.append(geno3b)
			else:
				exp_geno.append(geno3a)

			if geno4a not in offspring:
				exp_geno.append(geno4b)
			else:
				exp_geno.append(geno4a)

			EXP_HYB_GENO.append(exp_geno)

	print EXP_HYB_GENO

	###FIFTH: PREPARING OBSERVED GENOTYPE CALLS

	OBS_GENO_CNTS = [] #excludes unexpected genotypes (based on parental genotype)
	obs_geno_cnts = [] #includes unexpected genotypes (based on parental genotype)
	obs_cnts = [] #includes unexpected genotypes

	for idx, snp in enumerate(HYB_SEQ):
		snp1= [x for x in snp if "./." not in x] #check list for missing './.' genotype calls and then remove
		#print snp1             

		if idx in SNP_TEST_IDX: #assess genotype frequency of snps in which both parents are NOT homozygous
			#print idx
			#print idx, snp1, len(snp1)
			obs_geno_cnts.append(snp1)

	#removing unexpected genotypes: only include hyb genotypes that are expected from parental genotype calls
	for snp, geno in zip(obs_geno_cnts, EXP_HYB_GENO):
		trim_geno = [item for item in snp if item in geno]
		OBS_GENO_CNTS.append(trim_geno)
		obs_cnts.append(dict(Counter(trim_geno)))	

	###SIXTH: PREPARE EXPECTED GENTOYPE COUNTS/PROPORTION --> calculated based on sample size of observed data

	np.random.seed(10) #changed seed and the results still hold for multiple runs 

	EXP_GENO_CNTS = [] 
	for geno,obs in zip(EXP_HYB_GENO,OBS_GENO_CNTS):
		#print list(geno), ['%s' %geno[0],'%s' %geno[1],'%s' %geno[2],'%s' %geno[3]],"%i" %len(obs), int("%i" %len(obs))
		if len(obs) <= 20: #SET GENOTYPE CALL RATE FILTER
			exp_prop = np.random.choice(a=['%s' %geno[0],'%s' %geno[1],'%s' %geno[2],'%s' %geno[3]], p= [0.25, 0.25, 0.25, 0.25], size = 20)
			EXP_genos.append(list(exp_prop))
		if len(obs) >20: #SET GENOTYPE CALL RATE FILTER
			exp_prop = np.random.choice(a=['%s' %geno[0],'%s' %geno[1],'%s' %geno[2],'%s' %geno[3]], p= [0.25, 0.25, 0.25, 0.25], size = int("%i" %len(obs)))
			#print exp_prop, len(list(exp_prop))
			EXP_GENO_CNTS.append(list(exp_prop))


	print EXP_GENO_CNTS
	print OBS_GENO_CNTS

	#########CHI-SQUARE TEST for deviations from expected Mendelian segregation ratios######################################

	N_DF_CritVal_ChiSq_Pval = []
	N=0

	for exp, obs in zip(EXP_GENO_CNTS,OBS_GENO_CNTS):
		if len(obs) < 20: #SET GENOTYPE CALL RATE FILTER
			N_DF_CritVal_ChiSq_Pval.append([len(obs),"NA", "NA", "NA", "NA"])
			continue
		if len(exp) == 0 or len(obs) == 0 or len(obs) == 1:
			N_DF_CritVal_ChiSq_Pval.append([len(obs), "NA","NA", "NA", "NA"])
			continue

		exp_cnt = dict(Counter(exp))
		obs_cnt = dict(Counter(obs))
		#print len(exp_cnt), obs_cnt

		geno_cat = exp_cnt.keys() #genotype categories
		geno_cat_obs = obs_cnt.keys()

		#If the number of genotype categories in expected and observed is not the same, ignore it
		if len(geno_cat) != len(geno_cat_obs):
			N_DF_CritVal_ChiSq_Pval.append([len(obs),"NA", "NA", "NA", "NA"])
			continue

		if len(geno_cat) == len(geno_cat_obs):
			cnt_exp = np.asarray(exp_cnt.values()) #count of expected genotypes
			cnt_obs = np.asarray(obs_cnt.values()) #count of observed genotypes

			#PREP EXPECTED CONTINGENCY TABLE
			expected = pd.DataFrame(np.reshape(cnt_exp,(len(cnt_exp),1)))
			expected.columns = ["Geno"]
			expected.index = geno_cat

			#PREP OBSERVED CONTINGENCY TABLE
			observed = pd.DataFrame(np.reshape(cnt_obs,(len(cnt_obs),1)))
			observed.columns = ["Geno"]
			observed.index = geno_cat
			#print expected, observed
			#print obs_cnt.values(), sum(obs_cnt.values())  

			#CHI-SQUARED TEST
			chi_squared_stat = (((observed-expected)**2)/expected).sum().sum()
			#print "chisq_stat:", chi_squared_stat

			if len(geno_cat) == 2: #for number of geno categories = 2 for df = 2-1
				n_size = sum(obs_cnt.values())
				df = len(geno_cat)-1
				crit = stats.chi2.ppf(q = 0.95, # Find the critical value for 95% confidence*
					      df = 1)   # *
				#print "crit_val:", crit

				p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,  # Find the p-value
						     df=1)
				#print "Pval:", p_value
				N_DF_CritVal_ChiSq_Pval.append(['%i' %n_size, '%i' %df, crit,chi_squared_stat, p_value])

			if len(geno_cat) == 3: #for number of geno categories = 3 for df = 3-1
				n_size = sum(obs_cnt.values())
				df = len(geno_cat)-1
				crit = stats.chi2.ppf(q = 0.95, # Find the critical value for 95% confidence*
					      df = 2)   # *
				#print "crit_val:", crit

				p_value = 1 - stats.chi2.cdf(x=chi_squared_stat,  # Find the p-value
						     df=2)
				#print "Pval:", p_value
				N_DF_CritVal_ChiSq_Pval.append(['%i' %n_size, '%i' %df, crit,chi_squared_stat, p_value])
		else:
			N_DF_CritVal_ChiSq_Pval.append([len(obs), "NA","NA", "NA", "NA"])


	print N_DF_CritVal_ChiSq_Pval
	TEST_SNP_IDX = np.asarray(SNP_TEST_IDX).reshape((len(SNP_TEST_IDX),1))

	###PREPARE OUTPUT FILE

	snp_pos_tested = tab[SNP_TEST_IDX,0:2]
	par_snp_test = PAR_SEQ[SNP_TEST_IDX,:]

	exp_off_geno_cnt = [] #expected genotype frequency
	for geno, obs in zip(EXP_GENO_CNTS,OBS_GENO_CNTS):
		#print Counter(geno), Counter(obs)
		exp_off_geno_cnt.append(dict(Counter(geno)))

	OBS_cnts = np.asarray(obs_cnts)

	print len(snp_pos_tested), len(exp_off_geno_cnt), len(obs_cnts)

	#Modify chromosome IDs for making figures in R
	ChromMod_test = []
	for snp in snp_pos_tested:
		chrom = snp[0].replace('EG9_Chr','')
		if chrom[-1] == 'a':
			chrom_mod = chrom.replace('a','')
			ChromMod_test.append([chrom_mod])
		if chrom[-1] == 'b':
			chrom_mod = chrom.replace('b','.5')
			ChromMod_test.append([chrom_mod])
		if chrom[-1] != 'a' and chrom[-1] != 'b':
			ChromMod_test.append([chrom])


	chisq_stats = np.hstack((TEST_SNP_IDX, ChromMod_test, snp_pos_tested, np.asarray(N_DF_CritVal_ChiSq_Pval), np.asarray(par_snp_test), np.reshape(exp_off_geno_cnt, (len(exp_off_geno_cnt),1)), np.reshape(OBS_cnts,(len(OBS_cnts),1))))

	print chisq_stats[0], len(chisq_stats[0])
	print chisq_stats

	header=np.reshape(np.array(("IDX","CHR_MOD","CHROM","POS","N","DF","Crit_Val","ChiSq_val","Pval", 'Parental1_genotype', 'Parent2_genotype', "EXP_Counts", "OBS_Counts" )),(1,13))

	fin = np.vstack((header,chisq_stats))
	s = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'

	np.savetxt(args[5],fin,fmt=s)


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Missing command line arguments to script. Please read documentation:')
        print(__doc__)
        sys.exit(1)
    else:
        main(sys.argv)

