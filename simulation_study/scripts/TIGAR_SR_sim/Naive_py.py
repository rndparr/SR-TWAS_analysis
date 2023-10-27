#/usr/bin/env python

###############################################################
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys

from time import time

import numpy as np
import pandas as pd

# getting zetas
from scipy.optimize import minimize

# estimator classes
from sklearn.base import BaseEstimator
from sklearn.ensemble import StackingRegressor
from sklearn.model_selection import cross_val_score

# R2, Pval
import statsmodels.api as sm

###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='SR-TWAS script.')
# Specify tool directory
parser.add_argument('--SR_TWAS_dir', type=str)
# other arguments
parser.add_argument('--chr', type=str, dest='chrm')
parser.add_argument('--cvR2', type=int, choices=[0, 1], default=1)
parser.add_argument('--cvR2_threshold', type=float, default=0.005)
parser.add_argument('--format', type=str, dest='data_format', choices=['GT', 'DS'], default='GT')
parser.add_argument('--gene_exp', type=str, dest='geneexp_path')
parser.add_argument('--genofile', type=str, dest='geno_path')
parser.add_argument('--genofile_type', type=str, choices=['vcf', 'dosage'])
parser.add_argument('--hwe', type=float, default=0.00001)
parser.add_argument('--maf_diff', type=float, default=0)
parser.add_argument('--missing_rate', type=float, default=0.2)
# parser.add_argument('--naive', type=int, choices=[0, 1], default=1)
parser.add_argument('--out_dir', type=str)
parser.add_argument('--log_file', type=str, default='')
# parser.add_argument('--out_prefix', type=str, default='')
# parser.add_argument('--out_info_file', type=str, default='')
# parser.add_argument('--out_weight_file', type=str, default='')
parser.add_argument('--out_naive_prefix', type=str, default='')
parser.add_argument('--out_naive_info_file', type=str, default='')
parser.add_argument('--out_naive_weight_file', type=str, default='')
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path')
parser.add_argument('--weight_threshold', type=float, default=0)
parser.add_argument('--weights', dest='w_paths', nargs='+', default=[], type=str)
parser.add_argument('--window', type=int, default=1000000)
args = parser.parse_args()
sys.path.append(args.SR_TWAS_dir)

###############################################################
# set output file names
# if not args.out_prefix:
# 	args.out_prefix = 'CHR' + args.chrm + '_SR_train'

# if not args.out_info_file:
# 	args.out_info_file = args.out_prefix + '_GeneInfo.txt'

# if not args.out_weight_file:
# 	args.out_weight_file = args.out_prefix + '_eQTLweights.txt'


#############################################################
# Print input arguments to log
print(
'''********************************
Input Arguments:
Gene Annotation and Expression file: {geneexp_path}
Training sampleID file: {sampleid_path}
Chromosome: {chrm}
K (number of trained input models): {K}
cis-eQTL weight files:{w_paths_str}
Training genotype file: {geno_path}
Genotype file used for training is type: {genofile_type}
Genotype data format: {data_format}
Gene training region SNP inclusion window: +-{window}
Excluding SNPs if missing rate exceeds: {missing_rate}
{maf_diff_str1}xcluding SNPs matched between eQTL weight file and training genotype file {maf_diff_str2}
HWE p-value threshold for SNP inclusion: {hwe}
{cvR2_str1} Stacked Regression model by 5-fold cross validation{cvR2_str2}.
Number of threads: {thread}
Output directory: {out_dir}
Output training info file for SR-TWAS: {naive_info}
Output trained weights file for SR-TWAS: {naive_weight}
********************************'''.format(
	**args.__dict__,
	cvR2_str1 = {0:'Skipping evaluation of', 1:'Evaluating'}[args.cvR2],
	cvR2_str2 = {0:'', 1:' with inclusion threshold Avg.CVR2 >' + str(args.cvR2_threshold)}[args.cvR2],
	maf_diff_str1 = {0:'Not e', 1:'E'}[do_maf_diff],
	maf_diff_str2 = {0:'by MAF difference.', 1:'if MAF difference exceeds: |' + str(args.maf_diff) + '|'}[do_maf_diff],
	w_paths_str = '\n  '.join(args.w_paths),
	K = Kweights,
	naive_info = out_naive_info_path,
	naive_weight = out_naive_weight_path))

# tg.print_args(args.__dict__)

##########################
# STUFF FOR TEST FILES
##########################
# Gene Expression header, sampleIDs
print('Reading genotype, expression file headers, sample IDs.\n')
sampleID, sample_size, geno_info, exp_info = tg.sampleid_startup(**args.__dict__)

# Read in gene expression info
print('Reading gene expression data.\n')
GeneExp, TargetID, n_targets = tg.read_gene_annot_exp(**exp_info)

# get info for the weight files
weights_info = tg.weight_k_files_info(**args.__dict__)

ES_cols = [ES(k) for k in range(Kweights)]

# set up output files

### SR-TWAS OUTPUT
out_weight_cols = ['CHROM','POS','ID','REF','ALT','TargetID','MAF','p_HWE','ES']

# dictionary of dtypes for info output
info_wk_dtypes = merge_dicts([{W(k)+'_N_SNP':np.int64, W(k)+'_CVR2':np.float64, W(k)+'_R2':np.float64, W(k)+'_PVAL':np.float64} for k in range(Kweights)])
info_wk_cols = list(info_wk_dtypes.keys())

info_zeta_dtypes = {'Z'+str(k): np.float64 for k in range(Kweights)}
info_zeta_cols = list(info_zeta_dtypes.keys())

info_dtypes = {**info_zeta_dtypes, **info_wk_dtypes}

out_info_cols = ['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName','sample_size','N_SNP','N_EFFECT_SNP','CVR2', 'R2', 'PVAL'] + info_zeta_cols + info_wk_cols

### NAIVE OUTPUT
print('Creating file: ' + out_naive_weight_path + '\n')
pd.DataFrame(columns=out_weight_cols).to_csv(
	out_naive_weight_path,
	header=True,
	index=None,
	sep='\t',
	mode='w')

out_naive_info_cols = ['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName','sample_size','N_SNP','N_EFFECT_SNP','CVR2', 'R2', 'PVAL'] + info_wk_cols
info_naive_dtypes = info_wk_dtypes
print('Creating file: ' + out_naive_info_path + '\n')
pd.DataFrame(columns=out_naive_info_cols).to_csv(
	out_naive_info_path,
	sep='\t',
	header=True,
	index=None,
	mode='w')

print('********************************\n')

def get_Geno():
	target = TargetID[0]
	Expr = GeneExp.iloc[[0]]

	start = str(max(int(Expr.GeneStart) - args.window, 0))
	end = str(int(Expr.GeneEnd) + args.window)

	# Query files; function stops here if no data in genotype and/or no weight file data
	tabix_target_ks, empty_target_ks = tg.tabix_query_files(start, end, **args.__dict__)

	# READ IN AND PROCESS GENOTYPE DATA 
	# file must be bgzipped and tabix
	Geno = tg.read_tabix(start, end, sampleID, **geno_info)

	# filter out variants that exceed missing rate threshold
	Geno = tg.handle_missing(Geno, sampleID, args.missing_rate)

	# maf
	Geno = tg.calc_maf(Geno, sampleID, 0, op=operator.ge)

	# get, filter p_HWE
	Geno = tg.calc_p_hwe(Geno, sampleID, args.hwe)

	# center data
	Geno = tg.center(Geno, sampleID)

	return start, end, tabix_target_ks, empty_target_ks, Geno


##############################################################
# thread function
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Expr = GeneExp.iloc[[num]]

	# center data
	Expr = tg.center(Expr, sampleID)

	# # Query files; function stops here if no data in genotype and/or no weight file data
	tabix_target_ks = tabix_target_ks_orig
	empty_target_ks = empty_target_ks_orig

	# # READ IN AND PROCESS GENOTYPE DATA 
	Geno = Geno_orig.copy()

	# read in weight data for target from all K weight files
	Weights = pd.DataFrame()
	target_ks = []

	for k in tabix_target_ks:
		# read in from tabix, may be empty
		W_k = tg.read_tabix(start, end, target=target, raise_error=False, **weights_info[k])

		while not W_k.empty:
			## CURRENTLY IGNORES ENSEMBLE VERSION
			# ## CONSIDER CHECKING IF ('.' in target) and requiring ensemble version if specified in the gene annotation file; otherwise use whatever ensemble version in the weight file
			# if np.any(W_k.TargetID.str.contains('.')):
			# 	W_k['TargetID'] = W_k.TargetID.str.split('.', expand=True)[0]

			W_k['snpIDflip'] = tg.get_snpIDs(W_k, flip=True)

			W_k = W_k[np.any(W_k[['snpID','snpIDflip']].isin(Geno.snpID.values), axis=1)].reset_index(drop=True)

			# if not in Geno.snpIDs, assumed flipped; if flipped, 1 - MAF and -ES
			flip = np.where(W_k.snpID.isin(Geno.snpID), True, False)

			if not np.all(flip):
				# set correct snpID, MAF, ES
				W_k['snpID'] = np.where(flip, W_k.snpID, W_k.snpIDflip)
				if args.maf_diff:
					W_k[MAF(k)] = np.where(flip, W_k[MAF(k)], 1 - W_k[MAF(k)])
				W_k[ES(k)] = np.where(flip, W_k[ES(k)], -W_k[ES(k)])

			W_k = W_k.drop(columns=['CHROM','POS','REF','ALT','TargetID','snpIDflip'])

			target_ks.append(k)

			if Weights.empty:
				Weights = W_k
			else:
				Weights = Weights.merge(W_k, how='outer', on=['snpID'])

			break

		else:
			empty_target_ks.append(k)

	if Weights.empty:
		print('No cis-eQTL weights with snp overlap in genotype data for target.\n')
		return None

	# ks may be added in read in step, need to sort for correct output later
	empty_target_ks.sort()

	# add nan columns for weight files without data
	if args.maf_diff:
		empty_wk_cols = flatten([[ES(k),MAF(k)] for k in empty_target_ks])
	else: 
		empty_wk_cols = [ES(k) for k in empty_target_ks]
	Weights[empty_wk_cols] = np.nan

	# get overlapping snps, filter Geno
	Geno = Geno[Geno.snpID.isin(Weights.snpID.values)].reset_index(drop=True)

	# merge datasets (inner in case missing snpIDs after filtering)
	Train = Geno.merge(Weights,
		left_on='snpID', 
		right_on='snpID', 
		how='inner').set_index('snpID')

	# filter by MAF diff - if an MAF# differs too much, set ESk to nan
	if args.maf_diff:
		for k in range(Kweights):
			maf_diff = np.abs(Train[MAF(k)].values - Train['MAF'].values)
			Train.loc[maf_diff > args.maf_diff, ES(k)] = np.nan

	# filter out snps where all weights are nan
	all_missing_weights = Train[ES_cols].count(axis=1) == 0
	Train = Train[~all_missing_weights]
	n_snps = Train.index.size

	if not n_snps:
		print('No valid SNPs.')
		# print('All SNP MAFs for training data and testing data differ by a magnitude greater than ' + str(args.maf_diff) + '.\n')
		return None

	# do stacked regression
	print('Performing Stacked Regression.')
	X = Train[sampleID].T
	Y = Expr[sampleID].values.ravel()

	# SET UP MODELS
	if len(target_ks) == 1:
		# only one weight model with usable data, don't need to do stacking
		print('Only trained model ' + str(target_ks[0]) + ' has usable weights for target.')
		reg = WeightEstimator(Train[ES(target_ks[0])], W(target_ks[0])).fit(X, Y)
		# Train['Naive'] = Train[ES(target_ks[0])]
		# naive_reg = reg
	else:
		weight_estimators = [(str(k), WeightEstimator(Train[ES(k)], W(k))) for k in target_ks]
		reg = WeightStackingRegressor(estimators=weight_estimators).fit(X, Y)
		# Train['Naive'] = np.nanmean(Train[ES_cols].values, axis=1)
		# Train['Naive'] = np.mean(np.nan_to_num(Train[ES_cols].values), axis=1)
		# naive_reg = WeightEstimator(Train['Naive'], 'Naive').fit(X, Y)
	# treat missing weight data for base model as all-zero
	Train['Naive'] = np.mean(np.nan_to_num(Train[ES_cols].values), axis=1)
	naive_reg = WeightEstimator(Train['Naive'], 'Naive').fit(X, Y)

	# 5-FOLD CV
	if args.cvR2:
		naive_avg_r2_cv = naive_reg.avg_r2_cv(X, Y)

		print('Average 5-fold CV R2 for naive: {:.4f}'.format(naive_avg_r2_cv))

		if naive_avg_r2_cv < args.cvR2_threshold:
			print('Average R2 < ' + str(args.cvR2_threshold) + '; Skipping naive training for target.\n')

	else: 		
		naive_avg_r2_cv = 0

	# needed for naive method
	wk_out_vals = final_wk_out_vals(reg.out_vals(X, Y), Kweights)

	# OUTPUT
	if (naive_avg_r2_cv >= args.cvR2_threshold) or (not args.cvR2):	
		naive_avg_r2_cv = naive_reg.avg_r2_cv(X, Y)
		naive_r2, naive_pval = naive_reg.r2_pval(X, Y)

		# output trained weights from naive method
		Naive_Weight = Train[['CHROM','POS','REF','ALT','MAF','p_HWE']].copy()
		Naive_Weight['ID'] = Naive_Weight.index
		Naive_Weight['TargetID'] = target
		Naive_Weight['ES'] = Train['Naive']
		# filter for non-zero effect size
		Naive_Weight = Naive_Weight[Naive_Weight['ES'] != 0]
		# reorder columns
		Naive_Weight = Naive_Weight[out_weight_cols]
		Naive_Weight.to_csv(
			out_naive_weight_path,
			sep='\t',
			header=None,
			index=None,
			mode='a')

		# output training info
		Naive_Info = Expr[['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()
		Naive_Info['sample_size'] = Y.size
		Naive_Info['N_SNP'] = n_snps
		Naive_Info['N_EFFECT_SNP'] = Naive_Weight.ES.size
		Naive_Info['CVR2'] = naive_avg_r2_cv
		Naive_Info['R2'] = naive_r2
		Naive_Info['PVAL'] = naive_pval
		Naive_Info = Naive_Info.join(pd.DataFrame.from_records({**wk_out_vals}, index=[0]).astype(info_naive_dtypes))[out_naive_info_cols]
		Naive_Info.to_csv(
			out_naive_info_path,
			sep='\t',
			header=None,
			index=None,
			mode='a')

	print('')

###############################################################
# thread process
if __name__ == '__main__':
	print('Starting training for ' + str(n_targets) + ' target genes.\n')
	start, end, tabix_target_ks_orig, empty_target_ks_orig, Geno_orig = get_Geno()
	pool = multiprocessing.Pool(args.thread)
	pool.imap(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')

	# sort tabix output
	sort_tabix_output(out_naive_weight_path, args.out_dir + '/' + args.out_naive_weight_file)


###############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

###############################################################

sys.stdout.close()
