#/usr/bin/env python

###############################################################
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys

import time

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
start_time = time.time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='SR-TWAS script.')

parser.add_argument('--annot', type=str, dest='annot_path')


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
parser.add_argument('--naive', type=int, choices=[0, 1], default=1)
parser.add_argument('--out_dir', type=str)
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--out_info_file', type=str, default='')
parser.add_argument('--out_weight_file', type=str, default='')
parser.add_argument('--out_naive_prefix', type=str, default='')
parser.add_argument('--out_naive_info_file', type=str, default='')
parser.add_argument('--out_naive_weight_file', type=str, default='')
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path')
parser.add_argument('--weight_threshold', type=float, default=0)
parser.add_argument('--weights', dest='w_paths', nargs='+', default=[], type=str)
parser.add_argument('--window', type=int, default=1000000)
parser.add_argument('--no_log', action='store_true')
args = parser.parse_args()
sys.path.append(args.SR_TWAS_dir)

###############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm + '_SR_train'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_info_file:
	args.out_info_file = args.out_prefix + '_GeneInfo.txt'

if not args.out_weight_file:
	args.out_weight_file = args.out_prefix + '_eQTLweights.txt'

if not args.out_naive_prefix:
	args.out_naive_prefix = 'CHR' + args.chrm + '_Naive_train'

if not args.out_naive_info_file:
	args.out_naive_info_file = args.out_naive_prefix + '_GeneInfo.txt'

if not args.out_naive_weight_file:
	args.out_naive_weight_file = args.out_naive_prefix + '_eQTLweights.txt'

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, 'SR_CHR' + args.chrm)
else:
	out_sub_dir = args.out_dir

# Make output, log directories
os.makedirs(out_sub_dir, exist_ok=True)
os.makedirs(os.path.join(args.out_dir, 'logs'), exist_ok=True)

# set stdout to log
if not args.no_log:
	sys.stdout = open(os.path.join(args.out_dir, 'logs', args.log_file), 'w')

# Check tabix command
try:
	subprocess.check_call(['which','tabix'], stdout=subprocess.DEVNULL)
except:
	raise SystemExit('Error: Required tool TABIX is not available.\n')


# Check gene expression file
try:
	os.path.isfile(args.geneexp_path)
except:
	SystemExit('Error: Gene expression file does not exist.')

# Check genotype file 
try:
	os.path.isfile(args.geno_path)
except:
	SystemExit('Error: Training genotype file does not exist.')

# Check training sample ID file
try:
	os.path.isfile(args.sampleid_path)
except:
	SystemExit('Error: Training sample ID file does not exist.')

###############################################################
import TIGARutils as tg


def W(k):
	return 'W'+str(k)

def ES(k):
	return 'ES'+str(k)

def MAF(k):
	return 'MAF'+str(k)

def flatten(nested_list):
	return [j for i in nested_list for j in i]

def merge_dicts(dict_list):
	return_dict = {}
	for d in dict_list:
		return_dict.update(d)
	return return_dict

# return R^2
def get_r2(y, predY, pval=False):
	lm = sm.OLS(y, sm.add_constant(predY)).fit()
	if pval:
		return lm.rsquared, lm.f_pvalue   
	return lm.rsquared

def _mse(y, ypred):
	mat = y - ypred
	return np.dot(mat, mat)

# format final output
def final_wk_out_vals(target_k_outvals, K):
	out_dict = merge_dicts([{W(k)+'_N_SNP': 0, W(k)+'_CVR2': 0, W(k)+'_R2': 0, W(k)+'_PVAL': 1} for k in range(K)])
	# if naive:
	# 	out_dict.update({'Naive_N_SNP': 0, 'Naive_CVR2': 0, 'Naive_R2': 0, 'Naive_PVAL': 1})
	out_dict.update(target_k_outvals)
	return out_dict

def final_zetas(target_k_zetas, K):
	out_dict = {'Z'+str(k):0 for k in range(K)}
	out_dict.update(target_k_zetas)
	return out_dict

# estimator for individual training models
class WeightEstimator(BaseEstimator):
	_estimator_type = 'regressor'

	def __init__(self, raw_weights, name):
		self.raw_weights = raw_weights
		self.name = name

	def fit(self, X=None, y=None):
		self.coef_ = self.raw_weights.dropna()
		self.snpID = self.coef_.index.values
		self.n_features_in_  = self.coef_.size
		return self

	def predict(self, X):
		return np.dot(X[self.snpID], self.coef_)

	def score(self, X, y):
		return get_r2(y, self.predict(X))

	def r2_pval(self, X, y):
		return get_r2(y, self.predict(X), pval=True)

	def avg_r2_cv(self, X, y):
		return sum(cross_val_score(self, X, y)) / 5

	def out_vals(self, X, y):
		r2, pval = self.r2_pval(X, y)
		return {self.name+'_N_SNP': self.n_features_in_, 
			self.name+'_CVR2': self.avg_r2_cv(X, y), 
			self.name+'_R2': r2, 
			self.name+'_PVAL': pval}

# final estimator for stacking regressor
class ZetasEstimator(BaseEstimator):
	_estimator_type = 'regressor'

	def __init__(self, min_method=None, tol=None):
		super().__init__()
		self.min_method = min_method
		self.tol = tol

	def _loss_func(self, zeta, X, y):
		if (len(zeta) == 1):
			Zeta = np.array([*zeta, 1 - np.sum(zeta)])
		else:
			Zeta = np.array(zeta)
		predY = np.dot(X, Zeta)
		R2 = get_r2(y, predY)
		return 1 - R2

	def fit(self, X, y, sample_weights=None):
		## FIX THIS TO DEAL WITH ALL 0 COLUMNS
		K = np.shape(X)[1]
		# if only one valid model set zeta = 1
		if (K == 1):
			self.coef_ = np.array([1])
			return self
		elif (K == 2):
			# initialize zeta list; all models weighted equally
			zeta_0 = np.full(K-1, 1/K)
			bnds = tuple([(0, 1)] * (K-1))
			# minimize loss function
			self.fit_res_ = minimize(
				self._loss_func,
				zeta_0,
				args=(X, y), 
				bounds=bnds,
				tol=self.tol,
				method=self.min_method)
			zeta = self.fit_res_.x
			self.coef_ = np.array([*zeta, 1 - np.sum(zeta)])
		else:
			zeta_0 = np.full(K, 1/K)
			bnds = tuple([(0, 1)] * K)
			cons = ({'type': 'eq', 'fun': lambda x:  1 - sum(x)})			
			# minimize loss function
			self.fit_res_ = minimize(
				self._loss_func,
				zeta_0,
				args=(X, y), 
				bounds=bnds,
				tol=self.tol,
				method=self.min_method,
				constraints=cons)
			zeta = self.fit_res_.x
			self.coef_ = np.array(zeta)
		return self

	def predict(self, X):
		return np.dot(X, self.coef_)

	def score(self, X, y):
		return get_r2(y, self.predict(X))

# stacking regressor
class WeightStackingRegressor(StackingRegressor):
	def __init__(self, estimators, final_estimator=ZetasEstimator(), *, cv=None,
		n_jobs=None, passthrough=False, verbose=0):
		super().__init__(
			estimators=estimators,
			final_estimator=final_estimator,
			cv=cv,
			n_jobs=n_jobs,
			passthrough=passthrough,
			verbose=verbose)

	def fit(self, X, y, sample_weights=None):
		self.ks_ = [int(est[0]) for est in self.estimators]
		self = super().fit(X, y, sample_weights)
		self.zetas_ = {'Z'+self.estimators[k][0]: self.final_estimator_.coef_[k] for k in self.ks_}
		return self

	def score(self, X, y):
		return get_r2(y, self.predict(X))

	def r2_pval(self, X, y):
		return get_r2(y, self.predict(X), pval=True)

	def avg_r2_cv(self, X, y):
		return sum(cross_val_score(self, X, y)) / 5

	def out_vals(self, X, y):
		return merge_dicts([est[1].fit().out_vals(X,y) for est in self.estimators])

def sort_tabix_output(temp_path, out_path):
	try:
		call_args = [
			tg.get_abs_path(args.SR_TWAS_dir) + '/sort_tabix_output.sh', 
			temp_path, 
			out_path]

		subprocess.check_call(
			call_args,
			cwd=tg.get_abs_path(args.out_dir),
			stdout=subprocess.DEVNULL)

	except subprocess.CalledProcessError as err:
		raise err

#############################################################
# check input arguments
Kweights = len(args.w_paths)
if Kweights < 2:
	raise SystemExit('Must specify at least 2 weight files for stacked regression.\n')

if args.genofile_type == 'vcf':
	if (args.data_format != 'GT') and (args.data_format != 'DS'):
		raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')
		
elif args.genofile_type == 'dosage':
	args.data_format = 'DS'

else:
	raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')

# out_weight_path = args.out_dir + '/temp_' + args.out_weight_file
# out_info_path = args.out_dir + '/' + args.out_info_file

# out_naive_weight_path = args.out_dir + '/temp_' + args.out_naive_weight_file
# out_naive_info_path = args.out_dir + '/' + args.out_naive_info_file

out_weight_path = out_sub_dir + '/temp_' + args.out_weight_file
out_info_path = out_sub_dir + '/' + args.out_info_file

out_naive_weight_path = out_sub_dir + '/temp_' + args.out_naive_weight_file
out_naive_info_path = out_sub_dir + '/' + args.out_naive_info_file

do_maf_diff = 0 if not args.maf_diff else 1

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
Output training info file for SR-TWAS: {out_info}
Output trained weights file for SR-TWAS: {out_weight}
{naive_info}
{naive_weight}
********************************'''.format(
	**args.__dict__,
	cvR2_str1 = {0:'Skipping evaluation of', 1:'Evaluating'}[args.cvR2],
	cvR2_str2 = {0:'', 1:' with inclusion threshold Avg.CVR2 >' + str(args.cvR2_threshold)}[args.cvR2],
	maf_diff_str1 = {0:'Not e', 1:'E'}[do_maf_diff],
	maf_diff_str2 = {0:'by MAF difference.', 1:'if MAF difference exceeds: |' + str(args.maf_diff) + '|'}[do_maf_diff],
	w_paths_str = '\n  '.join(args.w_paths),
	K = Kweights,
	out_info = out_info_path,
	out_weight = out_weight_path,
	naive_info = {0:'', 1:'Output training info file for naive method: ' + out_naive_info_path}[args.naive],
	naive_weight = {0:'', 1:'Output trained weights file for naive method: ' + out_naive_weight_path}[args.naive]))

tg.print_args(args)

##########################
# STUFF FOR TEST FILES
##########################
# Gene Expression header, sampleIDs
print('Reading genotype, expression file headers, sample IDs.\n')
sampleID, sample_size, geno_info, exp_info = tg.sampleid_startup(**args.__dict__)

# Read in gene expression info
print('Reading gene expression data.\n')
GeneExp, TargetID, n_targets = tg.read_gene_annot_exp2(annot_path=args.annot_path, **exp_info)

# get info for the weight files
weights_info = tg.weight_k_files_info2(target=TargetID[0], **args.__dict__)

ES_cols = [ES(k) for k in range(Kweights)]

# set up output files

### SR-TWAS OUTPUT
# print('Creating file: ' + out_weight_path + '\n')
out_weight_cols = ['CHROM','POS','ID','REF','ALT','TargetID','MAF','p_HWE','ES']
# pd.DataFrame(columns=out_weight_cols).to_csv(
# 	out_weight_path,
# 	header=True,
# 	index=None,
# 	sep='\t',
# 	mode='w')

# dictionary of dtypes for info output
info_wk_dtypes = merge_dicts([{W(k)+'_N_SNP':np.int64, W(k)+'_CVR2':np.float64, W(k)+'_R2':np.float64, W(k)+'_PVAL':np.float64} for k in range(Kweights)])
info_wk_cols = list(info_wk_dtypes.keys())

info_zeta_dtypes = {'Z'+str(k): np.float64 for k in range(Kweights)}
info_zeta_cols = list(info_zeta_dtypes.keys())

info_dtypes = {**info_zeta_dtypes, **info_wk_dtypes}

print('Creating file: ' + out_info_path + '\n')
out_info_cols = ['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName','sample_size','N_SNP','N_EFFECT_SNP','CVR2', 'R2', 'PVAL'] + info_zeta_cols + info_wk_cols

pd.DataFrame(columns=out_info_cols).to_csv(
	out_info_path,
	sep='\t',
	header=True,
	index=None,
	mode='w')

if args.naive:
	### NAIVE OUTPUT
	# print('Creating file: ' + out_naive_weight_path + '\n')
	# pd.DataFrame(columns=out_weight_cols).to_csv(
	# 	out_naive_weight_path,
	# 	header=True,
	# 	index=None,
	# 	sep='\t',
	# 	mode='w')

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

def get_Geno_Expr():
	target = TargetID[0]
	Expr = GeneExp.iloc[[0]]

	start = str(max(int(Expr.GeneStart) - args.window, 0))
	end = str(int(Expr.GeneEnd) + args.window)

	# # Query files; function stops here if no data in genotype and/or no weight file data
	# tabix_target_ks, empty_target_ks = tg.tabix_query_files(start, end, **args.__dict__)

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
	Expr = tg.center(Expr, sampleID)

	# return start, end, tabix_target_ks, empty_target_ks, Geno

	return start, end, Geno, Expr

##############################################################
# thread function
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Expr = Expr_orig.copy()
	Expr['TargetID'] = target

	# # Query files; function stops here if no data in genotype and/or no weight file data
	# tabix_target_ks = tabix_target_ks_orig
	# empty_target_ks = empty_target_ks_orig
	tabix_target_ks = [k for k in range(Kweights)]
	empty_target_ks = []

	# # READ IN AND PROCESS GENOTYPE DATA 
	Geno = Geno_orig.copy()

	# read in weight data for target from all K weight files
	Weights = pd.DataFrame()
	target_ks = []

	for k in tabix_target_ks:
		# read in from tabix, may be empty
		W_k = tg.read_tabix(start, end, target=target, raise_error=False, **weights_info[k])

		# time.sleep(1)

		try:
			W_k['snpIDflip'] = tg.get_snpIDs(W_k, flip=True)

			W_k = W_k[np.any(W_k[['snpID','snpIDflip']].isin(Geno.snpID.values), axis=1)]
			W_k = W_k.reset_index(drop=True)

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

		except:
			if W_k.empty:
				empty_target_ks.append(k)

		finally:
			if not W_k.empty:
				if Weights.empty:
					Weights = W_k
				else:
					Weights = Weights.merge(W_k, how='outer', on=['snpID'])
			# print(W_k)
			W_k = pd.DataFrame()


	if Weights.empty:
		print('No cis-eQTL weights with snp overlap in genotype data for target.\n')
		return None

	# ks may be added in read in step, need to sort for correct output later
	empty_target_ks.sort()

	# print(Weights)

	# add nan columns for weight files without data
	if args.maf_diff:
		empty_wk_cols = flatten([[ES(k),MAF(k)] for k in empty_target_ks])
	else: 
		empty_wk_cols = [ES(k) for k in empty_target_ks]
	Weights[empty_wk_cols] = np.nan

	# print(Weights)

	# get overlapping snps, filter Geno
	# print('Geno: ' + str(Geno.index.size))
	Geno = Geno[Geno.snpID.isin(Weights.snpID.values)].reset_index(drop=True)
	Weights = Weights[Weights.snpID.isin(Geno.snpID.values)].reset_index(drop=True)

	# print('Geno: ' + str(Geno.index.size))
	# merge datasets (inner in case missing snpIDs after filtering)
	# print('Geno dtypes')
	# print(Geno.dtypes)
	# print('Weights dtypes')
	# print(Weights.dtypes)
	Train = Geno.merge(Weights,
		left_on='snpID', 
		right_on='snpID', 
		how='inner')
	# print('Train dtypes')
	# print(Train.dtypes)
	Train = Train.set_index('snpID')
	# print('Train: ' + str(Train.index.size))
	# print('Train dtypes')
	# print(Train.dtypes)

	# filter by MAF diff - if an MAF# differs too much, set ESk to nan
	# if args.maf_diff:
	# 	for k in range(Kweights):
	# 		maf_diff = np.abs(Train[MAF(k)].values - Train['MAF'].values)
	# 		Train.loc[maf_diff > args.maf_diff, ES(k)] = np.nan
	# print(Train)
	# filter out snps where all weights are nan0[]
	all_missing_weights = Train[ES_cols].count(axis=1) == 0
	# print(all_missing_weights)
	Train = Train[~all_missing_weights]
	# print(Train)
	n_snps = Train.index.size

	if not n_snps:
		# print(Weights)
		# print(Train)
		print('No valid SNPs.\n')
		# print('All SNP MAFs for training data and testing data differ by a magnitude greater than ' + str(args.maf_diff) + '.\n')
		return None

	# do stacked regression
	print('Performing Stacked Regression.')
	X = Train[sampleID].T
	Y = Expr[sampleID].values.ravel()

	## SET UP MODELS
	if len(target_ks) == 1:
		# only one weight model with usable data, don't need to do stacking
		print('Only trained model ' + str(target_ks[0]) + ' has usable weights for target.')
		reg = WeightEstimator(Train[ES(target_ks[0])], W(target_ks[0])).fit(X, Y)
		reg.zetas_ = {'Z'+str(target_ks[0]):1}

	else:
		# do stacking
		weight_estimators = [(str(k), WeightEstimator(Train[ES(k)], W(k))) for k in target_ks]
		reg = WeightStackingRegressor(estimators=weight_estimators).fit(X, Y)

	## 5-FOLD CV
	if args.cvR2:
		print('Running 5-fold CV.')
		avg_r2_cv = reg.avg_r2_cv(X, Y)

		print('Average R2 for 5-fold CV: {:.4f}'.format(avg_r2_cv))

		if avg_r2_cv < args.cvR2_threshold:
			print('Average R2 < ' + str(args.cvR2_threshold) + '; Skipping SR training for target.')

	else: 		
		print('Skipping evaluation by 5-fold CV average R2...')
		avg_r2_cv = 0

	# needed for naive method
	wk_out_vals = final_wk_out_vals(reg.out_vals(X, Y), Kweights)

	## OUTPUT
	if (avg_r2_cv >= args.cvR2_threshold) or (not args.cvR2):	
		# output values
		stacked_r2, stacked_pval = reg.r2_pval(X, Y)
		zetas = final_zetas(reg.zetas_, Kweights)

		# output Trained weights
		Weight_Out = Train[['CHROM','POS','REF','ALT','MAF','p_HWE']].copy()
		Weight_Out['ID'] = Weight_Out.index
		Weight_Out['TargetID'] = target
		Weight_Out['ES'] = np.dot(Train[ES_cols].fillna(0), list(zetas.values()))

		# filter for non-zero effect size
		Weight_Out = Weight_Out[Weight_Out['ES'] != 0]

		# reorder columns
		Weight_Out = Weight_Out[out_weight_cols]

		pd.DataFrame(columns=out_weight_cols).to_csv(
			out_sub_dir + '/temp_' + args.out_weight_file + '_' + target + '.txt',
			header=True,
			index=None,
			sep='\t',
			mode='w')
		Weight_Out.to_csv(
			out_sub_dir + '/temp_' + args.out_weight_file + '_' + target + '.txt',
			sep='\t',
			header=None,
			index=None,
			mode='a')
		sort_tabix_output(
			out_sub_dir + '/temp_' + args.out_weight_file + '_' + target + '.txt', 
			out_sub_dir + '/' + args.out_weight_file + '_' + target + '.txt')
		
		# output training info
		Info = Expr[['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()
		Info['sample_size'] = Y.size
		Info['N_SNP'] = n_snps
		Info['N_EFFECT_SNP'] = Weight_Out.ES.size
		Info['CVR2'] = avg_r2_cv
		Info['R2'] = stacked_r2
		Info['PVAL'] = stacked_pval
		Info = Info.join(pd.DataFrame.from_records({**zetas, **wk_out_vals}, index=[0]).astype(info_dtypes))[out_info_cols]
		Info.to_csv(
			out_info_path,
			sep='\t',
			header=None,
			index=None,
			mode='a')


	## NAIVE MODEL
	if args.naive:
		# SET UP MODELS
		# if len(target_ks) == 1:
		# 	Train['Naive'] = Train[ES(target_ks[0])]
		# 	naive_reg = reg
		# else:
		# 	# Train['Naive'] = np.nanmean(Train[ES_cols].values, axis=1)
		# 	Train['Naive'] = np.mean(np.nan_to_num(Train[ES_cols].values), axis=1)
		# 	naive_reg = WeightEstimator(Train['Naive'], 'Naive').fit(X, Y)
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

			pd.DataFrame(columns=out_weight_cols).to_csv(
				out_sub_dir + '/temp_' + args.out_naive_weight_file + '_' + target + '.txt',
				header=True,
				index=None,
				sep='\t',
				mode='w')
			Naive_Weight.to_csv(
				out_sub_dir + '/temp_' + args.out_naive_weight_file + '_' + target + '.txt',
				sep='\t',
				header=None,
				index=None,
				mode='a')
			sort_tabix_output(
				out_sub_dir + '/temp_' + args.out_naive_weight_file + '_' + target + '.txt', 
				out_sub_dir + '/' + args.out_naive_weight_file + '_' + target + '.txt')

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
	start, end, Geno_orig, Expr_orig = get_Geno_Expr()
	pool = multiprocessing.Pool(args.thread)
	pool.imap(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')

	# sort tabix output
	# sort_tabix_output(out_weight_path, args.out_dir + '/' + args.out_weight_file)

	# if args.naive:
	# 	sort_tabix_output(out_naive_weight_path, args.out_dir + '/' + args.out_naive_weight_file)


###############################################################
# time calculation
elapsed_sec = time.time() - start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

###############################################################

if not args.no_log:
	sys.stdout.close()
