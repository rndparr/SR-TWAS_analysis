#!/usr/bin/env python

#############################################################
# Import packages needed
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys

from time import time

import numpy as np
import pandas as pd

# import grid search for model selection
from sklearn.model_selection import KFold
# For Elastic Net Regression
from sklearn.linear_model import ElasticNetCV
# For OLS regression in cross validation
from scipy import stats
import statsmodels.api as sm

from warnings import filterwarnings
filterwarnings('ignore')

#############################################################
# time calculation
start_time = time()

#############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='EN Training')

# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

parser.add_argument('--chr', type=str, dest='chrm')
## 0 do not run cvR2, 1 run cvR2
parser.add_argument('--cvR2', type=int, choices=[0, 1], default=1)
parser.add_argument('--cvR2_threshold', type=float, default=0.005)

parser.add_argument('--alpha', type=float, default=0.5, metavar='0.5', nargs='+', help='ratio of L1 and L2 for EN model training (default: 0.5)')
parser.add_argument('--use_alpha', type=int, choices=[0, 1], default=0, 
		help='use specified alpha value? (0: [0.1, 0.5, 0.9, 1] [default], 1: --alpha value)')
parser.add_argument('--cv', type=int, dest='cv', default=5, 
		metavar=5, help='k for k-fold cross validation for EN model training (default: 5)')


parser.add_argument('--format', type=str, dest='data_format', choices=['GT', 'DS'], default='GT')
parser.add_argument('--gene_exp', type=str, dest='geneexp_path')
parser.add_argument('--genofile', type=str, dest='geno_path')
parser.add_argument('--genofile_type', type=str, choices=['vcf', 'dosage'])
parser.add_argument('--hwe', type=float, default=0.00001)
parser.add_argument('--job_suf', type=str, default='')
parser.add_argument('--maf', type=float, default=0.01)
parser.add_argument('--missing_rate', type=float, default=0.2)
parser.add_argument('--out_dir', type=str)
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--out_info_file', type=str, default='')
parser.add_argument('--out_weight_file', type=str, default='')
parser.add_argument('--sub_dir', type=str, default='')
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path')
parser.add_argument('--window', type=int, default=1000000)
args = parser.parse_args()
sys.path.append(args.TIGAR_dir)

###############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm + '_EN_train'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_info_file:
	args.out_info_file = args.out_prefix + '_GeneInfo.txt'

if not args.out_weight_file:
	args.out_weight_file = args.out_prefix + '_eQTLweights.txt'

if not args.job_suf:
	args.job_suf = '_CHR' + args.chrm

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, args.sub_dir)
else:
	out_sub_dir = args.out_dir

# Check tabix command
try:
	subprocess.check_call(['which','tabix'], stdout=subprocess.DEVNULL)
except:
	raise SystemExit('Error: Required tool TABIX is not available.\n')

# Check gene expression file
try:
	os.path.isfile(args.geneexp_path)
except:
	raise SystemExit('Error: Gene expression file does not exist.')

# Check genotype file 
try:
	os.path.isfile(args.geno_path)
except:
	raise SystemExit('Error: Training genotype file does not exist.')

# Check training sample ID file
try:
	os.path.isfile(args.sampleid_path)
except:
	raise SystemExit('Error: Training sample ID file does not exist.')

# Make output, log directories
os.makedirs(out_sub_dir, exist_ok=True)
os.makedirs(os.path.join(args.out_dir, 'logs'), exist_ok=True)

# set stdout to log
sys.stdout = open(os.path.join(args.out_dir, 'logs', args.log_file), 'w')

#############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

args.alpha = [0.1, 0.5, 0.9, 1] if not args.use_alpha else args.alpha

# set output file names
if not args.out_prefix:
	if args.chrm is not None:
		args.out_prefix = 'CHR' + args.chrm + '_EN_train'
	else:
		args.out_prefix = 'EN_train'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_info_file:
	args.out_info_file = args.out_prefix + '_GeneInfo.txt'

if not args.out_weight_file:
	args.out_weight_file = args.out_prefix + '_eQTLweights.txt'

# sub-directory in out directory
if args.sub_dir:
	args.out_sub_dir = os.path.join(args.out_dir, 'EN_CHR' + args.chrm)
else:
	args.out_sub_dir = args.out_dir
args.out_sub_dir = tg.get_abs_path(args.out_sub_dir)

args.tmp_weight_path = args.out_sub_dir + '/temp_' + args.out_weight_file
args.out_weight_path = args.out_sub_dir + '/' + args.out_weight_file
args.out_info_path = args.out_sub_dir + '/' +  args.out_info_file

#############################################################
# print to log
print('''********************************
Input Arguments
{chrm_str}
Gene Annotation and Expression file: {geneexp_path}
Training sampleID file: {sampleid_path}
Training genotype file: {geno_path}
Genotype file used for training is type: {genofile_type}
Genotype data format: {data_format}
Gene training region SNP inclusion window: +-{window}
Excluding SNPs if missing rate exceeds: {missing_rate}
MAF threshold for SNP inclusion: {maf}
HWE p-value threshold for SNP inclusion: {hwe}
{cvR2_str1} Elastic-Net model by 5-fold cross validation{cvR2_str2}.
Number of cross-validation folds used to tune Elastic-Net penalty parameter (lambda): {cv}
Ratio for L1 & L2 penalty used by Elastic-Net regression (alpha): {alpha}
Number of threads: {thread}
Output directory: {out_dir}
Output training weights file: {out_weight_path}
Output training info file: {out_info_path}\n
********************************\n'''.format(
	**args.__dict__,
	cvR2_str1 = {0:'Skipping evaluation of', 1:'Evaluating'}[args.cvR2],
	cvR2_str2 = {0:'', 1:' with inclusion threshold Avg.CVR2 >' + str(args.cvR2_threshold)}[args.cvR2],
	chrm_str = '\nChromosome: ' + args.chrm if not args.chrm is None else ''
	))

tg.print_args(args)

#############
# FUNCTIONS

def elastic_net(train, test=None, k=args.cv, Alphas=args.alpha):
	train = train.copy()
	trainX = train.iloc[:,0:-1]
	trainY = train.iloc[:,-1]

	if test is not None:
		test = test.copy()
		testX = test.iloc[:,0:-1]
		testY = test.iloc[:,-1]

	else:
		testX = trainX
		testY = trainY

	reg = ElasticNetCV(
		l1_ratio=Alphas,
		fit_intercept=False,
		alphas=np.arange(0, 1.01, 0.05),
		selection='random',
		cv=k).fit(trainX,trainY)
	# reg = ElasticNetCV(
	# 	l1_ratio=Alphas,
	# 	fit_intercept=False,
	# 	selection='random',
	# 	cv=k).fit(trainX,trainY)

	Alpha = reg.l1_ratio_
	Lambda = reg.alpha_
	cvm = np.min(reg.mse_path_)
	beta = reg.coef_

	# print('alphas: ' + str(reg.alphas_))

	predY = reg.predict(testX)

	lm = sm.OLS(testY, sm.add_constant(predY)).fit()

	Rsquared = lm.rsquared

	if test is not None:
		return Rsquared

	Pvalue = lm.f_pvalue

	return beta, Rsquared, Pvalue, Alpha, Lambda, cvm


# function to do the ith cross validation step
def do_cv(i, Geno_Exp_df, cv_trainID, cv_testID):
	Geno_Exp_df = Geno_Exp_df.copy()
	train_geno_exp = Geno_Exp_df.loc[cv_trainID[i]].dropna()
	test_geno_exp = Geno_Exp_df.loc[cv_testID[i]].dropna()

	cv_rsquared = elastic_net(train_geno_exp, test_geno_exp)
	return cv_rsquared

def sort_tabix_output(temp_path, out_path):
	try:
		call_args = [
			tg.get_abs_path(args.TIGAR_dir) + '/sort_tabix_output.sh', 
			temp_path, 
			out_path]

		subprocess.check_call(
			call_args,
			cwd=tg.get_abs_path(args.out_dir),
			stdout=subprocess.DEVNULL)

	except subprocess.CalledProcessError as err:
		raise err




#############################################################
# Startup for training jobs: get column header info, sampleIDs
print('Reading genotype, expression file headers, sample IDs.\n')
sampleID, sample_size, geno_info, exp_info = tg.sampleid_startup(**args.__dict__)

# Read in gene expression info
print('Reading gene expression data.\n')
GeneExp, TargetID, n_targets = tg.read_gene_annot_exp(**exp_info)

# PREP CROSS VALIDATION SAMPLES - Split sampleIDs for cross validation
if args.cvR2:
	print('Splitting sample IDs randomly for 5-fold cross validation by average R2...\n')
	kf = KFold(n_splits=5)
	kf_splits = [(sampleID[x], sampleID[y]) for x,y in kf.split(sampleID)]
	CV_trainID, CV_testID = zip(*kf_splits)

else:
	print('Skipping splitting samples for 5-fold cross validation...\n')
	CV_trainID, CV_testID = None, None

# PREP OUTPUT - print output headers to files
# print('Creating file: ' + args.tmp_weight_path + '\n')
print('Creating file: ' + args.out_weight_path + '_header.txt' + '\n')
weight_cols = ['CHROM','POS','snpID','REF','ALT','TargetID','MAF','p_HWE','ES']
pd.DataFrame(columns=weight_cols).to_csv(
	args.out_weight_path + '_header.txt.gz',
	header=True,
	index=None,
	sep='\t',
	mode='w',
	compression='gzip')


print('Creating file: ' + args.out_info_path + '\n')
info_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainPVALUE','TrainR2','k-fold','alpha','Lambda','cvm','CVR2_threshold']
pd.DataFrame(columns=info_cols).to_csv(
	args.out_info_path,
	header=True,
	index=None,
	sep='\t',
	mode='w')

print('********************************\n\nStarting Elastic-Net training for ' + str(n_targets) + ' target genes.\n')

#############################################################


def get_Geno():
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

	# return start, end, tabix_target_ks, empty_target_ks, Geno

	return start, end, Geno['snpID'].size, Geno


##############################################################
# thread function
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Expr = GeneExp.iloc[[num]]

	# center data
	Expr = tg.center(Expr, sampleID)

	# # READ IN AND PROCESS GENOTYPE DATA 
	Geno = Geno_orig.copy()

	# merge geno, expression files, transpose
	Geno_Exp = pd.concat([
		Geno.set_index(['snpID'])[sampleID],
		Expr.set_index(['TargetID'])[sampleID]
		]).T

	# 5-FOLD CROSS-VALIDATION
	# Evaluate whether Elastic Net model valid
	if args.cvR2:
		print('Running 5-fold CV.')
		# print('Starting with Elastic-Net penalty parameter')
		do_cv_args = [Geno_Exp, CV_trainID, CV_testID]
		k_fold_R2 = [do_cv(i, *do_cv_args) for i in range(5)]

		avg_r2_cv = sum(k_fold_R2) / 5

		print('Average R2 for 5-fold CV: {:.4f}'.format(avg_r2_cv))

		if avg_r2_cv < args.cvR2_threshold:
			print('Average R2 < ' + str(args.cvR2_threshold) + '; Skipping Elastic-Net training for TargetID: ' + target + '\n')
			return None

	else:
		avg_r2_cv = 0
		print('Skipping evaluation by 5-fold CV average R2...')

	# FINAL MODEL TRAINING
	print('Running Elastic-Net training.')
	# initialize Weight dataframe
	Weight = Geno[['CHROM','POS','snpID','REF','ALT','p_HWE','MAF']].copy()

	Weight['TargetID'] = target

	# do elastic net training
	Weight['ES'], R2, Pvalue, Alpha, Lambda, cvm = elastic_net(Geno_Exp)

	# filter
	Weight = Weight[Weight['ES'] != 0]
	n_effect_snp = Weight.ES.size

	# reorder columns for output
	Weight = Weight[weight_cols]

	# output
	pd.DataFrame(columns=weight_cols).to_csv(
		out_sub_dir + '/temp_' + args.out_weight_file + '_' + target + '.txt',
		header=True,
		index=None,
		sep='\t',
		mode='w')
	Weight.to_csv(
		out_sub_dir + '/temp_' + args.out_weight_file + '_' + target + '.txt',
		header=None,
		index=None,
		sep='\t',
		mode='a')
	sort_tabix_output(
		out_sub_dir + '/temp_' + args.out_weight_file + '_' + target + '.txt',
		out_sub_dir + '/' + args.out_weight_file + '_' + target + '.txt')	

	# output training information, result from elastic net
	Info = Expr[['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()

	Info['sample_size'] = sample_size
	Info['n_snp'] = n_snp
	Info['n_effect_snp'] = n_effect_snp
	Info['CVR2'] = avg_r2_cv
	Info['TrainPVALUE'] = Pvalue if not np.isnan(Pvalue) else 'NaN'
	Info['TrainR2'] = R2 if n_effect_snp else 0
	Info['k_fold'] = args.cv
	Info['alpha'] = Alpha
	Info['lambda'] = Lambda
	Info['cvm'] = cvm
	Info['CVR2_threshold'] = args.cvR2_threshold if args.cvR2 else 0

	Info.to_csv(
		args.out_info_path,
		header=None,
		index=None,
		sep='\t',
		mode='a')

	print('Target Elastic-Net training completed.\n')



if __name__ == '__main__':
	# print('Starting EN training for ' + str(n_targets) + ' target genes.\n')
	start, end, n_snp, Geno_orig = get_Geno()
	pool = multiprocessing.Pool(args.thread)
	pool.imap(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')

	# sort tabix output
	# sort_tabix_output(out_weight_path, args.out_dir + '/' + args.out_weight_file)



############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

sys.stdout.close()



