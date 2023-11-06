#!/usr/bin/env python

###############################################################
# import packages needed
import argparse
import operator
import os
import multiprocessing
import subprocess
import sys

from time import time

import numpy as np
import pandas as pd

###############################################################
# time calculation
start_time=time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Prediction')
# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

parser.add_argument('--chr', type=str, dest='chrm')
parser.add_argument('--format', type=str, dest='data_format', choices=['GT', 'DS'], default='GT')
parser.add_argument('--gene_anno', type=str, dest='annot_path')
parser.add_argument('--genofile', type=str, dest='geno_path')
parser.add_argument('--genofile_type', type=str, choices=['vcf', 'dosage'])
parser.add_argument('--maf_diff', type=float, default=0.2)
parser.add_argument('--missing_rate', type=float, default=0.2)
parser.add_argument('--out_dir', type=str)
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--out_pred_file', type=str, default='')
parser.add_argument('--sub_dir', type=str, default='')
parser.add_argument('--test_sampleID', type=str, dest='sampleid_path')
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--weight', type=str, dest='w_path')
parser.add_argument('--window', type=int, default=1000000)
args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
import TIGARutils as tg

#######################################################################

###############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm + '_Pred'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_pred_file:
	args.out_pred_file = args.out_prefix + '_GReX.txt'

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, args.sub_dir)
else:
	out_sub_dir = args.out_dir

# Make output, log directories
os.makedirs(out_sub_dir, exist_ok=True)
os.makedirs(os.path.join(args.out_dir, 'logs'), exist_ok=True)

# set stdout to log
sys.stdout = open(os.path.join(args.out_dir, 'logs', args.log_file), 'w')

# Check tabix command
try:
	subprocess.check_call(['which','tabix'], stdout=subprocess.DEVNULL)
except:
	raise SystemExit('Error: Required tool TABIX is not available.\n')

# Check gene expression file
try:
	os.path.isfile(args.annot_path)
except:
	raise SystemExit('Error: Gene annotation file does not exist.')

# Check genotype file 
try:
	os.path.isfile(args.geno_path)
except:
	raise SystemExit('Error: Training genotype file does not exist.')

# Check training sample ID file
try:
	os.path.isfile(args.sampleid_path)
except:
	raise SystemExit('Error: Test sample ID file does not exist.')

# Check weight file
try:
	os.path.isfile(args.w_path)
except:
	raise SystemExit('Error: Weight file does not exist.')


#######################################################################
# Input Arguments for GReX Prediction

# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the testing genotype file (bgzipped and tabixed) 
# --genofile_type: Genotype file type: 'vcf' or 'dosage'
# --genofile_colnames: File with column heads of genotype file
# --format: Genotype format in VCF file that should be used: 'GT' (default) for genotype data or 'DS' for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --maf_diff: MAF difference threshold for matching SNPs from eQTL weight file and test genotype file. If SNP MAF difference is greater than maf_diff (default 0.2), , the SNP will be excluded
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)

###############################################################
# check input arguments
if args.genofile_type == 'vcf':
	if (args.data_format != 'GT') and (args.data_format != 'DS'):
		raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')
		
elif args.genofile_type == 'dosage':
	args.data_format = 'DS'

else:
	raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')

# out_pred_path = args.out_dir + '/' + args.out_pred_file
out_pred_path = out_sub_dir + '/' + args.out_pred_file

###############################################################
# Print input arguments
print(
'''********************************
Input Arguments
Gene annotation file specifying genes for prediction: {annot_path}
Prediction sampleID file: {sampleid_path}
Chromosome: {chrm}
cis-eQTL weight file: {w_path}
Prediction genotype file: {geno_path}
Genotype file used for prediction is type: {genofile_type}
Genotype data format: {data_format}
Gene prediction region SNP inclusion window: +-{window}
Excluding SNPs if missing rate exceeds: {missing_rate}
Excluding SNPs matched between eQTL weight file and prediction genotype file if MAF difference exceeds: |{maf_diff}|
Number of threads: {thread}
Output directory: {out_dir}
Output prediction results file: {out_path}
********************************'''.format(
	**args.__dict__,
	out_path = out_pred_path))

# tg.print_args(args)

###############################################################

# Load genotype column names of test genotype file
sampleID, sample_size, geno_info = tg.sampleid_startup(**args.__dict__)

# Load annotation file
print('Reading gene annotation file.')
Gene, TargetID, n_targets = tg.read_gene_annot_exp(**args.__dict__)

# get weight file info
# weight_info = tg.weight_file_info(add_cols=['MAF'], **args.__dict__)
try:
	weight_info = tg.weight_file_info2(target=TargetID[0], add_cols=['MAF'], **args.__dict__)
except:
	weight_info = tg.weight_file_info2(target=TargetID[0], **args.__dict__)

print('Creating output file: ' + out_pred_path + '\n')
out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName', *sampleID]
pd.DataFrame(columns=out_cols).to_csv(
	out_pred_path,
	sep='\t', 
	index=None, 
	header=True, 
	mode='w')

print('********************************\n')


def get_Geno():
	target = TargetID[0]
	Gene_info = Gene.iloc[[0]]

	start = str(max(int(Gene_info.GeneStart) - args.window, 0))
	end = str(int(Gene_info.GeneEnd) + args.window)

	print('Reading genotype data.')
	Geno = tg.read_tabix(start, end, sampleID, **geno_info)

	# filter out variants that exceed missing rate threshold
	Geno = tg.handle_missing(Geno, sampleID, args.missing_rate)

	# calculate MAF
	Geno = tg.calc_maf(Geno, sampleID, 0, op=operator.ge)
	Geno['MAF_test'] = Geno.MAF

	# get flipped snpIDs
	Geno_flip = Geno.copy()
	Geno_flip['snpID'] = tg.get_snpIDs(Geno_flip, flip=True)
	Geno_flip['MAF_test'] = 1 - Geno_flip.MAF
	Geno_flip[sampleID] = 2 - Geno_flip[sampleID]

	# append flipped dataframe and sort by position
	Geno = Geno.append(Geno_flip).sort_values(by=['POS'])

	# drop unneeded columns
	Geno = Geno.drop(columns=['CHROM','POS','REF','ALT','MAF'])

	# center data
	Geno = tg.center(Geno, sampleID)

	return start, end, Geno


###############################################################
# thread function
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Gene_info = Gene.iloc[[num]].reset_index(drop=True)

	print('Getting weight data for target.')
	if args.maf_diff:
		Weight = tg.read_tabix(start, end, target=target, do_zgrep=True, **weight_info)[['snpID','ES','MAF']]
	else:
		Weight = tg.read_tabix(start, end, target=target, do_zgrep=True, **weight_info)[['snpID','ES']]

	Geno = Geno_orig.copy()

	snp_overlap = np.intersect1d(Weight.snpID, Geno.snpID)

	if not snp_overlap.size:
		print('No overlapping test SNPs between weight and genotype file for TargetID: ' + target + '\n')
		return None

	# filter out non-matching snpID rows
	Weight = Weight[Weight.snpID.isin(snp_overlap)]
	Geno = Geno[Geno.snpID.isin(snp_overlap)].reset_index(drop=True)

	# merge Geno, Weight
	Pred = Geno.merge(
		Weight, 
		left_on='snpID', 
		right_on='snpID', 
		how='inner')

	if args.maf_diff:
		Pred['diff'] = np.abs(Pred['MAF'].astype('float') - Pred['MAF_test'].astype('float'))
		
		Pred = Pred[Pred['diff'] <= args.maf_diff].drop(columns=['MAF','MAF_test','diff']).reset_index(drop=True)

		if Pred.empty:
			print('All SNP MAFs for training data and testing data differ by a magnitude greater than ' + str(args.maf_diff) + ' for TargetID: ' + target + '\n')
			return None

	print('Predicting GReX.\nN SNPs=' + str(Pred.snpID.size))

	# output results
	Pred_GReX = pd.DataFrame(
		data=np.dot(Pred[sampleID].T, Pred['ES'].values),
		index=sampleID).T

	Result = pd.concat([Gene_info, Pred_GReX], axis=1)

	Result.to_csv(
		out_pred_path,
		sep='\t',
		index=None,
		header=None,
		mode='a')

	print('Target prediction completed.\n')   

###############################################################
# thread process
if __name__ == '__main__':
	print('Starting prediction for ' + str(n_targets) + ' target genes.\n')
	start, end, Geno_orig = get_Geno()
	pool = multiprocessing.Pool(args.thread)
	pool.imap(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')

###############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)
