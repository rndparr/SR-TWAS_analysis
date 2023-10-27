
# # BASH:
# # load modules, python environment
# module load tabix/0.2.6; conda activate myenv; export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/:$PYTHONPATH; ipython

######################
# # PYTHON:
# 
# import packages
import argparse
import operator
import sys
import numpy as np
import pandas as pd
TIGAR_dir = '/mnt/icebreaker/data2/home/rparrish/github/TIGAR/'
sys.path.append(TIGAR_dir)
import TIGARutils as tg

parser = argparse.ArgumentParser(description='get dosage formats')
args = parser.parse_args()

# from pandas.io.parsers import ParserError
# from TIGARutils import *


# def read_tabix(start, end, sampleID, chrm, path, file_cols, col_inds, cols, dtype, genofile_type=None, data_format=None, target_ind=5, target=None, weight_threshold=0, raise_error=True, **kwargs):

# 	# subprocess command
# 	command_str = ' '.join(['tabix', path, chrm + ':' + start + '-' + end])

# 	proc = subprocess.Popen(
# 		[command_str],
# 		shell=True,
# 		stdout=subprocess.PIPE,
# 		bufsize=1)

# 	# initialize bytearray
# 	proc_out = bytearray()

# 	# set correct filter function by file type
# 	if genofile_type == 'vcf':
# 		bformat = str.encode(data_format)
# 		filter_line = functools.partial(filter_vcf_line, bformat=bformat, col_inds=col_inds, split_multi_GT=data_format == 'GT')
# 	elif genofile_type == 'weight':
# 		btarget = str.encode(target)
# 		filter_line = functools.partial(filter_weight_line, btarget=btarget, target_ind=target_ind, col_inds=col_inds)
# 	else:
# 		filter_line = functools.partial(filter_other_line, col_inds=col_inds)

# 	# while subprocesses running, read lines into byte array
# 	while proc.poll() is None:
# 		line = proc.stdout.readline()
# 		if len(line) == 0:
# 			break
# 		proc_out += filter_line(line)
# 	# read in lines still remaining after subprocess completes
# 	for line in proc.stdout:
# 		proc_out += filter_line(line)

# 	if not proc_out and raise_error:
# 		print('No tabix data for target.\n')
# 		raise NoTargetDataError

# 	# read data into dataframe
# 	try:
# 		df = pd.read_csv(
# 			StringIO(proc_out.decode('utf-8')),
# 			sep='\t',
# 			low_memory=False,
# 			header=None,
# 			names=cols,
# 			dtype=dtype)
# 	except ParserError:
# 		df_chunks = pd.read_csv(
# 			StringIO(proc_out.decode('utf-8')),
# 			sep='\t',
# 			chunksize=1000,
# 			iterator=True,
# 			low_memory=False,
# 			header=None,
# 			names=cols,
# 			dtype=dtype)

# 		df = pd.concat([chunk for chunk in df_chunks])

# 	# filter out rows where all sampleID values are nan
# 	if len(sampleID):
# 		df = df[df[sampleID].count(axis=1) != 0].reset_index(drop=True)

# 	df = optimize_cols(df)

# 	# get snpID, filter out duplicates
# 	if (genofile_type != 'weight') or (not 'snpID' in cols):
# 		df['snpID'] = get_snpIDs(df)
# 		df = df.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

# 	# weight file handling
# 	if (genofile_type == 'weight'):
# 		## figure out a way to handle target ids with decimal places when both files dont necessarily have those
# 		# if not np.all(df['TargetID'] == target):
# 		# partial_targetids = np.unique(df['TargetID'])

# 		# VC_TWAS uses ES = b + beta
# 		if (not 'ES' in cols) and (('b' in cols) and ('beta' in cols)):
# 			df['ES'] = df['b'] + df['beta']

# 		if weight_threshold:
# 			# filter out weights below threshold
# 			df = df[operator.gt(np.abs(df['ES']), weight_threshold)].reset_index(drop=True)

# 			if df.empty and raise_error:
# 				print('No test SNPs with cis-eQTL weights with magnitude that exceeds specified weight threshold for TargetID: ' + target + '.\n')
# 				raise NoTargetDataError

# 	# remove rows with non-valid GT values; ie those from multiallelic rows
# 	if (data_format == 'GT'):
# 		valid_GT = ['.|.', '0|0', '0|1', '1|0', '1|1', 
# 		'./.', '0/0', '0/1', '1/0', '1/1']
# 		df = df[np.all(df[sampleID].isin(valid_GT), axis=1)].reset_index(drop=True)

# 	# process dosage, vcf file sample df
# 	if (data_format == 'GT') or (data_format == 'DS'):
# 		df = reformat_sample_vals(df, data_format, sampleID)

# 	if df.empty and raise_error:
# 		print('No valid tabix data for target.\n')
# 		raise NoTargetDataError

# 	return df

args.data_format = 'GT'
args.genofile_type = 'vcf'
args.chrm = '19'
args.window = 1000000
args.maf = 0.05
args.missing_rate = 0.2
args.hwe = 0.00001

######################
# # # GTEx V8 Info:
# args.sampleid_path = '/mnt/YangFSS/data/rparrish/GTEx_V8/RefCovSampIDs.txt'
# args.geno_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/GTEx_ABCA7_raw.vcf.gz'
# out_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/GTEx_ABCA7_raw.dosage'
######################
# # ROSMAP Info:
# args.sampleid_path = '/mnt/YangFSS/data2/rparrish/ROSMAP_WGS/sampleid/ROSMAP_WGS_sampleid.txt'
args.sampleid_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/ROSMAP_sampleid.txt'
args.geno_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/ROSMAP_ABCA7_raw.vcf.gz'
out_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/ROSMAP_ABCA7_raw.dosage'

######################
###
## CHROM   GeneStart   GeneEnd TargetID    GeneName
## 19  1040101 1065572 ENSG00000064687.12  ABCA7
###
# Genotype file info for both:
GeneStart = '1040101'
GeneEnd = '1065572'

## output should have header:
## CHROM     POS     ID  REF     ALT     Sample1     Sample2...

######################
# get columns from the genotype file

sampleID, sample_size, geno_info = tg.sampleid_startup(**args.__dict__)

# get positions for gene start and end with window
start = str(max(int(GeneStart) - args.window, 0))
end = str(int(GeneEnd) + args.window)

Geno = tg.read_tabix(start, end, sampleID, **geno_info)
# Geno.shape

# filter out variants that exceed missing rate threshold
Geno = tg.handle_missing(Geno, sampleID, args.missing_rate)
# Geno.shape

# get, filter maf
Geno = tg.calc_maf(Geno, sampleID, args.maf)
# Geno.shape

# get, filter p_HWE
Geno = tg.calc_p_hwe(Geno, sampleID, args.hwe)
# Geno.shape

# output header to file
out_cols = ['CHROM', 'POS', 'snpID', 'REF', 'ALT'] + sampleID.tolist()
pd.DataFrame(columns=out_cols).to_csv(
	out_path,
	sep='\t',
	header=True,
	index=None,
	mode='w')

# output genotype dosage file to path
Geno[out_cols].to_csv(
	out_path,
	sep='\t',
	header=None,
	index=None,
	mode='a')

#####
# exit python environment 
# exit

# # BASH:
# # set out_path depending on dataset
# dataset=GTEx
# dataset=ROSMAP
# out_path=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/${dataset}_ABCA7_raw.dosage; bgzip -f ${out_path} && tabix -f -b2 -e2 -S1 ${out_path}.gz

