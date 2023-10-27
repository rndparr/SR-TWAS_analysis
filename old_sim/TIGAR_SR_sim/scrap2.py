
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys

from time import time

import numpy as np
import pandas as pd

import scipy.stats as stats
from sklearn.model_selection import KFold
import statsmodels.api as sm

parser = argparse.ArgumentParser(description='DPR Training')

# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

parser.add_argument('--chr', type=str, dest='chrm')
## 0 do not run cvR2, 1 run cvR2
parser.add_argument('--cvR2', type=int, choices=[0, 1], default=1)
parser.add_argument('--cvR2_threshold', type=float, default=0.005)
## '1' (Variational Bayesian) [default]; '2' (MCMC)
parser.add_argument('--dpr', type=str, choices=['1', '2'], default=1)
## 'additive' (fixed + random); 'fixed' (fixed effects) [default]
parser.add_argument('--ES', type=str, choices=['additive', 'fixed'], default='fixed')
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
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path')
parser.add_argument('--window', type=int, default=1000000)
args = parser.parse_args()

args.geno_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/GTEx_ABCA7_raw.dosage.gz'
args.out_info_file = 'GTEx_train_info.txt'
args.hwe = 1e-05
args.cvR2_threshold = 0.005
args.data_format = 'DS'
args.maf = 0.01
args.dpr = '1'
args.sub_dir = 0
args.out_weight_file = 'GTEx_train_weight'
args.log_file = 'GTEx_train_log.txt'
args.missing_rate = 0.2
args.ES = 'fixed'
args.sampleid_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/GTEx_train_465_sampleid.txt'
args.geneexp_path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/expression/raw_dosage/GTEx_expr.txt'
args.genofile_type = 'dosage'
args.TIGAR_dir = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/TIGAR_SR_sim/'
args.out_prefix = 'CHR19_DPR_train'
args.out_dir = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/train//sims'
args.cvR2 = 1
args.thread = 1
args.window = 1000000
args.chrm = '19'
args.job_suf = 'GTEx'


sys.path.append(args.TIGAR_dir)
import TIGARutils as tg

path = '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/train/sims/GTEx_train_weight_0.001_0.1_0.1_1.txt.gz'

target = '0.001_0.1_0.1_1'

args.w_paths = ['/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/train/sims/GTEx_train_weight']
weights_info = tg.weight_k_files_info2(target=target, **args.__dict__)


weights_info


weights_info[0]['col_inds']

df = pd.read_csv(
	path,
	sep='\t',
	low_memory=False,
	header=0,
	usecols=weights_info[0]['col_inds'],
	names=weights_info[0]['cols'],
	dtype=weights_info[0]['dtype'],
	compression='gzip')


