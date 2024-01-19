#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=FALSE, digits=20)

# load libraries
library(doFuture)
library(foreach)

# load grguments
args=(commandArgs(TRUE))
cat(paste0('\nargs:\n'))
print(args)
cat()
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==1) {
	sim_dir = as.character(args[[1]])
	suffix = ''
	scenario = '0.1_0.1_0.5_0.1_0.1'
	sim_i = 667
} else {
	sim_dir = as.character(args[[1]])
	suffix = as.character(args[[2]])
	scenario = as.character(args[[3]])
	sim_i = as.integer(args[[4]])
}


#################
## functions
read_bgzip <- function(path) {
	raw <- read.delim(con <- gzfile(path, 'r'), check.names = FALSE); close(con)
	return(raw)
}

resample <- function(x, ...) x[sample.int(length(x), ...)]

#################

## set scenario
# scenario <- '0.1_0.1_0.5_0.1_0.1'
## set seed, pick a simulation
set.seed(7654321)
# sim_i <- 667

# sample(1:1000, 1)
# 177 ## no PrediXcan-GTEx model
# 667

# c('PrediXcan-GTEx',
# 'TIGAR-ROSMAP',
# 'TIGAR-ROSMAP_valid',
# 'Naive-TIGAR-ROSMAP_PrediXcan-GTEx',
# 'SR-TIGAR-ROSMAP_PrediXcan-GTEx',
# 'Avg-SRbasevalid')

#################
## get trained weights
# TIGAR-ROSMAP_valid_train_weight_0.1_0.1_0.5_0.1_0.1_190.txt.gz
weight_dir <- paste0(sim_dir, 'train/sims/')

gtex_weight_path <- paste0(weight_dir, 'PrediXcan-GTEx', suffix, '_train_weight_', scenario, '_', sim_i, '.txt.gz')
rosmap_weight_path <- paste0(weight_dir, 'TIGAR-ROSMAP', suffix, '_train_weight_', scenario, '_', sim_i, '.txt.gz')
rosmap_v_weight_path <- paste0(weight_dir, 'TIGAR-ROSMAP_valid', suffix, '_train_weight_', scenario, '_', sim_i, '.txt.gz')

gtex_weights_raw <- read_bgzip(gtex_weight_path)
rosmap_weights_raw <- read_bgzip(rosmap_weight_path)
rosmap_v_weights_raw <- read_bgzip(rosmap_v_weight_path)

gtex_weights_nsnps <- nrow(gtex_weights_raw)
rosmap_weights_nsnps <- nrow(rosmap_weights_raw)
rosmap_v_weights_nsnps <- nrow(rosmap_v_weights_raw)

## get genotype file rows
gt_dir <- paste0(sim_dir, 'genotype/')

gtex_gt_path <- paste0(gt_dir, 'GTEx_ABCA7_raw.dosage.gz')
rosmap_gt_path <- paste0(gt_dir, 'ROSMAP_ABCA7_raw.dosage.gz')

gtex_gt <- read_bgzip(gtex_gt_path)[, 1:5]
rosmap_gt <- read_bgzip(rosmap_gt_path)[, 1:5]

gtex_gt_nsnps <- nrow(gtex_gt)
rosmap_gt_nsnps <- nrow(rosmap_gt)


## get final weight vector with same length as available snps
gtex_weights <- c(gtex_weights_raw[, 'ES'], rep(0, gtex_gt_nsnps - gtex_weights_nsnps))
rosmap_weights <- c(rosmap_weights_raw[, 'ES'], rep(0, rosmap_gt_nsnps - rosmap_weights_nsnps))
rosmap_v_weights <- c(rosmap_v_weights_raw[, 'ES'], rep(0, rosmap_gt_nsnps - rosmap_v_weights_nsnps))

#################
### output directory
type1_error_dir <- paste0(sim_dir, 'type1error/')
out_dir <- paste0(type1_error_dir, scenario, '_', sim_i, suffix,'/')
out_dir_train <- paste0(out_dir, 'train/sims/')

## make scenario out_dirs
mkdir <- function(path){ ifelse(file.exists(path), '', dir.create(path))}

mkdir(out_dir_train)
mkdir(paste0(out_dir, 'pred/split_files/'))
mkdir(paste0(out_dir, 'pred/target_files/'))
mkdir(paste0(out_dir, 'pred/logs/'))
mkdir(paste0(out_dir, 'results/'))

#################
## do sims

set.seed(1233218)
for (i in 1:10^6) {
	gtex_out <- cbind(gtex_gt, 'TargetID'=i, 'ES'=resample(gtex_weights))
	rosmap_out <- cbind(rosmap_gt, 'TargetID'=i, 'ES'=resample(rosmap_weights))
	rosmap_v_out <- cbind(rosmap_gt, 'TargetID'=i, 'ES'=resample(rosmap_v_weights))

	gtex_out_path <- paste0(out_dir_train, 'PrediXcan-GTEx_train_weight_', i, '.txt')
	rosmap_out_path <- paste0(out_dir_train, 'TIGAR-ROSMAP_train_weight_', i, '.txt')
	rosmap_v_out_path <- paste0(out_dir_train, 'TIGAR-ROSMAP_valid_train_weight_', i, '.txt')

	# output
	write.table(
		gtex_out,
		gtex_out_path,
		quote=FALSE,
		append=FALSE,
		col.names=TRUE,
		row.names=FALSE,
		sep='\t')

	write.table(
		rosmap_out,
		rosmap_out_path,
		quote=FALSE,
		append=FALSE,
		col.names=TRUE,
		row.names=FALSE,
		sep='\t')

	write.table(
		rosmap_v_out,
		rosmap_v_out_path,
		quote=FALSE,
		append=FALSE,
		col.names=TRUE,
		row.names=FALSE,
		sep='\t')

}
#


#################
## output expression file
search_targetID <- paste0(scenario, '_', sim_i)
expr_path <- paste0(sim_dir, 'expression/ROSMAP_expr_overlap.txt'

# get row number
row_match <- grep(search_targetID, readLines(expr_path))

# get expr header
expr_header <- colnames(read.table(expr_path, header=TRUE, check.names=FALSE, nrows=1))

# get correct row, set header
expr <- setNames(read.table(text=scan(expr_path, '', skip=row_match - 1, nlines=1, sep='\n'), header=FALSE, sep='\t'),
	expr_header)

write.table(
	expr,
	paste0(out_dir, 'ROSMAP_expr.txt'),
	quote=FALSE,
	append=FALSE,
	col.names=TRUE,
	row.names=FALSE,
	sep='\t')


#################
## output annot file for doing the SR/Naive training
annot <- data.frame('TargetID'=1:10^6)
write.table(
	annot,
	paste0(type1_error_dir, 'TargetIDs.txt'),
	quote=FALSE,
	append=FALSE,
	col.names=TRUE,
	row.names=FALSE,
	sep='\t')


