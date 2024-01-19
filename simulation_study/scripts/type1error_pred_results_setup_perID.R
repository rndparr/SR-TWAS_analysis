#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=FALSE)

# load arguments
args=(commandArgs(TRUE))
cat(paste0('\nargs:\n'))
print(args)
cat()
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==2) {
	sim_dir = as.character(args[[1]])
	i = as.integer(args[[2]])
	scenario_sim_i_suffix = '0.1_0.1_0.5_0.1_0.1_667'
} else {
	sim_dir = as.character(args[[1]])
	i = as.integer(args[[2]])
	scenario_sim_i_suffix = as.character(args[[3]])
}

# load libraries
library(reshape2)
library(doFuture)
library(foreach)
library(iterators)

# safely set up parallel environment
# some parallel computing packages in R do not account for shared computing setting and default to taking some percentage of cores *on the machine* rather than cores alloted to the job
# using availableCores with the method set to the correct job scheduler prevents this and removes need to manually pass the number of cores to the job 
# may need to change method to a different job scheduler
registerDoFuture()
pln <- plan(multicore, workers=availableCores(methods='SGE'))
print(paste0('Starting job with ', availableCores(methods='SGE'), ' cores'))

######################################################
# directory for simulation files
out_dir <- paste0(sim_dir, '/type1_error/', scenario_sim_i_suffix, '/')

datasets <- c('TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg')

######################################################
# sampleid
test_sampleids_path <- paste0(sim_dir, 'sampleid/ROSMAP_test_800_sampleid.txt')

# read in test sampleid data
test_sampleids <- read.table(test_sampleids_path, header = FALSE)$V1

##########
## SIM_EXPR
# path to simulated expression
sim_expr_path <- paste0(out_dir, 'ROSMAP_expr.txt')

# read in number of columns, set column classes
sim_expr_ncols <- length(colnames(read.table(sim_expr_path, header = TRUE, check.names = FALSE, nrow=1)))
sim_expr_colClasses <- c(rep('NULL', 3), 'character', 'NULL', rep('numeric', sim_expr_ncols - 5))

# read-in
sim_expr <- read.table(sim_expr_path, header=TRUE, check.names=FALSE, colClasses=sim_expr_colClasses, skipNul=TRUE)
sampleid <- colnames(sim_expr)[-1]

# center expression data
sim_expr[,sampleid] <- sim_expr[,sampleid] - rowMeans(sim_expr[,sampleid], na.rm=TRUE)

# reduce to test sampleids and reshape
sim_expr <- sim_expr[, test_sampleids]
sim_expr <- melt(sim_expr, measure.vars=test_sampleids, variable.name='sample_id', value.name='sim_expr')

# # get variance ## NEEDS TO BE FOR ALL SAMPLES NOT JUST TEST_SAMPLEIDS
# expr_var <- var(sim_expr$sim_expr)

######################################################
pred_split_dir <- paste0(out_dir, 'pred/split_files/')
target_files_dir <- paste0(out_dir, 'pred/target_files/')
# mkdir <- function(path){ ifelse(file.exists(path), '', dir.create(path)) }
# mkdir(target_files_dir)

targfil_out_cols <- c('TargetID', 'sample_id', 'sim_expr', datasets)

# print header to its own file
if (i == 1) {
	cat(paste(targfil_out_cols[-1], collapse='\t', sep=''), file=paste0(target_files_dir, 'header.txt'), sep='\n')
}

# read in pred data for a dataset,test_sampleid
read_pred <- function(dataset, pred_split_dir, test_sampleid){
	pred_path <- paste0(pred_split_dir, dataset, '_', test_sampleid, '.txt')
	pred_dataset <- read.table(pred_path, header=FALSE, colClasses=c('character', 'numeric'), fill=TRUE)
	colnames(pred_dataset) <- c('TargetID', dataset)
	return(pred_dataset)
}

# output to per-target files
output_pred <- function(row) {
	cat(paste(row[targfil_out_cols[-1]], collapse='\t', sep=''), file=paste0(target_files_dir, row[['TargetID']], '.txt'), sep='\n', append=TRUE)
}

targets <- c()

test_sampleid <- test_sampleids[i]

sim_data <- Reduce(function(...) merge(..., all=TRUE), 
	sapply(datasets, read_pred, pred_split_dir=pred_split_dir, test_sampleid=test_sampleid, simplify=FALSE))

sim_data$TargetID <- as.integer(sim_data$TargetID)
sim_data <- sim_data[order(sim_data$TargetID), ]
sim_data$sample_id <- test_sampleid
sim_data$sim_expr <- sim_expr[sim_expr$sample_id == test_sampleid, 'sim_expr']
sim_data <- sim_data[, targfil_out_cols]

targets <- sort(unique(c(sim_data$TargetID, targets)))

# output to per-target files
foreach(sim_dat=iter(sim_data, by='row'), .combine=c) %dopar% {
	output_pred(sim_dat)
}

# output list of targets
if (i == 1) {
	cat(paste(targets, sep='\n'), file=paste0(target_files_dir, 'targets.txt') )
}

