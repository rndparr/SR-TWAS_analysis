#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=F, digits=20)

# read in passed arguments
args <- (commandArgs(TRUE))
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
	dataset = args[[1]]
	n_train = as.numeric(args[[2]])
	n_valid = as.numeric(args[[3]])
	n_test = as.numeric(args[[4]])
}
# # example input:
# n_train <- 465
# dataset <- 'ROSMAP'; n_valid <- 400; n_test <- 800
# dataset <- 'GTEx'; n_valid <- 0; n_test <- 0


######################
# function to return sorted samples
get_sample <- function(samples, n_select, remove_set=NULL){
	selected <- sample(setdiff(samples, remove_set), n_select)
	return(samples[samples %in% selected])
}

######################
# set directory
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/'

# input file
sample_path <- paste0(dir, dataset, '_sampleid.txt')

# read in sampleids (previously obtained from genotype data file)
sampleid <- read.table(sample_path, header=FALSE)$V1

# output files
sample_train_path <- paste0(dir, dataset, '_train_',  n_train, '_sampleid.txt')
sample_valid_path <- paste0(dir, dataset, '_valid_',  n_valid, '_sampleid.txt')
sample_test_path <- paste0(dir, dataset, '_test_',  n_test, '_sampleid.txt')

# randomly select subjects for training/validation/testing
set.seed(1234567)

# train sampleids
sample_train <- get_sample(sampleid, n_train)

write.table(sample_train, 
	sample_train_path, 
	quote = FALSE, 
	row.names = FALSE, 
	col.names = FALSE)

# validation sampleids
if(n_valid > 0){
	sample_valid <- get_sample(sampleid, n_valid, sample_train)

	write.table(sample_valid, 
		sample_valid_path, 
		quote = FALSE, 
		row.names = FALSE, 
		col.names = FALSE)	

} else {
	sample_valid <- NULL
}

# test sampleids
if (n_test > 0){
	sample_test <- get_sample(sampleid, n_test, c(sample_train, sample_valid))

	write.table(sample_test, 
		sample_test_path, 
		quote = FALSE, 
		row.names = FALSE, 
		col.names = FALSE)
}

