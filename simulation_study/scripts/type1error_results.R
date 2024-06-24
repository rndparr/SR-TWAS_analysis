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
	scenario_sim_i_suffix = '0.1_0.1_0.5_0.1_0.1_667'
} else {
	sim_dir = as.character(args[[1]])
	scenario_sim_i_suffix = as.character(args[[2]])
}


# load xtable
library(xtable)

######################################################
# directory for simulation files
wkdir <- paste0(sim_dir, '/type1_error/', scenario_sim_i_suffix, '/')

pred_out_cols <- c('i', 'TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg')
pred_out_col_classes <- c('character', rep('numeric', 6))

# read in data
dat <- read.table(paste0(wkdir, 'results/pred_sig_results.txt'), header=FALSE, colClasses=pred_out_col_classes)
colnames(dat) <- pred_out_cols

# function to get results for table
get_results <- function(sig_lvl) {
	colSums(dat[, -1] < sig_lvl, na.rm=TRUE)
}

# significance levels to use in table
sig_lvls <- setNames(c(10^-4, 10^-5, 2.5*10^-6, 10^-6), c('1e-4', '1e-5', '2.5e-6', '1e-6'))

# get table
dat_table <- t(sapply(sig_lvls, get_results)) / 10^6

write.table(
	dat_table,
	paste0(wkdir, 'results/results_table.txt'),
	quote=FALSE,
	row.names=TRUE,
	col.names=TRUE,
	sep='\t')

# latex table
print(xtable(dat_table, display=c('g',rep('g',6))), math.style.exponents=TRUE)



