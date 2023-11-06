#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
# .libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE, digits=20)

# load libraries
library(doFuture)
library(foreach)

# load grguments
args=(commandArgs(TRUE))
cat(paste0('args:\n'))
print(args)
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==4) {
	sim_dir = as.character(args[[1]])
	jobs = as.character(args[[2]])
	suffix = ''
	ncores = as.integer(args[[3]])
	N_sims = as.integer(args[[4]])
} else {
	sim_dir = as.character(args[[1]])
	jobs = as.character(args[[2]])
	suffix = as.character(args[[3]])
	ncores = as.integer(args[[4]])
	N_sims = as.integer(args[[5]])
}

# set up parallel environment
options(future.rng.onMisuse='ignore')
registerDoFuture()
if (ncores > 1){
	pln <- plan(multicore, workers=ncores)
	cat(paste0('Starting job with ', ncores, ' cores\n'))
} else {
	plan(sequential)
}

jobs <- strsplit(jobs, ',')[[1]]

## job format: 
# [rosmap_prop_causal]_[gtex_prop_causal]_[gtex_prop_causal_overlap]_[rosmap_he2]_[gtex_he2]

# # jobs should be comma separated list
# # suffix give unique name to jobs; could be date, batch, etc
# suffix <- '_batch1'

######################
## MAIN PAPER
# # half overlap snps, same heritability
# jobs <- '0.001_0.001_0.5_0.1_0.1,0.001_0.001_0.5_0.2_0.2,0.001_0.001_0.5_0.5_0.5,0.01_0.01_0.5_0.1_0.1,0.01_0.01_0.5_0.2_0.2,0.01_0.01_0.5_0.5_0.5,0.05_0.05_0.5_0.1_0.1,0.05_0.05_0.5_0.2_0.2,0.05_0.05_0.5_0.5_0.5,0.1_0.1_0.5_0.1_0.1,0.1_0.1_0.5_0.2_0.2,0.1_0.1_0.5_0.5_0.5'
# suffix <- ''

### SUPPLEMENT JOBS
# # same snps, same heritability
# jobs <- '0.001_0.001_1_0.1_0.1,0.001_0.001_1_0.2_0.2,0.001_0.001_1_0.5_0.5,0.01_0.01_1_0.1_0.1,0.01_0.01_1_0.2_0.2,0.01_0.01_1_0.5_0.5,0.05_0.05_1_0.1_0.1,0.05_0.05_1_0.2_0.2,0.05_0.05_1_0.5_0.5,0.1_0.1_1_0.1_0.1,0.1_0.1_1_0.2_0.2,0.1_0.1_1_0.5_0.5'
# suffix <- '_sameSNPsameHe'

# # same snps, half heritability
# jobs <- '0.001_0.001_1_0.1_0.05,0.001_0.001_1_0.2_0.1,0.001_0.001_1_0.5_0.25,0.01_0.01_1_0.1_0.05,0.01_0.01_1_0.2_0.1,0.01_0.01_1_0.5_0.25,0.05_0.05_1_0.1_0.05,0.05_0.05_1_0.2_0.1,0.05_0.05_1_0.5_0.25,0.1_0.1_1_0.1_0.05,0.1_0.1_1_0.2_0.1,0.1_0.1_1_0.5_0.25'
# suffix <- '_sameSNPhalfHe'
######################


######################
## FUNCTIONS
######################
# calculate/impute MAF
calc_maf <- function(x){
	maf <- sum(x, na.rm=TRUE) / (2 * (length(x) - sum(is.nan(x))))
	return(maf)
}

impute_maf <- function(df){
	if (any(is.na(df))) {
		df[,-1] <- t(apply(df[,-1], 1, function(x) ifelse(is.na(x), calc_maf(x), x)))
	}
	return(df)
}

# function to get sampleids for GTEx, ROSMAP from header
get_sampleids <- function(dataset){
	path <- paste0(sim_dir, 'genotype/', dataset, '_ABCA7_raw.dosage.gz')
	sampleids <- colnames(read.delim(con <- gzfile(path, 'r'), 
			check.names=FALSE, nrow=1))[-c(1:5)]; close(con)
	return(sampleids)
}

# read in genotype data
get_gt <- function(dataset){
	path <- paste0(sim_dir, 'genotype/', dataset, '_ABCA7_raw.dosage.gz')
	gt <- read.delim(con <- gzfile(path, 'r'), check.names=FALSE); close(con)
	return(gt)
}

# minor rounding issues can mean n=length(x)+1 or n=-1
safe_sample <- function(x, n){ 
	sample(x, min(length(x), max(0, n)))
}

# set sampleids list
sampleid <- list(
	'GTEx' = get_sampleids('GTEx'),
	'ROSMAP' = get_sampleids('ROSMAP')
	)

# set genotype list
gt <- list(
	'GTEx' = get_gt('GTEx'),
	'ROSMAP' = get_gt('ROSMAP')
	)

### BASED ON ANY OF THE TOTAL SNPS
both <- merge(gt[['GTEx']][, c('snpID', sampleid[['GTEx']])],
	gt[['ROSMAP']][,c('snpID', sampleid[['ROSMAP']])],
	by='snpID', all=TRUE)

out_dir <- paste0(sim_dir, 'expression/')

# cat('Expression files in this directory calculated from rosmap_pc proportion of the 5924 total SNPs found in data; causal SNPs may be missing from one of the datasets; if SNP not found in dataset, MAF imputation was used in generating the file.', file=paste0(out_dir, 'note.txt'), sep='\n')

### OUTPUT HEADERS OF THE ALL EXPRESSION FILES
info_out_cols <- c('CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName')

cat(paste(c(info_out_cols, sampleid[['GTEx']], sampleid[['ROSMAP']]), 
		collapse='\t', sep=''), 
	file=paste0(out_dir, 'both_expr', suffix, '.txt'), sep='\n')

cat(paste(c(info_out_cols, sampleid[['GTEx']]), 
		collapse='\t', sep=''), 
	file=paste0(out_dir, 'GTEx_expr', suffix, '.txt'), sep='\n')

cat(paste(c(info_out_cols, sampleid[['ROSMAP']]), 
		collapse='\t', sep=''), 
	file=paste0(out_dir, 'ROSMAP_expr', suffix, '.txt'), sep='\n')

# get gene expression
get_expr <- function(expr, he2) {
	if (he2 == 0) {
		expr <- rnorm(length(expr), mean = 0, sd = 1)
		error_term <- rnorm(length(expr), mean = 0, sd = 1)
		expr_raw <- t(expr + error_term)
	} else {
		# gamma^2 x var(X*B) = X*b * heritability
		gamma <- as.vector(sqrt(he2 / var(expr)))
		expr_1 <- gamma * expr
		# Error term = N(0, 1 - he^2)
		error_term <- rnorm(length(expr_1), mean = 0, 
			sd = sqrt(1 - he2))
		expr_raw <- t(expr_1 + error_term)
	}
	return(expr_raw)
}

######################
## DO SIMULATION
######################

set.seed(654321)

for (job in jobs){

	# get job params; then find causal snps to use
	job_list <- as.numeric(strsplit(job,'_')[[1]])
	
	rosmap_pc <- job_list[1]
	gtex_pc <- job_list[2]
	gtex_poverlap <- job_list[3]
	rosmap_he2 <- job_list[4]
	gtex_he2 <- job_list[5]

	# get snps
	rosmap_snps <- gt[['ROSMAP']][,'snpID']
	gtex_snps <- gt[['GTEx']][,'snpID']
	overlap_snps <- intersect(rosmap_snps, gtex_snps)

	## get number of snps in rosmap, number of causal snps based on rosmap_pc
	n_snp_rosmap <- length(rosmap_snps)
	n_snp_gtex <- length(gtex_snps)

	## number of true causal to get
	n_rosmap_true_causal <- ceiling(rosmap_pc * n_snp_rosmap)
	n_gtex_true_causal <- ceiling(gtex_pc * n_snp_gtex)

	## number of snps in overlap calculation
	n_gtex_true_causal_overlap <- floor(n_gtex_true_causal * gtex_poverlap)
	n_gtex_true_causal_other <- n_gtex_true_causal - n_gtex_true_causal_overlap

	n_rosmap_true_causal_overlap <- safe_sample(n_gtex_true_causal_overlap:n_rosmap_true_causal, 1)
	n_rosmap_true_causal_other <- n_rosmap_true_causal - n_rosmap_true_causal_overlap

	## get rosmap true causal snps
	rosmap_true_causal_overlap <- safe_sample(overlap_snps, n_rosmap_true_causal_overlap)
	rosmap_true_causal_other <- safe_sample(setdiff(rosmap_snps, rosmap_true_causal_overlap), n_rosmap_true_causal_other)
	rosmap_true_causal <- sort(c(rosmap_true_causal_overlap, rosmap_true_causal_other))

	## get gtex true causal snps
	gtex_true_causal_overlap <- safe_sample(intersect(rosmap_true_causal, gtex_snps), n_gtex_true_causal_overlap)
	gtex_true_causal_other <- safe_sample(setdiff(gtex_snps, gtex_true_causal_overlap), n_gtex_true_causal_other)
	gtex_true_causal <- sort(c(gtex_true_causal_overlap, gtex_true_causal_other))

	## overlap true causal
	overlap_true_causal <- sort(intersect(rosmap_true_causal, gtex_true_causal))
	rosmap_only_true_causal <- sort(setdiff(rosmap_true_causal, overlap_true_causal))
	gtex_only_true_causal <- sort(setdiff(gtex_true_causal, overlap_true_causal))

	# size
	n_overlap_true_causal <- length(overlap_true_causal)
	n_rosmap_only_true_causal <- length(rosmap_only_true_causal)
	n_gtex_only_true_causal <- length(gtex_only_true_causal)

	# IMPUTE MAF
	both_imputed <- impute_maf(both[both$snpID %in% union(union(rosmap_true_causal, gtex_true_causal), overlap_true_causal), ])
	rownames(both_imputed) <- both_imputed$snpID

	# separate matrices; only overlap_matrix guaranteed to have anything in it
	overlap_snp_mat <- t(as.matrix(both_imputed[both_imputed$snpID %in% overlap_true_causal, -1]))
	rosmap_only_snp_mat <- t(as.matrix(both_imputed[both_imputed$snpID %in% rosmap_only_true_causal, -1 ]))[sampleid[['ROSMAP']],]
	gtex_only_snp_mat <- t(as.matrix(both_imputed[both_imputed$snpID %in% gtex_only_true_causal, -1]))[sampleid[['GTEx']],]

	######################
	## SIMULATE EXPRESSION FOR SNPS N_sims TIMES
	######################
	# start time
	start_time <- Sys.time()

	# do N_sims simulations/different expression
	expr_sims <- foreach(i=1:N_sims, .combine=rbind) %dopar% {

		# simulate effect size \beta from N(0,1)
		beta_overlap <- as.matrix(rnorm(n_overlap_true_causal, mean=0, sd=1))
		beta_rosmap_only <- as.matrix(rnorm(n_rosmap_only_true_causal, mean=0, sd=1))
		beta_gtex_only <- as.matrix(rnorm(n_gtex_only_true_causal, mean=0, sd=1))

		# expr = X*\beta
		expr_overlap_0 <- overlap_snp_mat %*% beta_overlap
		expr_rosmap_only_0 <- rosmap_only_snp_mat %*% beta_rosmap_only
		expr_gtex_only_0 <- gtex_only_snp_mat %*% beta_gtex_only

		# combine overlap, and cohort_only expr_0
		expr_rosmap_0 <- expr_overlap_0[sampleid[['ROSMAP']],] + expr_rosmap_only_0
		expr_gtex_0 <- expr_overlap_0[sampleid[['GTEx']],] + expr_gtex_only_0

		# account for heritability, error
		rosmap_expr <- get_expr(expr_rosmap_0, rosmap_he2)
		gtex_expr <- get_expr(expr_gtex_0, gtex_he2)
		
		expr_raw <- cbind(rosmap_expr, gtex_expr)

		target_id <- paste(job, as.character(i), sep='_')

		expr <- data.frame('TargetID'=target_id, expr_raw, check.names=FALSE)

	} %seed% 6543217

	######################
	## OUTPUT
	######################
	# set up gene info dataframe
	info_df <- data.frame('CHROM'=19, 'GeneStart'=1040101, 'GeneEnd'=1065572, 'GeneName'='ABCA7')
	info_out_cols <- c('CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName')

	# combine info_df and expr_sim, put columns in order
	out_cols <- c(info_out_cols, sampleid[['GTEx']], sampleid[['ROSMAP']])
	out_expr <- data.frame(info_df, expr_sims, check.names=FALSE)[, out_cols]

	# output to file
	all_out_path <- paste0(out_dir, 'both_expr', suffix, '.txt')

	## append to all expression file
	write.table(
		out_expr,
		all_out_path,
		quote=FALSE,
		append=TRUE,
		row.names=FALSE,
		col.names=FALSE,
		sep='\t')

	# output to individual dataset files
	for (dataset in c('GTEx', 'ROSMAP')){
		# combine info_df and expr_sim, put columns in order
		out_cols <- c(info_out_cols, sampleid[[dataset]])
		out_expr <- data.frame(info_df, expr_sims, check.names=FALSE)[, out_cols]

		out_path <- paste0(out_dir, dataset, '_expr', suffix, '.txt')
		write.table(
			out_expr,
			out_path,
			quote=FALSE,
			append=TRUE,
			row.names=FALSE,
			col.names=FALSE,
			sep='\t')
	}

	# time elapsed
	end_time <- Sys.time()
	cat(paste('Computation time:', end_time - start_time, '\n'))

}
