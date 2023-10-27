options(stringsAsFactors=FALSE)


#################
## functions
read_bgzip <- function(path) {
	raw <- read.delim(con <- gzfile(path, 'r'), check.names = FALSE); close(con)
	return(raw)
}

resample <- function(x, ...) x[sample.int(length(x), ...)]

#################

## set scenario
scenario <- 'overlap_0.1_0.1_0.5_0.1_0.1'
# gtex_model <- ''
gtex_model <- '_EN'

## set seed, pick a simulation
set.seed(1234567)
sim_i <- sample(1:1000, 1)
# 508

#################
## get trained weights
weight_dir <- '/home/rparrish/YangFSSdata/SR_TWAS/sim/train/sims/'

gtex_weight_path <- paste0(weight_dir, 'GTEx_train', gtex_model, '_weight_', scenario, '_', sim_i, '.txt.gz')
rosmap_weight_path <- paste0(weight_dir, 'ROSMAP_train_weight_', scenario, '_', sim_i, '.txt.gz')

gtex_weights_raw <- read_bgzip(gtex_weight_path)
rosmap_weights_raw <- read_bgzip(rosmap_weight_path)

gtex_weights_nsnps <- nrow(gtex_weights_raw)
rosmap_weights_nsnps <- nrow(rosmap_weights_raw)


## get genotype file rows
gt_dir <- '/home/rparrish/YangFSSdata/SR_TWAS/sim/genotype/'

gtex_gt_path <- paste0(gt_dir, 'GTEx_ABCA7_raw.dosage.gz')
rosmap_gt_path <- paste0(gt_dir, 'ROSMAP_ABCA7_raw.dosage.gz')

gtex_gt <- read_bgzip(gtex_gt_path)[, 1:5]
rosmap_gt <- read_bgzip(rosmap_gt_path)[, 1:5]

gtex_gt_nsnps <- nrow(gtex_gt)
rosmap_gt_nsnps <- nrow(rosmap_gt)


## get final weight vector with same length as available snps
gtex_weights <- c(gtex_weights_raw[, 'ES'], rep(0, gtex_gt_nsnps - gtex_weights_nsnps))
rosmap_weights <- c(rosmap_weights_raw[, 'ES'], rep(0, rosmap_gt_nsnps - rosmap_weights_nsnps))


#################
### output directory
type1_error_dir <- '/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/'
out_dir <- paste0(type1_error_dir, scenario, '_', sim_i, '/')
out_dir_train <- paste0(out_dir, 'train/sims/')


#################
## do sims


# for (i in 1:10^6) {
for (i in 1:924296) {
	gtex_out <- cbind(gtex_gt, 'TargetID'=i, 'ES'=resample(gtex_weights))
	rosmap_out <- cbind(rosmap_gt, 'TargetID'=i, 'ES'=resample(rosmap_weights))

	gtex_out_path <- paste0(out_dir_train, 'GTEx', gtex_model, '_weights_', i, '.txt')
	rosmap_out_path <- paste0(out_dir_train, 'ROSMAP_weights_', i, '.txt')

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

}
#

set.seed(7654321)

for (i in 924297:10^6) {

	gtex_out <- cbind(gtex_gt, 'TargetID'=i, 'ES'=resample(gtex_weights))
	rosmap_out <- cbind(rosmap_gt, 'TargetID'=i, 'ES'=resample(rosmap_weights))

	gtex_out_path <- paste0(out_dir_train, 'GTEx', gtex_model, '_weights_', i, '.txt')
	rosmap_out_path <- paste0(out_dir_train, 'ROSMAP_weights_', i, '.txt')

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

}



#################
## output expression file
search_targetID <- '0.1_0.1_0.5_0.1_0.1_508'
expr_path <- '/home/rparrish/YangFSSdata/SR_TWAS/sim/expression/raw_dosage/ROSMAP_expr_overlap.txt'

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



##################################
# n_skip <- row_match - 1
# cbind(scan(expr_path, '', skip=n_skip, nlines=1, sep='\n'))

# gtex_out
# rosmap_out
# gtex_out_path <- paste0(out_dir, 'GTEx', gtex_model, '_weights_', i, '.txt')
# rosmap_out_path <- paste0(out_dir, 'ROSMAP_weights_', i, '.txt')
# # paste0(out_dir, 'GTEx', gtex_model, '_train_weight_', scenario, '_', sim_i, '_', i, '.txt')
# gtex_model
# 'GTEx_EN_train_weight_'
# bleh <- cbind(gtex_gt, 'ES'=resample(gtex_weights))
# resample(gtex_weights)

