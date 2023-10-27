#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=F, digits=5)


############################################
## UPDATE JOB NAMING
newjobvec <- function(x){
	if(length(x) == 3){
		# (causal_prop for both)_(GTEx_He2)_(ROSMAP_He2) ->
		# (ROSMAP_causal_prop)_(GTEx_causal_prop)_(GTEx_overlap)_(ROSMAP_He2)_(GTEx_He2)
		return(c(x[[1]], x[[1]], '1', x[[3]], x[[2]]))
	} else {
		return(x)
	}
}
newjobvec2 <- function(x){
	if(length(x) == 4){
		# (causal_prop for both)_(GTEx_He2)_(ROSMAP_He2) ->
		# (ROSMAP_causal_prop)_(GTEx_causal_prop)_(GTEx_overlap)_(ROSMAP_He2)_(GTEx_He2)
		return(c(x[[1]], x[[1]], '1', x[[3]], x[[2]], x[[4]]))
	} else {
		return(x)
	}
}

oldjobstr <- function(job){
	x <- strsplit(job, '_')[[1]]
	if(length(x) == 3){
		# (causal_prop for both)_(GTEx_He2)_(ROSMAP_He2) ->
		# (ROSMAP_causal_prop)_(GTEx_causal_prop)_(GTEx_overlap)_(ROSMAP_He2)_(GTEx_He2)
		return(job)
	} else {
		return(paste(x[[1]], x[[5]], x[[4]], sep='_'))
	}
}


############################################
pheno_h2_list <- c(0.05, 0.1, 0.175, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)
oldjobs <- c('0.001_0.1_0.1','0.001_0.2_0.2','0.001_0.5_0.5','0.01_0.1_0.1','0.01_0.2_0.2','0.01_0.5_0.5','0.05_0.1_0.1','0.05_0.2_0.2','0.05_0.5_0.5','0.1_0.1_0.1','0.1_0.2_0.2','0.1_0.5_0.5','0.001_0.05_0.1','0.001_0.1_0.2','0.001_0.25_0.5','0.001_0_0.1','0.001_0_0.2','0.001_0_0.5','0.01_0.05_0.1','0.01_0.1_0.2','0.01_0.25_0.5','0.01_0_0.1','0.01_0_0.2','0.01_0_0.5','0.05_0.05_0.1','0.05_0.1_0.2','0.05_0.25_0.5','0.05_0_0.1','0.05_0_0.2','0.05_0_0.5','0.1_0.05_0.1','0.1_0.1_0.2','0.1_0.25_0.5','0.1_0_0.1','0.1_0_0.2','0.1_0_0.5'
	,
	'0.1_0.1_0.5_0.2_0.2','0.1_0.1_0.5_0.5_0.5',
	'0.001_0.001_0.5_0.1_0.1', '0.001_0.001_0.5_0.2_0.2', '0.001_0.001_0.5_0.5_0.5', 
	'0.01_0.01_0.5_0.1_0.1', '0.01_0.01_0.5_0.2_0.2', '0.01_0.01_0.5_0.5_0.5', 
	'0.05_0.05_0.5_0.1_0.1', '0.05_0.05_0.5_0.2_0.2', '0.05_0.05_0.5_0.5_0.5', 
	'0.1_0.1_0.5_0.1_0.1'
	)

# jobs_info <- as.data.frame(do.call(rbind, strsplit(jobs, '_')))
jobs_info <- as.data.frame(do.call(rbind, lapply(strsplit(oldjobs, '_'), newjobvec)))
# colnames(jobs_info) <- c('causal_prop','GTEx_He2','ROSMAP_He2')
colnames(jobs_info) <- c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap','ROSMAP_He2','GTEx_He2')
jobs_info <- as.data.frame(sapply(jobs_info, as.numeric))
jobs_info$oldjob <- oldjobs
jobs_info$job <- with(jobs_info, paste(ROSMAP_causal_prop, GTEx_causal_prop,  GTEx_overlap, ROSMAP_He2, GTEx_He2, sep='_'))
# jobs_info <- jobs_info[,  c('oldjob','causal_prop','GTEx_He2','ROSMAP_He2')]
jobs_info$job_cat <- ''
jobs_info[jobs_info$GTEx_He2==jobs_info$ROSMAP_He2, 'job_cat'] <- 'same'
jobs_info[2*jobs_info$GTEx_He2==jobs_info$ROSMAP_He2, 'job_cat'] <- 'half'
jobs_info[jobs_info$GTEx_He2==0, 'job_cat'] <- 'zero'
jobs_info[jobs_info$GTEx_overlap!=1, 'job_cat'] <- 'overlap'

# same_he2 <- jobs_info[jobs_info$GTEx_He2==jobs_info$ROSMAP_He2, 'job']
# half_he2 <- jobs_info[2*jobs_info$GTEx_He2==jobs_info$ROSMAP_He2, 'job']
# zero_he2 <- jobs_info[jobs_info$GTEx_He2==0, 'job']
# pheno_h2_list <- c(0.5, 0.75)

############################################
# set directories
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/'
out_dir <- paste0(dir, 'plot/')
setwd(out_dir)

############################################
# set up data frames
# training info data
sr_dat <- data.frame()
train_dat <- data.frame()


train_info_paths <- c(
	paste0(dir, 'train/', 'ROSMAP_train_info.txt'),
	paste0(dir, 'train/', 'GTEx_train_EN_info.txt'),
	paste0(dir, 'train/', 'Naive_train_info_EN.txt'),
	paste0(dir, 'train/', 'SR_train_info_EN.txt'),
	paste0(dir, 'train/', 'GTEx_train_EN_info_123121.txt'),
	paste0(dir, 'train/', 'Naive_train_info_EN_123121.txt'),
	paste0(dir, 'train/', 'SR_train_info_EN_123121.txt'),
	paste0(dir, 'train/', 'ROSMAP_train_info_120621.txt'),
	paste0(dir, 'train/', 'GTEx_train_EN_info_120621.txt'),
	paste0(dir, 'train/', 'Naive_train_info_EN_120621.txt'),
	paste0(dir, 'train/', 'SR_train_info_EN_120621.txt'),
	paste0(dir, 'train/', 'ROSMAP_train_info_121721.txt'),
	paste0(dir, 'train/', 'GTEx_train_EN_info_121721.txt'),
	paste0(dir, 'train/', 'Naive_train_info_EN_121721.txt'),
	paste0(dir, 'train/', 'SR_train_info_EN_121721.txt'),
	paste0(dir, 'train/', 'GTEx_train_EN_info_overlap.txt'),
	paste0(dir, 'train/', 'Naive_train_info_EN_overlap.txt'),
	paste0(dir, 'train/', 'SR_train_info_EN_overlap.txt')
)

# 'ROSMAP','GTEx','Naive','SR','ROSMAP','GTEx','Naive','SR'

datasets <- c('ROSMAP', 'GTEx', 'Naive', 'SR', 'GTEx', 'Naive', 'SR','ROSMAP','GTEx','Naive','SR','ROSMAP','GTEx','Naive','SR','GTEx','Naive','SR')
n_j <- length(datasets)

for (j in 1:n_j) {
	# for (dataset in c('ROSMAP', 'GTEx', 'Naive', 'SR')) {

	train_dat_path <- train_info_paths[[j]]
	dataset <- datasets[[j]]

	# train_dat_path <- paste0(dir, 'train/', dataset, '_train_info', suffix, '.txt')
	# train_dat_path <- paste0(dir, 'train/', dataset, '_train_info.txt')

	if (dataset == 'ROSMAP') {
		col_classes <- rep('NULL', 12)
		col_classes[c(9,11)] <- rep('numeric', 2)
	} else if (dataset == 'GTEx') {
		col_classes <- rep('NULL', 16)
		col_classes[c(9,11)] <- rep('numeric', 2)
	} else if(dataset == 'Naive') {
		col_classes <- rep('NULL', 19)
		col_classes[c(9:10, 13:14, 17:18)] <- rep('numeric', 6)
	} else if(dataset == 'SR') {
		col_classes <- rep('NULL', 21)
		col_classes[c(9:10, 12:13, 15:16, 19:20)] <- rep('numeric', 8)
	}
	col_classes[4] <- 'character'

	train_dat_part <- read.table(train_dat_path, header=TRUE, check.names=FALSE, sep='\t', colClasses=col_classes, skipNul=TRUE)
	# train_dat_part <- read.table(train_dat_path, header=TRUE, check.names=FALSE, sep='\t', colClasses=col_classes, skipNul=TRUE,nrows=5)
	# bleh <- read.table(train_dat_path, header=TRUE, check.names=FALSE, sep='\t', nrow=5)

	if(dataset == 'Naive' || dataset == 'SR') {
		if(dataset == 'SR') {
			sr_dat_part <- read.table(train_dat_path, header=TRUE, check.names=FALSE, sep='\t')[, c(1:5, 12,13)]
			sr_dat_part$cohort <- dataset
			sr_dat <- rbind(sr_dat_part, sr_dat)
		} else if (dataset == 'Naive') {
			train_dat_part$Z0 <- 0.5
			train_dat_part$Z1 <- 0.5
		}
		train_dat_part <- rbind(
			data.frame('TargetID'=train_dat_part$TargetID,  
				'CVR2'=NA, 'TrainR2'=NA,
				'ValidCVR2'=train_dat_part$CVR2, 
				'ValidR2'=train_dat_part$R2, 
				'cohort'=dataset,
				'Z0'=train_dat_part$Z0,
				'Z1'=train_dat_part$Z1),
			data.frame(
				'TargetID'=train_dat_part$TargetID, 
				'CVR2'=NA, 'TrainR2'=NA,
				'ValidCVR2'=train_dat_part$W0_CVR2, 
				'ValidR2'=train_dat_part$W0_R2, 
				'cohort'='GTEx',
				'Z0'=NA,
				'Z1'=NA),
			data.frame(
				'TargetID'=train_dat_part$TargetID, 
				'CVR2'=NA, 'TrainR2'=NA,
				'ValidCVR2'= train_dat_part$W1_CVR2,
				'ValidR2'=train_dat_part$W1_R2, 
				'cohort'='ROSMAP',
				'Z0'=NA,
				'Z1'=NA))
	} else {
		colnames(train_dat_part) <- c('TargetID', 'CVR2', 'TrainR2')
		train_dat_part$cohort <- dataset
		train_dat_part$ValidCVR2 <- NA
		train_dat_part$ValidR2 <- NA
		train_dat_part$Z0 <- NA
		train_dat_part$Z1 <- NA		
	}
	train_dat_part <- train_dat_part[, c('TargetID', 'CVR2', 'TrainR2', 'ValidCVR2', 'ValidR2', 'cohort', 'Z0', 'Z1')]
	train_dat <- rbind(train_dat_part, train_dat)
}


# get job info
train_dat$oldjob <- sub('_[^_]+$', '', train_dat$TargetID)
target_info <- do.call(rbind, lapply(strsplit(train_dat$TargetID, '_'), newjobvec2))
colnames(target_info) <- c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap','ROSMAP_He2', 'GTEx_He2', 'i')

train_dat$job <- paste(target_info[,1], target_info[,2], target_info[,3], target_info[,4], target_info[,5], sep='_')
train_dat <- cbind(target_info, train_dat)
cols <- intersect(colnames(train_dat), colnames(jobs_info))
train_dat <- merge(jobs_info, train_dat, all.y=TRUE, by=c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap','ROSMAP_He2', 'GTEx_He2','oldjob','job'))

train_dat$cohort <- factor(train_dat$cohort, levels=c('GTEx','ROSMAP','Naive','SR'),
	labels=c('PrediXcan_GTEx', 'TIGAR_ROSMAP', 'Naive', 'SR-TWAS'))


# train_dat$CVR2 <- 100 * train_dat$CVR2
# train_dat$TrainR2 <- 100 * train_dat$TrainR2
# train_dat$ValidCVR2 <- 100 * train_dat$ValidCVR2
# train_dat$ValidR2 <- 100 * train_dat$ValidR2

############################################
# expr pred, pheno pred, power data

pred_dat <- data.frame()

for (suffix in c('_EN')) {
	pred_dat_path <- paste0(dir, 'power/all_pred_results', suffix, '.txt')
	pred_dat_part <- read.table(pred_dat_path, header=TRUE, check.names=FALSE, sep='\t')

	# if (!('old_job' %in% colnames(pred_dat_part))) {
	# 	pred_dat_part$old_job <- pred_dat_part$job
	# }
	pred_dat_part$oldjob <- sub('_[^_]+$', '', pred_dat_part$TargetID)

	# pred_dat_part$oldjob <- pred_dat_part$old_job
	pred_dat_part <- pred_dat_part[, c('job', 'i', 'cohort', 'ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2', 'expr_R2', 'pheno_R2', 'pheno_pval', 'oldjob', 'TargetID')]
	pred_dat <- rbind(pred_dat, pred_dat_part)
}
pred_dat$cohort <- factor(pred_dat$cohort, levels=c('GTEx','ROSMAP','Naive','SR'),
	labels=c('PrediXcan_GTEx', 'TIGAR_ROSMAP', 'Naive', 'SR-TWAS'))
pred_dat <- merge(jobs_info, pred_dat, all.y=TRUE, by=c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap','ROSMAP_He2', 'GTEx_He2','job','oldjob'))




# dput(colnames(pred_dat_part))
# pred_dat_path <- paste0(dir, 'power/all_pred_results.txt')
# # pred_dat_path <- paste0(dir, 'power/all_pred_results_test.txt')
# pred_dat <- read.table(pred_dat_path, header=TRUE, check.names=FALSE, sep='\t')
# pred_dat <- merge(jobs_info, pred_dat, all.y=TRUE)



# pred_dat$expr_R2 <- 100 * pred_dat$expr_R2
# pred_dat$pheno_R2 <- 100 * pred_dat$pheno_R2

############################################
# summary power data
power_dat <- data.frame()

for (suffix in c('_EN')) {
	power_dat_path <- paste0(dir, 'power/all_power_results', suffix, '.txt')
	power_dat_part <- read.table(power_dat_path, header=TRUE, check.names=FALSE, sep='\t')

	# if (!('old_job' %in% colnames(power_dat_part))) {
	# 	power_dat_part$old_job <- power_dat_part$job
	# }
	power_dat_part$oldjob <- sapply(power_dat_part$job, oldjobstr, USE.NAMES=FALSE)

	# power_dat_part$oldjob <- power_dat_part$old_job
	power_dat_part <- power_dat_part[, c('job', 'cohort', 'ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2', 'oldjob', 'power')]
	power_dat <- rbind(power_dat, power_dat_part)
}
power_dat$cohort <- factor(power_dat$cohort, levels=c('GTEx','ROSMAP','Naive','SR'),
	labels=c('PrediXcan_GTEx', 'TIGAR_ROSMAP', 'Naive', 'SR-TWAS'))
power_dat <- merge(jobs_info, power_dat, all.y=TRUE, by=c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap','ROSMAP_He2', 'GTEx_He2','job','oldjob'))
power_dat[power_dat$GTEx_overlap!=1, 'job_cat'] <- 'overlap'

# pred_dat[pred_dat$GTEx_overlap != 1,]
# power_dat[power_dat$GTEx_overlap != 1,]
# power_dat_path <- paste0(dir, 'power/', 'all_power_results.txt')
# # power_dat_path <- paste0(dir, 'power/', 'all_power_results_test.txt')
# power_dat <- read.table(power_dat_path, header=TRUE, check.names=FALSE, sep='\t')

# power_dat$cohort <- factor(power_dat$cohort, levels=c('GTEx','ROSMAP','Naive','SR'),
# 	labels=c('TIGAR_GTEx', 'TIGAR_ROSMAP', 'Naive', 'SR-TWAS'))

# power_dat <- merge(jobs_info, power_dat, all.y=TRUE)

############################################
# save plot data
save(jobs_info, train_dat, pred_dat, power_dat, file=paste0(dir, 'plot/plot_data_EN.Rdata'))



