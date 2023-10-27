#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE)

######################################################
# directory for simulation files
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/'


sig_level_dat <- data.frame(
	siglvl=c(10^-4, 10^-5, 10^-6, 2.5*10^-6),
	text=c('1e-4', '1e-5', '1e-6', '2.5e-6'),
	seed=c(8654321, 6654321, 5654321, 7654321),
	path_suf=c('_1e-4', '_1e-5', '_1e-6', ''))



dat <- read.table(paste0(dir, 'results/pred_sig_results.txt'), header=FALSE, colClasses=c('character', rep('numeric',4)))
colnames(dat) <- c('i', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')

get_results <- function(sig_lvl) {
	colSums(dat[, -1] < sig_lvl, na.rm=TRUE)
}


sig_lvls <- setNames(c(10^-4, 10^-5, 2.5*10^-6, 10^-6), c('1e-4', '1e-5', '2.5e-6', '1e-6'))

bleh <- sapply(sig_lvls, get_results)

t(bleh)

##############


sum_dat <- data.frame()
for (i in 1:4){
	rawdat <- read.table(paste0(dir, 'results/pred_sig_results', sig_level_dat[i, 'path_suf'], '.txt'), header=FALSE, colClasses=c('character', rep('numeric',4)))
	colnames(rawdat) <- c('i', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')
	
	sum_dat_i <- colSums(rawdat[,-1], na.rm=TRUE)
	sum_dat_i$siglvl <- sig_level_dat[i, 'siglvl']

	sum_dat <- rbind(sum_dat, sum_dat_i)
}

sum_dat[order(-sum_dat$siglvl), c('siglvl', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')]

dat2 <- read.table(paste0(dir, 'results/pred_sig_results.txt'), header=FALSE, colClasses=c('character', rep('numeric',4)))
colnames(dat2) <- c('i', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')



# dat <- read.table(paste0(dir, 'results/pred_sig_results.txt'), header=FALSE, colClasses=c('character', rep('numeric',4)))
# colnames(dat) <- c('i', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')

# head(dat)
# # ROSMAP  GTEx_EN Naive_EN    SR_EN 
# #      1        4        3        1


# colSums(dat[,-1], na.rm=TRUE)
# # ROSMAP  GTEx_EN Naive_EN    SR_EN 
# #      1        4        3        1

# colSums(is.na(dat[,-1]))
# # ROSMAP  GTEx_EN Naive_EN    SR_EN 
# #      0        1   167791   173181
