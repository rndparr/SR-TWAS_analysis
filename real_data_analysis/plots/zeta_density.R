#!/usr/bin/env Rscript


library(ggplot2)
library(reshape)

Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE)


# twas <- 'PD_TWAS'
# twas <- 'AD_TWAS2'

if (twas == 'PD_TWAS'){
	# PDTWAS
	SRTWAS_traindir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/GTEx_6tissue/'
	
	SR_path <- function(chr){ return(paste0(SRTWAS_traindir, 'CHR', chr, '_SR_train_GTEx_GeneInfo.txt'))}

} else if (twas == 'AD_TWAS2'){
	# ADTWAS
	SRTWAS_traindir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/Brain_ROSMAP_GTEx_PredDB/'

	SR_path <- function(chr){ return(paste0(SRTWAS_traindir, 'CHR', chr, '_SR_train_4basemodels_GeneInfo.txt'))}
}

##############


# set directory
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/112323/'
suffix <- ''


## Read in data
dat <- data.frame()
for (chr in 1:22) {
	# get train info data
		# get data; grep the target id from tabix output, read as table with column names from header
		raw <- read.csv(SR_path(chr), sep='\t')
		dat <- rbind(dat, raw)
}


if (twas == 'PD_TWAS'){
	# PDTWAS
	# zcols <- paste0('Z_', c('BRNACC','BRNCDT','BRNCTXA','BRNNCC','BRNPTM','BLOOD'))
	zcols <- c('BRNACC','BRNCDT','BRNCTXA','BRNNCC','BRNPTM','BLOOD')
	colnames(dat)[startsWith(colnames(dat), 'Z')]  <- zcols

} else if (twas == 'AD_TWAS2'){
	# ADTWAS2
	## checked order from submission script
	# zcols <- paste0('Z_', c('TIGAR_ROSMAP_DLPFC','TIGAR_GTEx_BRNCTXB','PrediXcan_GTEx_BRNCTXB','PrediXcan_ROSMAP_DLPFC'))
	zcols <- c('TIGAR_ROSMAP_DLPFC','TIGAR_GTEx_BRNCTXB','PrediXcan_GTEx_BRNCTXB','PrediXcan_ROSMAP_DLPFC')
	colnames(dat)[startsWith(colnames(dat), 'Z')] <- zcols
}

dat <- dat[, c('TargetID', zcols)]
dat2 <- reshape::melt(dat, measure.vars=zcols, variable_name='Z')



if (twas=='AD_TWAS2') {
	dat2$cohort <- ''
	dat2$tissue <- ''

	dat2[dat2$Z == 'TIGAR_ROSMAP_DLPFC', 'cohort'] <- 'TIGAR'
	dat2[dat2$Z == 'TIGAR_ROSMAP_DLPFC', 'tissue'] <- 'ROSMAP_DLPFC'

	dat2[dat2$Z == 'PrediXcan_ROSMAP_DLPFC', 'cohort'] <- 'PrediXcan'
	dat2[dat2$Z == 'PrediXcan_ROSMAP_DLPFC', 'tissue'] <- 'ROSMAP_DLPFC'

	dat2[dat2$Z == 'TIGAR_GTEx_BRNCTXB', 'cohort'] <- 'TIGAR'
	dat2[dat2$Z == 'TIGAR_GTEx_BRNCTXB', 'tissue'] <- 'GTEx_BRNCTXB'

	dat2[dat2$Z == 'PrediXcan_GTEx_BRNCTXB', 'cohort'] <- 'PrediXcan'
	dat2[dat2$Z == 'PrediXcan_GTEx_BRNCTXB', 'tissue'] <- 'GTEx_BRNCTXB'

	dat2$Z <- factor(dat2$Z, 
		levels=c('PrediXcan_ROSMAP_DLPFC', 'PrediXcan_GTEx_BRNCTXB', 'TIGAR_ROSMAP_DLPFC', 'TIGAR_GTEx_BRNCTXB'), 
		labels=c('PrediXcan_ROSMAP_DLPFC', 'PrediXcan_GTEx_BRNCTXB', 'TIGAR_ROSMAP_DLPFC', 'TIGAR_GTEx_BRNCTXB'))
}


## plot
pz <- ggplot(dat2, aes(x=value)) +
	geom_histogram(stat='bin', bins=15, alpha=0.75, color='#218576', fill='#218576') +
	theme_bw() +
	theme(legend.position=legend_pos, 
		legend.title=element_blank(),
		legend.margin=margin(0, 0, 0, 0)
		) +
	labs( 
		x = bquote('Zeta ('*zeta*')'), 
		fill = '', color='') + facet_wrap(~Z)

ggsave(paste0(out_dir, 'zeta_hist_facetwrap_SR_', twas, suffix, '.png'), pz, width=7, height=5)


