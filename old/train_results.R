#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE, digits=7)

#####################
## LIBRARIES
library(ggplot2)
library(gridExtra)
library(reshape2)

## directory
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'

## training model
train_dat <- 'ROSMAP'
train_dat <- 'ROS_MAP'

## models
models <- c('PredDB', 'GTEx', unlist(strsplit(train_dat,'_')), paste0('Naive_', train_dat), paste0('SR_', train_dat))
models2 <- setNames(list('PrediXcan_GTEx', 'TIGAR_GTEx', paste0('TIGAR_', unlist(strsplit(train_dat,'_'))), 'Naive', 'SR-TWAS'), c('PredDB', 'GTEx', unlist(strsplit(train_dat,'_')), paste0('Naive_', train_dat), paste0('SR_', train_dat)))

## paths
paths <- setNames(list(
	PredDB=function(chr) paste0('/mnt/YangFSS/data2/rparrish/PredDB_TIGAR_format/eQTL_weights/Brain_Frontal_Cortex_BA9/CHR', chr, '_PredictDB_GeneInfo.txt'),
	GTEx=function(chr) paste0('/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Brain_Frontal_Cortex_BA9/DPR_CHR', chr, '/CHR', chr, '_DPR_train_GeneInfo.txt'),
	ROSMAP=function(chr) paste0('/mnt/YangFSS/data2/rparrish/ROSMAP_WGS/TIGAR_train_', train_dat, '/CHR', chr, '_DPR_train_', train_dat, '_WGS_b38_GeneInfo.txt'),
	Naive=function(chr) paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/Brain_', train_dat, '_GTEx_PredDB/CHR', chr, '_naive_train_GeneInfo.txt'),
	SR=function(chr) paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/Brain_', train_dat, '_GTEx_PredDB/CHR', chr, '_SR_train_GeneInfo.txt')), models)


## Read in data
dat <- data.frame()
for (chr in 1:22) {
	# get train info data
	for (model in models){
		# get data; grep the target id from tabix output, read as table with column names from header
		raw <- read.csv(paths[[model]](chr), sep='\t')

		if (startsWith(model, 'SR') | startsWith(model, 'Naive')){
			colnames(raw) <- tolower(colnames(raw) )
			raw <- raw[,tolower(c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','N_SNP','N_Effect_SNP','CVR2','R2'))]
			colnames(raw) <- c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainR2')
		} else {
			raw <- raw[,c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainR2')]
		}

		# model factor name
		raw$model <- models2[[model]]
	
		dat <- rbind(dat, raw)
	}
}

# set factor order
dat$model <- factor(dat$model, levels=models2[models2 %in% unique(dat$model)])



p1_dat <- dat[, c('model','TargetID','CVR2','TrainR2')]
p1_dat_melt <- reshape2::melt(
		p1_dat, 
		measure.vars = c('CVR2', 'TrainR2'))
p1_dat_melt <- na.omit(p1_dat_melt)
p1_dat_melt$variable <- factor(
	p1_dat_melt$variable, 
	labels = c(bquote('Average 5-Fold CV' ~ italic(R)^2),
		bquote('Training' ~ italic(R)^2)
	),
	levels = c(
		'CVR2',
		'TrainR2'
	))

p1 <- ggplot(p1_dat_melt, 
				aes(x = value, color = model)) +  
	geom_line(stat='Density', alpha = 0.8, size=0.5) +
	labs( 
		x = bquote(italic(R)^2), 
		y = 'Density', 
		fill = '', color='') +
	facet_wrap(variable ~ .,
	 # . ~ variable ,
		scales = 'free', 
		# nrow=2,
		# ncol=2,
		ncol=1,
		labeller = label_parsed
		) +
	theme(legend.position = 'bottom')

	#  + 
	# scale_color_manual(name = '', labels = label_names, breaks = names(pal), values = pal) + 
	# scale_linetype_manual(name = '', labels = c('model training', 'base model in SR-TWAS training'), values = c('solid', 'dotted'))


ggsave(paste0(out_dir, train_dat, '_train_R2_density.png'), p1, width=7, height=7)


