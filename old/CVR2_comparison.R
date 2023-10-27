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
# train_dat <- 'ROS_MAP'

## models
# models <- c('PredDB', 'GTEx', unlist(strsplit(train_dat,'_')), paste0('Naive_', train_dat), paste0('SR_', train_dat))
# models2 <- setNames(list('PrediXcan_GTEx', 'TIGAR_GTEx', paste0('TIGAR_', unlist(strsplit(train_dat,'_'))), 'Naive', 'SR-TWAS'), c('PredDB', 'GTEx', unlist(strsplit(train_dat,'_')), paste0('Naive_', train_dat), paste0('SR_', train_dat)))

models <- c('Naive', 'SR')

## paths
paths <- setNames(list(
	PredDB=function(chr) paste0('/mnt/YangFSS/data2/rparrish/PredDB_TIGAR_format/eQTL_weights/Brain_Frontal_Cortex_BA9/CHR', chr, '_PredictDB_GeneInfo.txt'),
	GTEx=function(chr) paste0('/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Brain_Frontal_Cortex_BA9/DPR_CHR', chr, '/CHR', chr, '_DPR_train_GeneInfo.txt'),
	ROSMAP=function(chr) paste0('/mnt/YangFSS/data2/rparrish/ROSMAP_WGS/TIGAR_train_', train_dat, '/CHR', chr, '_DPR_train_', train_dat, '_WGS_b38_GeneInfo.txt'),
	Naive=function(chr) paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/Brain_', train_dat, '_GTEx_PredDB/CHR', chr, '_naive3_train_GeneInfo.txt'),
	SR=function(chr) paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/Brain_', train_dat, '_GTEx_PredDB/CHR', chr, '_SR_train_GeneInfo.txt')), models)


## Read in data
dat <- data.frame()
for (chr in 1:22) {
	# get train info data
	for (model in models){
		# get data; grep the target id from tabix output, read as table with column names from header
		raw <- read.csv(paths[[model]](chr), sep='\t')
		if (startsWith(model, 'SR') | startsWith(model, 'Naive')){
			colnames(raw) <- tolower(colnames(raw))
			search_cols <- tolower(c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','N_SNP','N_EFFECT_SNP','CVR2','TrainR2'))
			raw <- raw[, search_cols]
			colnames(raw) <- c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainR2')
		} else {
			raw <- raw[,c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainR2')]
		}

		# model factor name
		# raw$model <- models2[[model]]
		raw$model <- model
	
		dat <- rbind(dat, raw)
	}
}

# set factor order
dat$model <- as.factor(dat$model)

# get CVR2 max
CVR2_max <- max(dat$CVR2, na.rm=TRUE)
dat2 <- dcast(dat, GeneName + TargetID ~ model, value.var='CVR2')


# number of genes each model
table(dat$model)
# Naive    SR 
#  4563 21921 

# summary stats
by(dat[,c('CVR2','TrainR2','n_snp','n_effect_snp','model')], dat[,'model'], summary)
# dat[, "model"]: Naive
#       CVR2             TrainR2            n_snp       n_effect_snp   
#  Min.   :-1.35513   Min.   :0.02764   Min.   :   3   Min.   :  1.00  
#  1st Qu.: 0.01806   1st Qu.:0.16333   1st Qu.:1276   1st Qu.: 14.00  
#  Median : 0.06109   Median :0.25788   Median :1656   Median : 25.00  
#  Mean   : 0.10086   Mean   :0.28480   Mean   :1690   Mean   : 31.85  
#  3rd Qu.: 0.14951   3rd Qu.:0.38034   3rd Qu.:2036   3rd Qu.: 43.00  
#  Max.   : 0.67802   Max.   :0.97185   Max.   :5621   Max.   :199.00  
#    model     
#  Naive:4563  
#  SR   :   0  
# ------------------------------------------------------------ 
# dat[, "model"]: SR
#       CVR2            TrainR2            n_snp        n_effect_snp  
#  Min.   :0.01001   Min.   :0.02782   Min.   :   95   Min.   :   95  
#  1st Qu.:0.02344   1st Qu.:0.21641   1st Qu.: 5031   1st Qu.: 5031  
#  Median :0.03806   Median :0.29686   Median : 5985   Median : 5985  
#  Mean   :0.06528   Mean   :0.30892   Mean   : 6076   Mean   : 6076  
#  3rd Qu.:0.06745   3rd Qu.:0.38884   3rd Qu.: 6880   3rd Qu.: 6880  
#  Max.   :0.69530   Max.   :0.86800   Max.   :29383   Max.   :29383  
#    model      
#  Naive:    0  
#  SR   :21921

naive_genes <- dat[(dat$model == 'Naive') & !(is.na(dat$CVR2)), 'GeneName']
sr_genes <- dat[(dat$model == 'SR') & !(is.na(dat$CVR2)), 'GeneName']

# genes in both
length(intersect(naive_genes, sr_genes)) # 4378
# sr only genes
length(setdiff(sr_genes, naive_genes)) # 17523
# naive only genes
length(setdiff(naive_genes, sr_genes)) # 185


# #######
# bleh <- function(x){ median(x[,'CVR2'],na.rm=TRUE) }
# bleh <- function(x){print(x)}
# ## for finding max R2 value
# expr_R2_vals <- as.numeric(unlist(dat[,c(ROSMAP_train_cohort, 'GTEx', 'PredDB', 'Naive', 'SR')]))
# R2_max <- max(expr_R2_vals, na.rm=TRUE)

# ## melt data
# datm <- melt(dat, measure.vars = c('PredDB','GTEx',ROSMAP_train_cohort,'Naive'), variable.name='cohort', value.name='x')

# datm$cohort <- factor(datm$cohort,
# 	levels=c('PredDB','GTEx',ROSMAP_train_cohort,'Naive'),
# 	labels = c(
# 		'PrediXcan_GTEx',
# 		'TIGAR_GTEx', 
# 		paste0('TIGAR_', ROSMAP_train_cohort), 
# 		'Naive')
# 	)

## color by which R^2 > cutoff
# ggplot2 default color palette


dat2[is.na(dat2)] <- 0
dat2[which(dat2$Naive < 0), 'Naive'] <- 0
dat2[which(dat2$SR < 0), 'SR'] <- 0

r2_cut <- 0.005

dat2$color_sort <- NA
dat2[dat2$Naive <= r2_cut & dat2$SR <= r2_cut, 'color_sort'] <- 4 #Neither
dat2[dat2$Naive > r2_cut & dat2$SR > r2_cut, 'color_sort'] <- 3 #Both
dat2[dat2$Naive > r2_cut & dat2$SR <= r2_cut, 'color_sort'] <- 2 #Other
dat2[dat2$Naive <= r2_cut & dat2$SR > r2_cut, 'color_sort'] <- 1 #SR
dat2 <- dat2[order(dat2$color_sort, decreasing=TRUE),]


# default values for cohorts to use
# all_cohorts <- levels(datm$cohort)



R2_plot <- function(data, cohorts=all_cohorts, plot_title=bquote('Expression prediction'~italic(R)^{2}~'for'~.(test_dat)~'test data'), 
	theme=theme_gray() + theme(strip.background=element_blank(),
		strip.placement='outside'),
	legend.pos='bottom',
	strip_pos='bottom',
	size=1, fixed=TRUE){

	# dataf <- data[data$cohort %in% cohorts, ]

	R2_plt <- ggplot(data, aes(x=Naive, y=SR, color=color)) + 
		# geom_abline(intercept=0, slope=1, color='black', linetype='dashed') +
		geom_point(alpha=0.5, size=1.25*size) + 
		geom_abline(intercept=0, slope=1, color='black', linetype='dashed', alpha=1)
	if (fixed){
		R2_plt <- R2_plt + coord_fixed(ratio=1, 
		xlim=c(0, CVR2_max), 
		ylim=c(0, CVR2_max)) 
	} else {
		R2_plt <- R2_plt + xlim(0, CVR2_max) +
			ylim(0, CVR2_max)
	}
		R2_plt <- R2_plt + 
		# facet_wrap(. ~ cohort, 
		# 	# ncol=length(levels(data$cohort)), 
		# 	ncol=2,
		# 	strip.position=strip_pos, 
		# 	labeller=label_parsed) +
		labs(
			title=plot_title, 
			x='Naive',
			# y=bquote('Stacked Regression'~italic(R)^2),
			# y=bquote('SR-TWAS'~italic(R)^2),
			y='SR-TWAS',
			fill='') +
		theme +
		theme(legend.pos=legend.pos,
		legend.margin=margin(0, 0, 0, 0)
		) +
		theme +
		scale_color_identity(name=bquote('CV'~italic(R)^2 > .(r2_cut)), 
			labels=c('SR-TWAS Only', 'Naive Only', 'Both','Neither'),guide='legend') +
		guides(colour = guide_legend(override.aes = list(size=2*size, alpha=0.9)))

	return(R2_plt)
}

####

siz <- 6
ar <- 1


# function for palette
gg_color_hue <- function(n) {hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]}

pal <- setNames(gg_color_hue(4), 1:4)

dat2$color <- pal[dat2$color_sort]
dat2$color <- factor(dat2$color, levels=pal)

R2qq <- R2_plot(dat2, plot_title=NULL, theme=theme_bw(), strip_pos='top',size=1.5,legend.pos='right',fixed=FALSE)
ggsave(paste0(out_dir, 'CVR2_', train_dat, '.pdf'), R2qq, width=6.5, height=5, units='in')










#######


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


