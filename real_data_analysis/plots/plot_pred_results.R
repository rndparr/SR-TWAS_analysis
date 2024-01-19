#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=F, digits=6, scipen=10)

# libraries
library(ggplot2)
library(reshape2)

## directories, input depend on whith test data set
test_dat <- 'MAP'
ROSMAP_train_cohort <- list('MAP'='ROS', 'ROS'='MAP')[[test_dat]]

# directory
SR_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/'
out_dir <- paste0(SR_dir, 'plots/')

## read in results data
path <- paste0(SR_dir, 'SR_TWAS_pred/pred_results_', test_dat, '_test_data.txt')
dat <- read.table(path, header=TRUE, check.names=FALSE, fill=TRUE, sep='\t')
colnames(dat)[which(colnames(dat) == 'ROSMAP')] <- ROSMAP_train_cohort

## for finding max R2 value
expr_R2_vals <- as.numeric(unlist(dat[,c(ROSMAP_train_cohort, 'GTEx', 'PredDB', 'Naive', 'SR')]))
R2_max <- max(expr_R2_vals, na.rm=TRUE)

## melt data
datm <- melt(dat, measure.vars = c('PredDB','GTEx',ROSMAP_train_cohort,'Naive'), variable.name='cohort', value.name='x')

datm$cohort <- factor(datm$cohort,
	levels=c('PredDB','GTEx',ROSMAP_train_cohort,'Naive3'),
	labels = c(
		'PrediXcan_GTEx_BRNCTXB',
		'TIGAR_GTEx_BRNCTXB', 
		paste0('TIGAR_', ROSMAP_train_cohort, '_DLPFC'), 
		paste0('Naive_', test_dat, '_DLPFC'))
	)

## color by which R^2 > cutoff
# ggplot2 default color palette
r2_cut <- 0.005

datm$color_sort <- NA
datm[datm$x <= r2_cut & datm$SR <= r2_cut, 'color_sort'] <- 4 #Neither
datm[datm$x > r2_cut & datm$SR > r2_cut, 'color_sort'] <- 3 #Both
datm[datm$x > r2_cut & datm$SR <= r2_cut, 'color_sort'] <- 2 #Other
datm[datm$x <= r2_cut & datm$SR > r2_cut, 'color_sort'] <- 1 #SR
datm <- datm[order(datm$color_sort, decreasing=TRUE),]


# default values for cohorts to use
all_cohorts <- levels(datm$cohort)


# functions
R2_plot <- function(data, 
		cohorts=all_cohorts, 
		plot_title=bquote('Expression prediction'~italic(R)^{2}~'for'~.(test_dat)~'test data'), 
		theme=theme_gray() + theme(strip.background=element_blank(),
		strip.placement='outside'),
		legend.pos='bottom',
		strip_pos='bottom',
		size=1,
		fixed=TRUE){

	dataf <- data[data$cohort %in% cohorts, ]

	R2_plt <- ggplot(dataf, aes(x=x, y=SR, color=color)) + 
		geom_point(alpha=0.5, size=1.25*size) + 
		geom_abline(intercept=0, slope=1, color='black', linetype='dashed', alpha=1)
	if (fixed){
		R2_plt <- R2_plt + coord_fixed(ratio=1, 
		xlim=c(0, R2_max), 
		ylim=c(0, R2_max)) 
	} else {
		R2_plt <- R2_plt + xlim(0, R2_max) +
			ylim(0, R2_max)
	}
		R2_plt <- R2_plt + facet_wrap(. ~ cohort, 
			ncol=2,
			strip.position=strip_pos, 
			labeller=label_parsed) +
		labs(
			title=plot_title, 
			x=NULL,
			y=paste0('SR-TWAS_', test_dat, '_DLPFC'),
			fill='') +
		theme +
		theme(legend.pos=legend.pos,
		legend.margin=margin(0, 0, 0, 0)
		) +
		theme +
		scale_color_identity(name=bquote(italic(R)^2 > .(r2_cut)), 
			labels=c('SR-TWAS Only', 'Other Model Only', 'Both', 'Neither'),guide='legend') +
		guides(colour = guide_legend(override.aes = list(size=2*size, alpha=0.9)))

	return(R2_plt)
}

####
siz <- 6
ar <- 1

# set up palette
# green red orange palette
pal <- c(
	'1'='#A41720', #SR
	'2'='#F86726', #Other
	'3'='#218576', #Both
	'4'='#1A2423' #Neither
	)
datm$color <- pal[datm$color_sort]
datm$color <- factor(datm$color, levels=pal)

R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw()+ theme(legend.position='bottom',legend.margin=margin(0, 0, 0, 0)), strip_pos='top',size=1.5,legend.pos='right',fixed=FALSE) + ggtitle(paste0(test_dat, ' Test Data'))
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '.pdf'), R2qq + ggtitle(NULL), width=5.5, height=6, units='in')

