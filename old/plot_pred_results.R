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
# test_dat <- 'ROS'
ROSMAP_train_cohort <- list('MAP'='ROS', 'ROS'='MAP')[[test_dat]]

# cluster
SR_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/'

# local
# SR_dir <- '~/HGCC/YangFSSdata/SR_TWAS/'

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
	levels=c('PredDB','GTEx',ROSMAP_train_cohort,'Naive'),
	labels = c(
		'PrediXcan_GTEx_BRNCTXB',
		'TIGAR_GTEx_BRNCTXB', 
		paste0('TIGAR_', ROSMAP_train_cohort, '_DLPFC'), 
		paste0('Naive_', ROSMAP_train_cohort, '_DLPFC'))
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
R2_plot_bw <- function(data, cohorts=all_cohorts, plot_title=bquote('Expression prediction'~italic(R)^{2}~'for'~.(test_dat)~'test data')){

	dataf <- data[data$cohort %in% cohorts, ]

	R2_plt <- ggplot(dataf, aes(x=x, y=SR)) + 
		geom_point(color='black', alpha=0.3, size=1) + 
		geom_abline(intercept=0, slope=1, color='red') +
		coord_fixed(ratio=1, 
			xlim=c(0, R2_max), 
			ylim=c(0, R2_max)) + 
		facet_wrap(. ~ cohort, 
			# ncol=length(levels(data$cohort)), 
			ncol=2,
			strip.position='bottom', 
			labeller=label_parsed) +
		labs(
			title=plot_title,
			x=NULL,
			# y=bquote('SR-TWAS'~italic(R)^2),
			y=paste0('SR-TWAS_', ROSMAP_train_cohort, '_DLPFC'),
			fill='') +
		theme(strip.background=element_blank(),
			strip.placement='outside')
	
	return(R2_plt)
}

R2_plot <- function(data, cohorts=all_cohorts, plot_title=bquote('Expression prediction'~italic(R)^{2}~'for'~.(test_dat)~'test data'), 
	theme=theme_gray() + theme(strip.background=element_blank(),
		strip.placement='outside'),
	legend.pos='bottom',
	strip_pos='bottom',
	size=1, fixed=TRUE){

	dataf <- data[data$cohort %in% cohorts, ]

	R2_plt <- ggplot(dataf, aes(x=x, y=SR, color=color)) + 
		# geom_abline(intercept=0, slope=1, color='black', linetype='dashed') +
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
			# ncol=length(levels(data$cohort)), 
			ncol=2,
			strip.position=strip_pos, 
			labeller=label_parsed) +
		labs(
			title=plot_title, 
			x=NULL,
			# y=bquote('Stacked Regression'~italic(R)^2),
			# y=bquote('SR-TWAS'~italic(R)^2),
			y=paste0('SR-TWAS_', ROSMAP_train_cohort, '_DLPFC'),
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
# ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_blue_ashg2.pdf'), R2qq, width=6, height=6.5, units='in')
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '.pdf'), R2qq + ggtitle(NULL), width=5.5, height=6, units='in')

# R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw()+ theme(strip.text=element_text(color='white', face='bold',size=rel(1.1)), strip.text.x = element_text(color='white', face='bold'), strip.background=element_rect(fill='#09597b', colour='#09597b'), legend.position='bottom',legend.margin=margin(0, 0, 0, 0)), strip_pos='top',size=1.5,legend.pos='right',fixed=FALSE) + ggtitle(paste0(test_dat, ' Test Data'))
# # ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_blue_ashg2.pdf'), R2qq, width=6, height=6.5, units='in')
# ggsave(paste0(out_dir, 'pred_R2_', test_dat, '.pdf'), R2qq + ggtitle(NULL), width=6, height=5, units='in')



################
# function for palette
gg_color_hue <- function(n) {hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]}

pal <- setNames(gg_color_hue(4), 1:4)

datm$color <- pal[datm$color_sort]
datm$color <- factor(datm$color, levels=pal)

R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw(), strip_pos='top',size=1.5,legend.pos='right',fixed=FALSE)
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '.pdf'), R2qq, width=6, height=5, units='in')



############


pal <- c(
	'1'='#62a323', #SR
	'2'='#c4af25', #Other
	'3'='#80c3e9', #Both
	'4'='#09597b' #Neither
	)


###############

out_dir <- paste0(SR_dir, 'plots/ashg/')

# R2qq <- R2_plot(datm, theme=theme_bw(), strip_pos='top',size=1.5)
# ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_ashg.pdf'), R2qq, width=siz, height=siz*ar, units='in')

R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw(), strip_pos='top',size=1.5,legend.pos='right',fixed=FALSE) + ggtitle(paste0(test_dat, ' Test Data'))
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_ashg.pdf'), R2qq, width=6, height=5, units='in')


R2qq <- R2_plot(datm, theme=theme_bw(), strip_pos='top',size=1.5)
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_ashg.png'), R2qq, width=siz, height=siz*ar, units='in')

R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw(), strip_pos='top',size=1.5)
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_ashg.png'), R2qq, width=siz, height=siz*ar, units='in')

### BLUE
R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw()+ theme(strip.text = element_text(color='white', face='bold'), strip.background=element_rect(fill='#09597b', colour='#09597b'),), strip_pos='top',size=1.5,legend.pos='right',fixed=FALSE) + ggtitle(paste0(test_dat, ' Test Data'))
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_blue_ashg.pdf'), R2qq, width=6, height=5, units='in')


R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw()+ theme(strip.text=element_text(color='white', face='bold',size=rel(1.1)), strip.text.x = element_text(color='white', face='bold'), strip.background=element_rect(fill='#09597b', colour='#09597b'), legend.position='bottom',legend.margin=margin(0, 0, 0, 0)), strip_pos='top',size=1.5,legend.pos='right',fixed=FALSE) + ggtitle(paste0(test_dat, ' Test Data'))
# ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_blue_ashg2.pdf'), R2qq, width=6, height=6.5, units='in')
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_blue_ashg2.pdf'), R2qq + ggtitle(NULL), width=5.5, height=5, units='in')

#####


## old
# siz <- 10
# ar <- 4.5 / 8

## 3 panels
# siz <- 10
# ar <- 0.422932

# 4 panels
# siz <- 8
# ar <- 0.316486

siz <- 6
ar <- 1.1


R2qq_bw <- R2_plot_bw(datm)
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_bw.png'), R2qq_bw, width=siz, height=siz, units='in')


R2qq_bw <- R2_plot_bw(datm, plot_title=NULL)
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_bw_notitle.png'), R2qq_bw, width=siz, height=siz, units='in')


R2qq <- R2_plot(datm)
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '.png'), R2qq, width=siz, height=siz*ar, units='in')


R2qq <- R2_plot(datm, plot_title=NULL)
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_notitle.png'), R2qq, width=siz, height=siz*ar, units='in')

####
# ggplot2 colors
# pal <- c(
# 	'1'='#00BFC4', #SR
# 	'2'='#7CAE00', #Other
# 	'3'='#F8766D', #Both
# 	'4'='#C77CFF' #Neither
# 	)

pal <- c(
	'1'='#62a323', #SR
	'2'='#c4af25', #Other
	'3'='#80c3e9', #Both
	'4'='#09597b' #Neither
	)

datm$color <- pal[datm$color_sort]
datm$color <- factor(datm$color, levels=pal)


R2qq <- R2_plot(datm, plot_title=NULL, theme=theme_bw())
ggsave(paste0(out_dir, 'pred_R2_', test_dat, '_ashg.png'), R2qq, width=siz, height=siz*ar, units='in')


############

# function for palette
gg_color_hue <- function(n) {hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]}

# color by which R^2 > cutoff
r2_cut <- 0.005
datm$color <- ''
datm[datm$x <= r2_cut & datm$SR <= r2_cut, 'color'] <- '4'
datm[datm$x <= r2_cut & datm$SR > r2_cut, 'color'] <- '3'
datm[datm$x > r2_cut & datm$SR <= r2_cut, 'color'] <- '2'
datm[datm$x > r2_cut & datm$SR > r2_cut, 'color']  <- '1'
datm$color <- as.factor(datm$color)

# functions
R2_plot <- function(data, cohorts=all_cohorts, plot_title=bquote('Expression prediction'~italic(R)^{2}~'for'~.(test_dat)~'test data')){
	dataf <- data[data$cohort %in% cohorts, ]

	R2_plt <- ggplot(dataf, aes(x=x, y=SR, color=color)) + 
		# geom_abline(intercept=0, slope=1, color='black', linetype='dashed') +
		geom_point(alpha=0.5, size=1.25) + 
		geom_abline(intercept=0, slope=1, color='black', linetype='dashed', alpha=1) +
		coord_fixed(ratio=1, 
			xlim=c(0, R2_max), 
			ylim=c(0, R2_max)) + 
		facet_wrap(. ~ cohort, 
			# ncol=length(levels(data$cohort)), 
			ncol=2,
			strip.position='bottom', 
			labeller=label_parsed) +
		labs(
			title=plot_title, 
			x=NULL,
			# y=bquote('Stacked Regression'~italic(R)^2),
			# y=bquote('SR-TWAS'~italic(R)^2),
			y='SR-TWAS',
			fill='') +
		theme(strip.background=element_blank(),
			strip.placement='outside',
			legend.pos='bottom',
			legend.margin=margin(0, 0, 0, 0)
			) +
		# scale_color_manual(name=bquote(italic(R)^2 > .(r2_cut)), labels=setNames(c('Both', 'Single Cohort Only', 'ST Only', 'Neither'),  as.character(1:4)), values=gg_color_hue(4)) 
		# scale_color_manual(name=bquote(italic(R)^2 > .(r2_cut)),  breaks = c(3,2,1,4), labels=setNames(c('Both', 'Other Model Only', 'SR-TWAS Only', 'Neither'),  as.character(1:4)), values=gg_color_hue(4))
	scale_color_manual(name=bquote(italic(R)^2 > .(r2_cut)),  breaks = as.character(c(3,2,1,4)),
 labels=setNames(c('SR-TWAS Only', 'Other Model Only', 'Both','Neither'),  as.character(c(3,2,1,4))), values=gg_color_hue(4)) + 
	guides(colour = guide_legend(override.aes = list(size=2)))

	return(R2_plt)
}






# R2qq_bw <- R2_plot_bw(datm)

# ggsave(paste0(out_dir, 'all_pred_R2_', test_dat,'_bw.png'), R2qq_bw, width=siz, height=siz, units='in')


# R2qq_bw <- R2_plot_bw(datm, plot_title=NULL)

# ggsave(paste0(out_dir, 'all_pred_R2_', test_dat,'_bw_notitle.png'), R2qq_bw, width=siz, height=siz, units='in')


# R2qq <- R2_plot(datm)

# ggsave(paste0(out_dir, 'all_pred_R2_', test_dat,'.png'), R2qq, width=siz, height=siz*ar, units='in')


# R2qq <- R2_plot(datm, plot_title=NULL)

# ggsave(paste0(out_dir, 'all_pred_R2_', test_dat,'_notitle.png'), R2qq, width=siz, height=siz*ar, units='in')

#########
datm2 <- melt(dat, measure.vars = c(ROSMAP_train_cohort,'GTEx','PredDB','Naive','SR'), variable.name='cohort', value.name='x')

r2_cut <- 0.005


# only those greater greater than r2_cut
datm2_f1 <- datm2[datm2$x > r2_cut, ]
merge(aggregate(x ~ cohort, data=datm2_f1, median),
	aggregate(x ~ cohort, data=datm2_f1, mean),
	by='cohort')

# only targets where at least one of SR or other cohort is > r2_dat
targets_r2 <- datm[datm$color != '4', 'TargetID']
datm2_f2 <- datm2[datm2$TargetID %in% targets_r2, ]
merge(aggregate(x ~ cohort, data=datm2_f2, median),
	aggregate(x ~ cohort, data=datm2_f2, mean),
	by='cohort')

## MAP test
#    model    median       mean
# 1   GTEx 0.0000000 0.01352386
# 2  Naive 0.0169739 0.05240389
# 3 PredDB 0.0000000 0.00944142
# 4    ROS 0.0145148 0.04902433
# 5     SR 0.0172128 0.05514766

## ROS test
#    model    median       mean
# 1   GTEx 0.0000000 0.01287217
# 2    MAP 0.0139594 0.04628768
# 3  Naive 0.0160991 0.04973874
# 4 PredDB 0.0000000 0.00885992
# 5     SR 0.0161520 0.05199054

merge(
	merge(aggregate(x ~ cohort, data=datm2_f2, median),
		aggregate(x ~ cohort, data=datm2_f2, mean),
		by='cohort'),
	aggregate(x ~ cohort, data=datm2_f1, length),
	by='cohort')

# ## MAP test
#    model    median       mean    N
# 1   GTEx 0.0000000 0.01352386 2405
# 2  Naive 0.0169739 0.05240389 8364
# 3 PredDB 0.0000000 0.00944142  867
# 4    ROS 0.0145148 0.04902433 7913
# 5     SR 0.0172128 0.05514766 8425

# ## ROS test
#    model    median       mean    N
# 1   GTEx 0.0000000 0.01287217 2327
# 2    MAP 0.0139594 0.04628768 7793
# 3  Naive 0.0160991 0.04973874 8212
# 4 PredDB 0.0000000 0.00885992  858
# 5     SR 0.0161520 0.05199054 8163


##########
## MAX VALUE STUFF

# get max value columns
dat$max_GTEx_ROSMAP <- pmax(dat[['GTEx']], dat[[ROSMAP_train_cohort]], na.rm = TRUE)
dat$max_GTEx_ROSMAP <- pmax(dat[['GTEx']], dat[[ROSMAP_train_cohort]], na.rm = TRUE)
dat$max_GTEx_ROSMAP_Naive <- pmax(dat$max_GTEx_ROSMAP, dat[['Naive']], na.rm = TRUE)

# melt
datm <- melt(dat, measure.vars = c(ROSMAP_train_cohort,'GTEx','PredDB','Naive','max_GTEx_ROSMAP','max_GTEx_ROSMAP_Naive'), variable.name='cohort', value.name='x')
# datm$cohort <- factor(datm$cohort,
# 	labels = c(
# 		bquote('GTEx' ~ italic(R)^2), 
# 		bquote(.(ROSMAP_train_cohort) ~ italic(R)^2), 
# 		bquote('max(GTEx' ~ italic(R)^2 * ','~ .(ROSMAP_train_cohort) ~ italic(R)^2 * ')'),
# 		bquote('Naive' ~ italic(R)^2),
# 		bquote('max(GTEx' ~ italic(R)^2 * ',' ~ .(ROSMAP_train_cohort) ~ italic(R)^2 * ', Naive'~italic(R)^{2}*')')))
datm$cohort <- factor(datm$cohort,
	labels = c(
		paste0('TIGAR_', ROSMAP_train_cohort), 
		'TIGAR_GTEx', 
		'PredDB',
		bquote('Naive'),
		bquote('max(TIGAR_GTEx, TIGAR_'* .(ROSMAP_train_cohort) * ')'),
		bquote('max(GTEx,' ~ .(ROSMAP_train_cohort)~ ', Naive'~italic(R)^{2}*')')))

# default values for cohorts to use
all_cohorts <- levels(datm$cohort)
single_cohorts <- all_cohorts[-c(grep('max', all_cohorts))]

##########

# prep data
plot1_dat <- melt(dat, measure.vars=c('GTEx',ROSMAP_train_cohort,'Naive'), variable.name='cohort', value.name='x')
plot1_dat$cohort <- factor(plot1_dat$cohort,
	labels=c(bquote('GTEx' ~ italic(R)^{2}), bquote(.(ROSMAP_train_cohort) ~ italic(R)^{2}), bquote('Naive' ~ italic(R)^{2})))

# plot
plot_title <- bquote('Expression prediction'~italic(R)^{2}~'for'~.(test_dat)~'test data')
R2qq <- ggplot(plot1_dat, aes(x = x, y = SR)) + 
	geom_point(color = 'black', alpha=0.5) + 
	geom_abline(intercept=0, slope=1, color = 'red') +
	coord_fixed(ratio = 1, 
		xlim = c(0, max(expr_R2_vals, na.rm=TRUE)), 
		ylim = c(0, max(expr_R2_vals, na.rm=TRUE))) + 
	facet_wrap(. ~ cohort, #scales = 'free',
		ncol=3, strip.position='bottom', 
		labeller = label_parsed #label_bquote(cols = .(cohort)~R^2)
		# labeller = as_labeller(c(
		# 	GTEx = bquote('GTEx'~R^2), 
		# 	ROS = bquote('ROS'~R^2)))
		) +
	labs(title = plot_title, 
		x = NULL,
		y = bquote('Stacked Regression'~italic(R)^{2}),
		fill = '') +
	theme(strip.background = element_blank(),
		strip.placement = 'outside')

siz <- 10
ar <- 4.5 / 8

ggsave(paste0(out_dir, 'all_pred_R2_', test_dat,'.png'), R2qq, width=13.3, height=siz*ar, units='in')
# ggsave(paste0(app_dir , 'results/all_pred_R2_', test_dat,'.png'), R2qq, width=siz, height=siz*ar, units='in')


###
plot2_dat <- dat
plot2_dat$max_GTEx_ROSMAP <- pmax(plot2_dat[['GTEx']], plot2_dat[[ROSMAP_train_cohort]], na.rm = TRUE)

plot2_dat2 <- melt(plot2_dat, measure.vars = c('GTEx','ROS','max_GTEx_ROSMAP','Naive'), variable.name='cohort', value.name='x')

plot2_dat2$cohort <- factor(plot2_dat2$cohort,
	labels = c(
		bquote('GTEx' ~ italic(R)^2), 
		bquote(.(ROSMAP_train_cohort) ~ italic(R)^2), 
		bquote('max(GTEx' ~ italic(R)^2 * ','~ .(ROSMAP_train_cohort) ~ italic(R)^2 * ')'),
		bquote('Naive' ~ italic(R)^{2})))

# function for palette
gg_color_hue <- function(n) {hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]}

# color by which R^2 > cutoff
r2_cut <- 0.005
plot2_dat$color <- ''
plot2_dat[plot2_dat$x <= r2_cut & plot2_dat$ST <= r2_cut, 'color'] <- '4'
plot2_dat[plot2_dat$x <= r2_cut & plot2_dat$ST > r2_cut, 'color'] <- '3'
plot2_dat[plot2_dat$x > r2_cut & plot2_dat$ST <= r2_cut, 'color'] <- '2'
plot2_dat[plot2_dat$x > r2_cut & plot2_dat$ST > r2_cut, 'color']  <- '1'
plot2_dat$color <- as.factor(plot2_dat$color)




