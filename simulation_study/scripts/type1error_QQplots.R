#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=FALSE)


sim_dir=#insert sim_dir
scenario_sim_i_suffix = '0.1_0.1_0.5_0.1_0.1_667'


###############
# load libraries
library(fastqq)
library(ggplot2)
library(paletteer)
library(reshape2)



######################################################
# directory for simulation files
wkdir <- paste0(sim_dir, 'type1_error/', scenario_sim_i_suffix, '/')

header <- colnames(read.table(paste0(wkdir, 'results/header.txt'), header=TRUE, check.names=FALSE))

# c('i', 'TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg')
header[header=='SR'] <- 'SR-TWAS'
header[header=='Avg'] <- 'Avg-valid+SR'

dat_raw <- read.table(paste0(wkdir, 'results/pred_sig_results.txt'), header=FALSE, colClasses=c('character', rep('numeric', length(header)-1)))
colnames(dat_raw) <- header

dat <- melt(dat_raw, measure.vars=c('PrediXcan-GTEx', 'TIGAR-ROSMAP', 'TIGAR-ROSMAP_valid', 'Naive', 'SR-TWAS', 'Avg-valid+SR'),
	variable.name='model', value.name='pvalue', na.rm=TRUE)

dat_raw <- NULL


##############
labels_vec <- c('PrediXcan-GTEx', 'TIGAR-ROSMAP', 'Naive', 'TIGAR-ROSMAP_valid', 'SR-TWAS', 'Avg-valid+SR')
breaks_vec <- labels_vec

palette <- 'ggsci::default_jama' # c('black','orange','blue','red','green','purple')
pal0 <- setNames(c(paletteer_d(palette, n=6)), c('black','orange','blue','red','green','purple'))
color_pal <- setNames(c(paletteer_d(palette, n=6)), c('Naive','TIGAR-ROSMAP','PrediXcan-GTEx','SR-TWAS','TIGAR-ROSMAP_valid', 'Avg-valid+SR'))[labels_vec]
shape_pal <- setNames(c(0, 2, 6, 3, 1, 5), c('PrediXcan-GTEx', 'TIGAR-ROSMAP', 'TIGAR-ROSMAP_valid', 'Naive', 'SR-TWAS', 'Avg-valid+SR'))[labels_vec]


##############
get_model_dat <- function(model){
	ci <- 0.95

	obs <- -log10(sort(dat[dat$model == model, 'pvalue']))
	n_obs <- length(obs)
	expect <- -log10(ppoints(n_obs))

	pdat0 <- drop_dense(obs, expect)
	colnames(pdat0) <- c('obs', 'expect')

	n_obs <- nrow(pdat0)
	max_expect <- max(pdat0$expect)

	print(paste0('model: ', model, ', n_obs: ', n_obs))

	pdat <- data.frame(
		model=model,
		obs=pdat0$obs,
		expect=pdat0$expect,
		clower=-log10(qbeta(p=(1 - ci) / 2,
			shape1=seq(n_obs),
			shape2=rev(seq(n_obs)))),
		cupper=-log10(qbeta(
			p=(1 + ci) / 2,
			shape1=seq(n_obs),
			shape2=rev(seq(n_obs)))),
		eq_max_expect=0
	)
	pdat[pdat$expect==max_expect, 'eq_max_expect'] <- 1
	# dat <<-  dat[dat$model != model,]

	return(pdat)
}


##############
plot_dat <- data.frame()

for (model in header[-1]){
	model_dat <- get_model_dat(model)

	if (nrow(plot_dat)==0){
		plot_dat <- model_dat
	} else {
		plot_dat <- rbind(plot_dat, model_dat)
	}
}

plot_dat$model <- factor(plot_dat$model, levels=c('PrediXcan-GTEx', 'TIGAR-ROSMAP', 'TIGAR-ROSMAP_valid', 'Naive', 'SR-TWAS', 'Avg-valid+SR'))

head(plot_dat)


##############
qp <- ggplot(plot_dat, 
		aes(x=expect, 
			y=obs, 
			color=model)) + 
	geom_point(size=1, alpha=0.6) +	
	geom_segment(
		data=plot_dat[plot_dat$eq_max_expect==1,],
		aes(x=0,
			xend=expect,
			y=0,
			yend=expect),
		alpha=0.5,
		color='grey30',
		lineend='round') +
	scale_color_manual(
		name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +
	facet_wrap(~ model, nrow=2) +
	theme_bw() + 
	theme(
		aspect.ratio=1,
		legend.position='bottom') + 
	guides(color=guide_legend(override.aes=list(alpha=0.85, size=3))) + 
	labs(x=bquote("Expected "-"log"[10]("p-value")), 
		y=bquote("Observed "-"log"[10]("p-value")))


ggsave(paste0(wkdir, 'plots/', 'type1_error_qqplot.pdf'), qp, width=6.25, height=6.25)


