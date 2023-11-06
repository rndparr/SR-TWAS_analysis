#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=FALSE, digits=5)

###############
# read in passed arguments
args <- (commandArgs(TRUE))
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==2) {
	sim_dir = as.character(args[[1]])
	suffix = ''
	ggplot2_lib = as.character(args[[2]])
} else {
	sim_dir = as.character(args[[1]])
	suffix = as.character(args[[2]])
	ggplot2_lib = as.character(args[[3]])
}

# load libraries
library(ggplot2, lib.loc=ggplot2_lib)
library(grid)
library(gtable)
library(gridExtra)
library(paletteer)
library(reshape2)


################
# PLOT FUNCTIONS
facet_labs <- function(p, tlab=NULL, rlab=NULL){
	gt <- ggplotGrob(p)

	if (!is.null(tlab)){
		strip <-c(subset(gt$layout, grepl('strip-t', gt$layout$name), select=t:r))
		gt <- gtable_add_rows(gt, unit(2, 'lines'), max(strip$t)-1)
		gt <- gtable_add_grob(gt, 
			textGrob(tlab, rot=0, gp=gpar(cex=0.9)), t=max(strip$t), l=min(strip$l), r=max(strip$r), b=max(strip$b))
	}

	if (!is.null(rlab)){
		strip <-c(subset(gt$layout, grepl('strip-r', gt$layout$name), select=t:r))
		gt <- gtable_add_cols(gt, unit(2, 'lines'), max(strip$r))
		gt <- gtable_add_grob(gt, 
		  textGrob(rlab, rot=-90, gp=gpar(cex=0.9)),  t=min(strip$t), l=max(strip$r)+1, b=max(strip$b))
	}

	return(gt)
}


# EXPRESSION FACET LABEL
get_expr_facet_lab <- function(job_cat) { 
	expr_facet_lab_var1 <- list('same'='', 'half'='ROSMAP ', 'zero'='ROSMAP ', 'overlap'='')[[job_cat]]
	expr_facet_lab_var2 <- list(
		'same'=bquote(italic(h)[italic(e)*",GTEx"]^{2}==italic(h)[italic(e)*",ROSMAP"]^{2}), 
		'half'=bquote(italic(h)[italic(e)*",GTEx"]^{2}==0.5%*%italic(h)[italic(e)*",ROSMAP"]^{2}), 
		'zero'=bquote(italic(h)[italic(e)*",ROSMAP"]^{2}*'; '*italic(h)[italic(e)*",GTEx"]^{2}==0),
		'overlap'=bquote(italic(h)[italic(e)*",GTEx"]^{2}==italic(h)[italic(e)*",ROSMAP"]^{2})
		)[[job_cat]]
	expr_facet_lab <- bquote(.(expr_facet_lab_var1)*'Expression Heritability ('*.(expr_facet_lab_var2)*')')
	return(expr_facet_lab)
}

################
# set directories
out_dir <- paste0(sim_dir, 'plot/')
setwd(out_dir)

# get jobs
load(paste0(sim_dir, 'power/data/jobs', suffix, '.Rdata'))

# phenotype
pheno_h2_list <- c(0.05, 0.1, 0.175, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)

filetypes <- c('jpeg')


subtitle_suf <- ''
subtitle_suf <- '_subtitle'

# load data
load(paste0(sim_dir, 'plot/plot_data_scenario5', suffix, '.Rdata'))


############################################
# set up palette
## https://pmassicotte.github.io/paletteer_gallery/
labels_vec <- c('PrediXcan-GTEx', 'TIGAR-ROSMAP', 'Naive', 'TIGAR-ROSMAP_valid', 'SR-TWAS', 'Avg-Base+SR')
breaks_vec <- labels_vec

palette <- 'ggsci::default_jama' # c('black','orange','blue','red','green','purple')
pal0 <- setNames(c(paletteer_d(palette, n=6)), c('black','orange','blue','red','green','purple'))
color_pal <- setNames(c(paletteer_d(palette, n=6)), c('Naive','TIGAR-ROSMAP','PrediXcan-GTEx','SR-TWAS','TIGAR-ROSMAP_valid', 'Avg-Base+SR'))[labels_vec]
shape_pal <- setNames(c(0, 2, 6, 3, 1, 5), c('PrediXcan-GTEx', 'TIGAR-ROSMAP', 'TIGAR-ROSMAP_valid', 'Naive', 'SR-TWAS', 'Avg-Base+SR'))[labels_vec]

legend_pos <- 'bottom'


####################
# POWER
####################
# subtitle
subtitle_var1 <- get_subtitle_var1(job_cat)
expr_facet_lab <- get_expr_facet_lab(job_cat)

dat1 <- na.omit(power_dat)

power_plot_raw <- ggplot(data=dat1, 
	aes(x=Hp2, y=power, color=cohort, shape=cohort)) + 
	geom_point(alpha=0.6, size=2) + 
	geom_line(alpha=0.9, size=0.5) + 
	theme_bw() +
	theme(legend.position=legend_pos, 
		legend.title=element_blank(),
		legend.margin=margin(0, 0, 0, 0)) +
	scale_color_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +
	scale_shape_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=shape_pal) +
	labs(
		title=NULL,
		subtitle=subtitle_var1,
		color=NULL,
		x=bquote('Phenotype Heritability ('*italic(h)[italic(p)]^{2}*')'),
		y='Power')

power_plot1 <- power_plot_raw + 
	facet_grid(ROSMAP_causal_prop ~ ROSMAP_He2)
power_plot1 <- facet_labs(power_plot1, 
		expr_facet_lab, 
		bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'))

power_plot2 <- power_plot_raw + 
	facet_grid(ROSMAP_He2 ~ ROSMAP_causal_prop, scales='free_y')
power_plot2 <- facet_labs(power_plot2,  
		bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'),
		expr_facet_lab)

for (filetype in filetypes) {
	ggsave(paste0(out_dir, 'power', subtitle_suf, suffix, '.', filetype), power_plot1, width=w, height=h)
	ggsave(paste0(out_dir, 'power_2', subtitle_suf, suffix, '.', filetype), power_plot2, width=w, height=h2)
}


#########################
## EXPRESSION PREDICTION TEST R2
#########################
boxplot_dat <- pred_dat

# subtitle
subtitle_var1 <- get_subtitle_var1(job_cat)
expr_facet_lab <- get_expr_facet_lab(job_cat)

boxplot_p_raw <- ggplot(data=boxplot_dat, aes(x=cohort, y=expr_R2, color=cohort)) + 
	geom_boxplot() +
	stat_summary(fun=mean, geom='point', position=position_dodge(.9), size=2.5, shape=18) +
	theme_bw() +
	theme(legend.position=legend_pos, 
		legend.title=element_blank(),
		legend.margin=margin(0, 0, 0, 0),
		legend.box='vertical',
		strip.text=element_text(
			color='black', 
			size=rel(0.9)
			),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title.x = element_text(color='white')
		) +
	scale_color_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +
	scale_fill_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +
	labs(
		subtitle=subtitle_var1,
		x=NULL,
		y=bquote('Expression Prediction'~italic(R)^2)
		)

boxplot_p1 <- boxplot_p_raw + 
	facet_grid(ROSMAP_causal_prop ~ ROSMAP_He2)
boxplot_p1 <- facet_labs(boxplot_p1, 
		expr_facet_lab, 
		bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'))

boxplot_p2 <- boxplot_p_raw + 
	facet_grid(ROSMAP_He2 ~ ROSMAP_causal_prop, scales='free_y')
boxplot_p2 <- facet_labs(boxplot_p2, 
		bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'),
		expr_facet_lab)

for (filetype in filetypes) {
	ggsave(paste0(out_dir, 'exprR2', subtitle_suf, suffix,'.', filetype), boxplot_p1, width=w, height=h)
	ggsave(paste0(out_dir, 'exprR2_2', subtitle_suf, suffix,'.', filetype), boxplot_p2, width=w, height=h2)
}

#########################
### CVR2 PLOTS
#########################
# can't do CVR2,R2 with Avg data because there's no evaluation step
cvr2_r2_x <- 'Naive'
cvr2_r2_y <- 'SR-TWAS'


# subtitle
subtitle_var1 <- get_subtitle_var1(job_cat)
expr_facet_lab <- get_expr_facet_lab(job_cat)

# data
cvr2_dat1 <- train_dat[(train_dat$cohort %in% c(cvr2_r2_x, cvr2_r2_y)), c('TargetID', 'ROSMAP_causal_prop', 'ROSMAP_He2', 'ValidCVR2', 'cohort')]

cvr2_dat2 <- dcast(cvr2_dat1, TargetID + ROSMAP_causal_prop +  ROSMAP_He2 ~ cohort, value.var='ValidCVR2')
CVR2_max <- max(cvr2_dat1$ValidCVR2, na.rm=TRUE)


# remove NA, values < 0
cvr2_dat2[is.na(cvr2_dat2)] <- 0
cvr2_dat2[which(cvr2_dat2[[cvr2_r2_x]] < 0), cvr2_r2_x] <- 0
cvr2_dat2[which(cvr2_dat2[[cvr2_r2_y]] < 0), cvr2_r2_y] <- 0

r2_cut <- 0.005

# set colors
cvr2_dat2$color_sort <- NA
cvr2_dat2[cvr2_dat2[[cvr2_r2_x]] <= r2_cut & cvr2_dat2[[cvr2_r2_y]] <= r2_cut, 'color_sort'] <- 4 #Neither
cvr2_dat2[cvr2_dat2[[cvr2_r2_x]] > r2_cut & cvr2_dat2[[cvr2_r2_y]] > r2_cut, 'color_sort'] <- 3 #Both
cvr2_dat2[cvr2_dat2[[cvr2_r2_x]] > r2_cut & cvr2_dat2[[cvr2_r2_y]] <= r2_cut, 'color_sort'] <- 2 #Other
cvr2_dat2[cvr2_dat2[[cvr2_r2_x]] <= r2_cut & cvr2_dat2[[cvr2_r2_y]] > r2_cut, 'color_sort'] <- 1 #SR
cvr2_dat2 <- cvr2_dat2[order(cvr2_dat2$color_sort, decreasing=TRUE),]

cvr2_dat2$color <- as.character(cvr2_dat2$color_sort)

# do plot
CVR2_plot_raw <- ggplot(cvr2_dat2, aes(x=.data[[cvr2_r2_x]], y=.data[[cvr2_r2_y]], color=color)) + 
	geom_point(alpha=0.5, size=1.875) + 
	geom_abline(intercept=0, slope=1, color='black', linetype='dashed', alpha=1) +
	xlim(0, CVR2_max) +
	ylim(0, CVR2_max) + 
	theme_bw() +
	theme(legend.pos=legend_pos,
		legend.margin=margin(0, 0, 0, 0)) +
	guides(colour = guide_legend(override.aes = list(size=3, alpha=0.9))) +
	scale_color_manual(name=bquote('CV'~italic(R)^2 > .(r2_cut)),
		labels=c(paste(cvr2_r2_y, 'Only'), paste(cvr2_r2_x, 'Only'), 'Both', 'Neither'),
		breaks=c('1', '2', '3', '4'),
		values=setNames(pal0[c('red','black','blue','orange')], c('1', '2', '3', '4'))) +
	labs(
		title=NULL,
		subtitle=subtitle_var1,
		x=cvr2_r2_x,
		y=cvr2_r2_y,
		fill='')

CVR2_plot1 <- CVR2_plot_raw + 
	facet_grid(ROSMAP_causal_prop ~ ROSMAP_He2)
CVR2_plot1 <- facet_labs(CVR2_plot1, 
	expr_facet_lab, 
	bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'))

CVR2_plot2 <- CVR2_plot_raw + 
	facet_grid(ROSMAP_He2 ~ ROSMAP_causal_prop, scales='free_y')
CVR2_plot2 <- facet_labs(CVR2_plot2, 
	bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'),
	expr_facet_lab)

for (filetype in filetypes) {
	ggsave(paste0(out_dir, 'CVR2_', cvr2_r2_x, '_', cvr2_r2_y, subtitle_suf, suffix, '.', filetype), CVR2_plot1, width=w, height=h)
	ggsave(paste0(out_dir, 'CVR2_2_', cvr2_r2_x, '_' ,cvr2_r2_y, subtitle_suf, suffix, '.', filetype), CVR2_plot2, width=w, height=h2)
}

#########################
## R2 PLOT
#########################
# subtitle
subtitle_var1 <- get_subtitle_var1(job_cat)
expr_facet_lab <- get_expr_facet_lab(job_cat)

# data
r2_dat1 <- train_dat[(train_dat$cohort %in% c(cvr2_r2_y, cvr2_r2_x)), c('TargetID', 'ROSMAP_causal_prop', 'ROSMAP_He2', 'ValidR2', 'cohort')]

r2_dat2 <- dcast(r2_dat1, TargetID + ROSMAP_causal_prop +  ROSMAP_He2 ~ cohort, value.var='ValidR2')
R2_max <- max(r2_dat1$ValidR2, na.rm=TRUE)

# remove NA, values < 0
r2_dat2[is.na(r2_dat2)] <- 0
r2_dat2[which(r2_dat2[[cvr2_r2_x]] < 0), cvr2_r2_x] <- 0
r2_dat2[which(r2_dat2[[cvr2_r2_y]] < 0), cvr2_r2_y] <- 0

r2_cut <- 0.005

# set colors
r2_dat2$color_sort <- NA
r2_dat2[r2_dat2[[cvr2_r2_x]] <= r2_cut & r2_dat2[[cvr2_r2_y]] <= r2_cut, 'color_sort'] <- 4 #Neither
r2_dat2[r2_dat2[[cvr2_r2_x]] > r2_cut & r2_dat2[[cvr2_r2_y]] > r2_cut, 'color_sort'] <- 3 #Both
r2_dat2[r2_dat2[[cvr2_r2_x]] > r2_cut & r2_dat2[[cvr2_r2_y]] <= r2_cut, 'color_sort'] <- 2 #Other
r2_dat2[r2_dat2[[cvr2_r2_x]] <= r2_cut & r2_dat2[[cvr2_r2_y]] > r2_cut, 'color_sort'] <- 1 #SR
r2_dat2 <- r2_dat2[order(r2_dat2$color_sort, decreasing=TRUE),]

r2_dat2$color <- as.character(r2_dat2$color_sort)

# do plot
R2_plot_raw <- ggplot(r2_dat2, aes(x=.data[[cvr2_r2_x]], y=.data[[cvr2_r2_y]], color=color)) + 
	geom_point(alpha=0.5, size=1.875) + 
	geom_abline(intercept=0, slope=1, color='black', linetype='dashed', alpha=1) +
	xlim(0, R2_max) +
	ylim(0, R2_max) +
	theme_bw() +
	theme(legend.pos=legend_pos,
		legend.margin=margin(0, 0, 0, 0)) +
	guides(colour = guide_legend(override.aes = list(size=3, alpha=0.9))) +
	scale_color_manual(name=bquote(italic(R)^2 > .(r2_cut)),
		labels=c(paste(cvr2_r2_y, 'Only'), paste(cvr2_r2_x, 'Only'), 'Both', 'Neither'),
		breaks=c('1', '2', '3', '4'),
		values=setNames(pal0[c('red','black','blue','orange')], c('1', '2', '3', '4'))) +
	labs(
		title=NULL,
		subtitle=subtitle_var1,
		x=cvr2_r2_x,
		y=cvr2_r2_y,
		fill='')

R2_plot1 <- R2_plot_raw + 
	facet_grid(ROSMAP_causal_prop ~ ROSMAP_He2)
R2_plot1 <- facet_labs(R2_plot1, 
	expr_facet_lab, 
	bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'))

R2_plot2 <- R2_plot_raw + 
	facet_grid(ROSMAP_He2 ~ ROSMAP_causal_prop, scales='free_y')
R2_plot2 <- facet_labs(R2_plot2, 
	bquote('Proportion of Causal SNPs ('*italic(p)['causal']*')'),
	expr_facet_lab)

for (filetype in filetypes) {
	ggsave(paste0(out_dir, 'R2_', cvr2_r2_x, '_', cvr2_r2_y, subtitle_suf, suffix, '.', filetype), R2_plot1, width=w, height=h)
	ggsave(paste0(out_dir, 'R2_2_', cvr2_r2_x, '_', cvr2_r2_y, subtitle_suf, suffix, '.', filetype), R2_plot2, width=w, height=h2)
}



