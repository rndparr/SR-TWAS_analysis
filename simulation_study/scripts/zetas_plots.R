#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=FALSE)

###############
# read in passed arguments
args <- (commandArgs(TRUE))
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==1) {
	sim_dir = as.character(args[[1]])
	setting_suffixes = c('', '_sameSNPsameHe', '_sameSNPhalfHe')
	setting_names = c('0.5_overlap', 'sameSNPsameHe', 'sameSNPhalfHe')

} else {
	sim_dir = as.character(args[[1]])
	setting_suffixes = as.character(args[[2]])
	setting_names = as.character(args[[3]])
	
	setting_suffixes <- strsplit(setting_suffixes, ',')[[1]]
	setting_names <- strsplit(setting_names, ',')[[1]]
}

n_settings <- length(setting_suffixes)

# load libraries
library(ggplot2)
library(reshape)
library(grid)
library(gtable)
library(gridExtra)
library(paletteer)
library(reshape2)

################
# set directories
out_dir <- paste0(sim_dir, 'plot/')
setwd(out_dir)

# output file path vars
subtitle_suf <- ''
filetypes <- c('png')

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
get_expr_facet_lab <- function(j) { 
	if (gtex_he2_factor == 'other') {
		expr_facet_lab_var1 <- 'ROSMAP '
		expr_facet_lab_var2 <- bquote(italic(h)[italic(e)*",ROSMAP"]^{2})
	} else if (gtex_he2_factor == 1) {
		expr_facet_lab_var1 <- ''
		expr_facet_lab_var2 <- bquote(italic(h)[italic(e)*",GTEx"]^{2}==italic(h)[italic(e)*",ROSMAP"]^{2})

	} else if (gtex_he2_factor == 0) {
		expr_facet_lab_var1 <- 'ROSMAP '
		expr_facet_lab_var2 <- bquote(italic(h)[italic(e)*",ROSMAP"]^{2}*'; '*italic(h)[italic(e)*",GTEx"]^{2}==0)

	} else {
		expr_facet_lab_var1 <- 'ROSMAP '
		expr_facet_lab_var2 <- bquote(italic(h)[italic(e)*",GTEx"]^{2}==.(gtex_he2_factor)%*%italic(h)[italic(e)*",ROSMAP"]^{2})
	}
	expr_facet_lab <- bquote(.(expr_facet_lab_var1)*'Expression Heritability ('*.(expr_facet_lab_var2)*')')
	return(expr_facet_lab)
}

get_pcausal_facet_lab <- function(){
	if (gtex_pc_factor == 1) {
		pcausal_facet_lab_var1 <- ''
		pcausal_facet_lab_var2 <- bquote(italic(p)['causal'])

	} else if(gtex_pc_factor == 0) {
		pcausal_facet_lab_var1 <- 'ROSMAP '
		pcausal_facet_lab_var2 <- bquote(italic(p)['ROSMAP']*'; '*italic(p)['GTEx']==0)

	} else {
		pcausal_facet_lab_var1 <- 'ROSMAP '
		pcausal_facet_lab_var2 <- bquote(italic(p)['GTEx']==.(gtex_pc_factor)%*%italic(p)['ROSMAP'])
	}
	pcausal_facet_lab <- bquote(.(pcausal_facet_lab_var1)*'Proportion of Causal SNPs ('*.(pcausal_facet_lab_var2)*')')
	return(pcausal_facet_lab)
}

get_subtitle <- function() {
	if (subtitle_suf == '') {
		return(NULL)
	} else {
		ret <- paste0('GTEx causal SNPs chosen for ',  sprintf("%.0f%%", gtex_poverlap * 100), ' overlap with ROSMAP causal SNPs.')
		return(ret)
	} 
}

#########################
# Plot size
w <- 8.5
h <- 8 + ifelse(subtitle_suf == '', 0, 0.5)
h2 <- 6.25 + ifelse(subtitle_suf == '', 0, 0.5)

################
## plot settings
gtex_pc_factor <- 1
gtex_he2_factor <- 'other'

subtitle_str <- get_subtitle()
expr_facet_lab <- get_expr_facet_lab()
pcausal_facet_lab <- get_pcausal_facet_lab()

# set up palette
palette <- 'ggsci::default_jama' # c('black','orange','blue','red','green','purple')
breaks_vec <- c('0.5_overlap', 'sameSNPsameHe', 'sameSNPhalfHe')

labels_vec <- c(
	bquote(atop(italic(h)[italic(e)*",GTEx"]^{2}==italic(h)[italic(e)*",ROSMAP"]^{2}*';'*phantom(000000)~phantom(000000)~phantom(000000)~phantom(000000),phantom(000000)~'50% causal SNP overlap')),
	bquote(atop(italic(h)[italic(e)*",GTEx"]^{2}==italic(h)[italic(e)*",ROSMAP"]^{2}*';'*phantom(000000)~phantom(000000)~phantom(000000)~phantom(000000),~'same causal SNPs')), 
	bquote(atop(italic(h)[italic(e)*",GTEx"]^{2}==0.5%*%italic(h)[italic(e)*",ROSMAP"]^{2}*';'*phantom(000000),~'same causal SNPs')))

color_pal <- setNames(paletteer_d(palette, n=3), c('0.5_overlap', 'sameSNPsameHe', 'sameSNPhalfHe'))

# other plot settings
legend_pos <- 'bottom'

####################
# LOAD DATA
####################

# read in data for all scenarios
for (i in 1:n_settings){
	load(paste0(sim_dir, 'plot/plot_data', setting_suffixes[[i]], '.Rdata'))
	train_dat$name <- setting_names[[i]]

	zetas_dat <- rbind(zetas_dat, train_dat[train_dat$cohort == 'SR-TWAS',])
}

# set up dataframe for plotting
zetas_dat <- zetas_dat[, c('Z0', 'Z1', 'ROSMAP_He2', 'ROSMAP_causal_prop', 'suffix')]
zetas_dat$name <- factor(zetas_dat$name)

rosmap_zetas_dat2 <- reshape::melt(zetas_dat, measure.vars='Z0', id.vars=c('suffix', 'ROSMAP_He2', 'ROSMAP_causal_prop'))
gtex_zetas_dat2 <- reshape::melt(zetas_dat, measure.vars='Z1', id.vars=c('suffix', 'ROSMAP_He2', 'ROSMAP_causal_prop'))


############
## ROSMAP Zetas
zetas_plot_raw <- ggplot(rosmap_zetas_dat2, aes(x=value, group=suffix, color=suffix, fill=suffix)) +
	geom_histogram(alpha=0.15, 
		position='identity', 
		bins=20) +
	theme_bw() +
	theme(legend.position=legend_pos, 
		legend.title=element_blank(),
		legend.margin=margin(0, 0, 0, 0)) +
	scale_color_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +
	scale_fill_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +	
	labs(
		title=NULL,
		subtitle=subtitle_str,
		color=NULL,
		x=bquote('TIGAR-ROSMAP Zeta ('*zeta*')'),
		y='count')

zetas_plot1 <- zetas_plot_raw + 
	facet_grid(ROSMAP_causal_prop ~ ROSMAP_He2)
zetas_plot1 <- facet_labs(zetas_plot1, 
		expr_facet_lab, pcausal_facet_lab)

zetas_plot2 <- zetas_plot_raw + 
	facet_grid(ROSMAP_He2 ~ ROSMAP_causal_prop)
zetas_plot2 <- facet_labs(zetas_plot2,  
		pcausal_facet_lab, expr_facet_lab)

for (filetype in filetypes) {
	ggsave(paste0(out_dir, 'rosmap_zetas', subtitle_suf, suffix, '.', filetype), zetas_plot1, width=w, height=h)
	ggsave(paste0(out_dir, 'rosmap_zetas_2', subtitle_suf, suffix, '.', filetype), zetas_plot2, width=w, height=h2)
}


## GTEx Zetas
zetas_plot_raw <- ggplot(gtex_zetas_dat2, aes(x=value, group=suffix, color=suffix, fill=suffix)) +
	geom_histogram(alpha=0.15, 
		position='identity', 
		bins=20) +
	theme_bw() +
	theme(legend.position=legend_pos, 
		legend.title=element_blank(),
		legend.margin=margin(0, 0, 0, 0)) +
	scale_color_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +
	scale_fill_manual(name='',
		labels=labels_vec,
		breaks=breaks_vec,
		values=color_pal) +	
	labs(
		title=NULL,
		subtitle=subtitle_str,
		color=NULL,
		x=bquote('PrediXcan-GTEx Zeta ('*zeta*')'),
		y='count')

zetas_plot1 <- zetas_plot_raw + 
	facet_grid(ROSMAP_causal_prop ~ ROSMAP_He2)
zetas_plot1 <- facet_labs(zetas_plot1, 
		expr_facet_lab, pcausal_facet_lab)

zetas_plot2 <- zetas_plot_raw + 
	facet_grid(ROSMAP_He2 ~ ROSMAP_causal_prop)
zetas_plot2 <- facet_labs(zetas_plot2,  
		pcausal_facet_lab, expr_facet_lab)

for (filetype in filetypes) {
	ggsave(paste0(out_dir, 'gtex_zetas', subtitle_suf, suffix, '.', filetype), zetas_plot1, width=w, height=h)
	ggsave(paste0(out_dir, 'gtex_zetas_2', subtitle_suf, suffix, '.', filetype), zetas_plot2, width=w, height=h2)
}
