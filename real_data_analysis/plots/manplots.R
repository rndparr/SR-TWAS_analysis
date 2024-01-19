#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

### load libraries
library(xtable)
library(ggplot2)
library(ggrepel)
library(lattice)
library(gridExtra)
library(GenomicRanges)
library(reshape2)
library(ggbreak)


####
### TWAS to do
# twas <- 'PD_TWAS'
# twas <- 'AD_TWAS2'

### directories
## HGCC
in_dir <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/', twas, '/')
load(file=paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/overlap_tables_data_', twas, '.RData'))
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/output/'

###############

### function
manhattan_plot <- function(data, # dataframe (columns: 'CHR','POS','Pvalue','label_text')
	## genes:
		point_colors=c('#44B5AD','#3A948E','#36807A','#2f615d'), # color for non-sig genes
		sig_color1='#FD923F', # color for sig. genes without labels
		sig_color2='#D92B26', # color for sig. genes with labels
		point_alpha=0.9, # transparency for genes
	## significance level:
		sig_level=2.5e-6, # significance level value
		sig_level_line_col='black', # line color
		sig_linetype='dashed', # linetype
		sig_line_size=1, # line size
	## plot theme, other plot variables:
		chr_vec=1:22, # chromosomes to plot
		chr_gap=500, # gap between chromosomes in x-axis
		theme=theme_bw(), # ggplot2 theme (can use custom theme)
		plot_bg_col=NULL, # background color if different from theme
		# panel_border=NULL,#element_blank(), # ggplot panel border (default: blank)
		panel_border=element_rect(fill=NA, colour='#333333', size=0.175),
		strip_background=element_rect(colour='black', size=0.175),
		# panel_border=element_rect(fill='#333333', colour='#333333', size=0.175),
		text_size=10, # text size
	## point labelling:
		geom_label_size=2, # label text size
		label_fill='white', # label background color
		label_col='black', # label border color
		label_seg_col='black', # color of line from label to point
		min_segment_length=0.01, # minimum length of line from label to point
		segment_size=0.2, # line from label to point
		label_force=2, # force of repulsion between overlapping text labels
		point_padding=1e-06, # padding around genes
		seed=NA, # optional seed for generating label positions
		max_iter=15000 # number of iterations to use to generate label positions
	){
	# setup dataframe for plotting; get plot positions from chromosome positions
	plot_dat <- NULL # dataframe
	endPos <- 0  # place on x-axis where last chromosome ended
	x_axis_chr_breaks <- NULL # chromosome label positions
	x_axis_chr_labels <- NULL # list of chromosomes to label
	for (chr in chr_vec) {
		# get data for chr
		temp <- data[data$CHR==chr, ]
		if (nrow(temp) > 0) {
			# append chromosome to list of chromosomes to label
			x_axis_chr_labels <- c(x_axis_chr_labels, chr)
			# get unique positions for this chr
			uniq_pos <- sort(unique(temp$POS))
			uniq_pos <- setNames(order(uniq_pos), uniq_pos)
			# set POS to order value of the unique positions
			temp$POS <- uniq_pos[as.character(temp$POS)]
			# get plot positions for genes on this chr
			temp$plotPos <- (temp$POS - min(temp$POS, na.rm=TRUE) ) + endPos + 1
			# get new end position based on max position
			endPos <- max(temp$plotPos, na.rm=TRUE) + chr_gap
			# append label position
			x_axis_chr_breaks <- c(x_axis_chr_breaks, mean(temp$plotPos, na.rm=TRUE) )
			# add rows to plot_dat
			plot_dat <- rbind(plot_dat, temp)
		}
	}
	# set min, max values for axes
	min_x <- min(plot_dat$plotPos)
	max_x <- max(plot_dat$plotPos)
	max_y <- max(-log10(plot_dat$Pvalue))
	# plot
	p <- ggplot(data=plot_dat, 
				aes(x=plotPos, y=-log10(Pvalue), label=label_text)) + 
		# non-sig. genes:
		geom_point(data=subset(plot_dat, Pvalue >= sig_level),
			aes(x=plotPos, y=-log10(Pvalue), color=factor(CHR)),
			size=1, alpha=point_alpha) + 
		scale_color_manual(values=rep(point_colors, 22)) +
		# sig. genes
		geom_point(data=subset(plot_dat, Pvalue < sig_level),
			aes(x=plotPos, y=-log10(Pvalue), fill=factor(CHR)),
			size=ifelse(subset(plot_dat, Pvalue < sig_level)$label_text=='', 1.25, 1.5),
			color=ifelse(subset(plot_dat, Pvalue < sig_level)$label_text=='', sig_color1, sig_color2),
			alpha=point_alpha) +
		# add labels
		geom_label_repel(data=subset(plot_dat, Pvalue < sig_level), 
			min.segment.length=min_segment_length,
			segment.size=segment_size,
			segment.color=label_seg_col,
			box.padding=1.1,
			size=geom_label_size, 
			alpha=1,
			ylim=c(-log10(sig_level), max_y),
			xlim=c(min_x, max_x),
			force=label_force,
			point.padding=point_padding,
			max.iter=max_iter,
			colour=label_col,
			fill=label_fill,
			seed=seed) +
		# significance level line
		geom_hline(yintercept=-log10(sig_level), 
			linetype=sig_linetype, 
			size=sig_line_size, 
			color=sig_level_line_col) +
		# remove legend
		guides(color='none', fill='none') + 
		# set axes titles
		labs(x='Chromosome', y=bquote(-"log"[10]("p-value"))) + 
		# x-axis labels, breaks
		scale_x_continuous(breaks=x_axis_chr_breaks, 
			labels=x_axis_chr_labels, 
			expand=c(0.01, 0)) + 
		# don't clip to extent of plot panel
		coord_cartesian(clip='off') +
		# pad y-axis
		scale_y_continuous(expand=c(0.05, 0)) +
		theme +
		# convenient theme options for a manhattan plot
		theme(
			text=element_text(size=text_size),
			# text=element_text(size=text_size, face='bold'), 
			# axis.title=element_text(face='bold', size=text_size),
			axis.text.x=element_text(size=text_size-1, 
				face='plain', angle=-90, 
				vjust=0.5, hjust=0),
			axis.text.y=element_text(#face='bold', 
					size=text_size-1),
			panel.grid.major.x=element_blank(),
			panel.grid.minor.x=element_blank(),
			panel.border=panel_border,
			plot.background=element_rect(fill=plot_bg_col),
			plot.tag=element_text(face='bold', size=text_size+5, family='Helvetica'),
			strip.background=strip_background
			)

}

# helper function to add labels from dataframe
do_man_plot <- function(dataset=NULL, data=dat, drop_most_sig=FALSE, id_names=NULL, extra_cols=NULL,...){
	if(!is.null(dataset)){
		plot_dat <- data[data[['dataset']] == dataset, c('Chrom','Start','Gene','Pvalue',extra_cols)]
		colnames(plot_dat) <- c('CHR','POS','ID','Pvalue',extra_cols)
	} else {
		plot_dat <- data[, c('Chrom','Start','Gene','Pvalue','dataset',extra_cols)]
		colnames(plot_dat) <- c('CHR','POS','ID','Pvalue','dataset',extra_cols)
	}

	# if drop most significant
	if (drop_most_sig){
		plot_dat <- plot_dat[-with(plot_dat, which(Pvalue == min(Pvalue))),]
	}

	if(!is.null(id_names)){
		plot_dat$label_text <- plot_dat$ID
		plot_dat[with(plot_dat, !(ID %in% id_names)), 'label_text'] <- ''
	} else {
		plot_dat$label_text <- ''
	}

	return(manhattan_plot(plot_dat, ...))
}


# function to automatically create and save plots with proper labels and sizing for base models, ensemble models, and all models
output_manplots <- function(group1_models, group1_labels, group2_models, group2_labels, group1_name='base', group2_name='naivesr', suffix='', filetype='pdf', group1_combo_labels=NULL, group2_combo_labels=NULL, group1_plot_add=NULL, group2_plot_add=NULL, ...){

	all_models <- c(group1_models, group2_models)

	if (is.null(group1_combo_labels)){
		group1_combo_labels <- group1_labels
	}

	if (is.null(group2_combo_labels)){
		group2_combo_labels <- group2_labels
	}

	all_labels <- c(group1_combo_labels, group2_combo_labels)

	if (length(all_models) == 5) {
		base_height <- 5
	} else if (length(all_models) == 6) {
		base_height <- 8
	} else if (length(all_models) == 7) {
		base_height <- 10
	} else if (length(all_models) == 8) {
		base_height <- 12
	} else {
		base_height <- 6.5
	}

	# Base models
	if (length(group1_models) > 0) {
		temp <- dat[dat$dataset %in% group1_models,]
		temp$dataset <- droplevels(temp$dataset)
		temp$dataset <- factor(temp$dataset, levels=group1_models, labels=group1_labels)
		p1 <- do_man_plot(data=temp, id_names=id_names, geom_label_size=1.5,...)
	
		if (!is.null(group1_plot_add)){
			p1 <- p1 + group1_plot_add
		}
		ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group1_name, suffix,'.', filetype), p1 + facet_grid(dataset ~ .) , width=4.5, height=6.5)

	}

	# Group2
	if (length(group2_models) > 0) {
		temp <- dat[dat$dataset %in% group2_models,]
		temp$dataset <- droplevels(temp$dataset)
		temp$dataset <- factor(temp$dataset, levels=group2_models, labels=group2_labels )
		p2 <- do_man_plot(data=temp, id_names=id_names, geom_label_size=1.5,...) + facet_grid(dataset ~ .)

		if (!is.null(group2_plot_add)){
			p2 <- p2 + group2_plot_add
		}

		if (length(group2_models) == 2) {
			group2_height <- 0.7*5	
		} else if (length(group2_models) == 3) {
			group2_height <- 0.7*7
		} else if (length(group2_models) == 4) {
			group2_height <- 0.7*9
		}
		ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group2_name, suffix,'.', filetype), p2, width=5, height=group2_height)
	}

	# all models
	if ((length(group2_models) > 0) & (length(group1_models) > 0)){
		temp <- dat[dat$dataset %in% all_models,]
		temp$dataset <- droplevels(temp$dataset)
		temp$dataset <- factor(temp$dataset, levels=all_models, labels=all_labels)
		p3 <- do_man_plot(data=temp, id_names=id_names, geom_label_size=1.5,...) + facet_grid(dataset ~ .)

		ggsave(paste0(out_dir, twas, '_', method, '_manplots', suffix, '_all.', filetype), p3 + facet_grid(dataset ~ .) , width=5, height=0.7*(base_height+6))
	}

}


###########
## AD-TWAS
id_names <- c('AC092849.1', 'AL110118.2', 'FOSB', 'HLA-DRA', 'RN7SL225P', 'SRD5A3P1')

output_manplots(
	group1_models=c('PrediXcan_GTEx', 'TIGAR_GTEx', 'PrediXcan_ROSMAP', 'TIGAR_ROSMAP'),
	group1_labels=c('PrediXcan\nGTEx\nBRNCTXB', 'TIGAR\nGTEx\nBRNCTXB', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC'),
	group1_name='base',
	group2_models=c('PrediXcan_valid','TIGAR_valid', 'SR', 'Avg_valid_SR'), 
	group2_labels=c('PrediXcan\nROSMAP\nSMA','TIGAR\nROSMAP\nSMA','SR-TWAS\nROSMAP\nSMA','Avg valid+SR\nROSMAP\nSMA'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='',
	filetype='pdf',
	theme=theme_bw())


###########
## PD-TWAS
id_names <- c('AC005082.12', 'ADORA2B', 'CD38', 'LA16c-431H6.7', 'MAPK8IP1P2', 'MMRN1', 'MYLPF', 'NDUFAF2', 'PARL')

# tissue abbreviations:
# 'BRNACC', 'GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx'
# 'BRNCDT', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx'
# 'BRNCTXA', 'GTEx_Brain_Cortex_LD_GTEx'
# 'BRNNNCC', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx'
# 'BRNPTM', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx'
# 'BLOOD', 'GTEx_Whole_Blood_LD_GTEx'

# BRNSNG Brain_Substantia_nigra

output_manplots(
	group1_models=c('BRNACC', 'BRNCDT', 'BRNCTXA', 'BRNNNCC', 'BRNPTM', 'BLOOD'),
	group1_labels=c('BRNACC', 'BRNCDT', 'BRNCTXA', 'BRNNNCC', 'BRNPTM', 'BLOOD'),
	group1_combo_labels=c('TIGAR\nGTEx\nBRNACC', 'TIGAR\nGTEx\nBRNCDT', 'TIGAR\nGTEx\nBRNCTXA', 'TIGAR\nGTEx\nBRNNNCC', 'TIGAR\nGTEx\nBRNPTM', 'TIGAR\nGTEx\nBLOOD'),
	group1_name='base',
	group1_plot_add = scale_y_continuous(sec.axis=sec_axis(~ ., name='TIGAR GTEx Base Models', breaks=NULL, labels=NULL)),
	group2_models=c('PrediXcan_valid','TIGAR_valid', 'SR', 'Avg_valid_SR'), 
	group2_labels=c('PrediXcan\nGTEx\nBRNSNG','TIGAR\nGTEx\nBRNSNG', 'SR-TWAS\nGTEx\nBRNSNG', 'Avg valid+SR\nGTEx\nBRNSNG'),
	group2_combo_labels=c('PrediXcan\nGTEx\nBRNSNG','TIGAR\nGTEx\nBRNSNG', 'SR-TWAS\nGTEx\nBRNSNG', 'Avg valid+SR\nGTEx\nBRNSNG'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='',
	filetype='pdf',
	theme=theme_bw())

