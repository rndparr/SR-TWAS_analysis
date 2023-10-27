#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

### load libraries
library(xtable)
library(ggplot2)
library(ggrepel)
# library(QTLRel)
library(lattice)
library(gridExtra)
library(GenomicRanges)
# library(extrafont)
library(reshape2)
library(ggbreak)


####
### TWAS to do
# twas <- 'AD_TWAS'
# twas <- 'PD_TWAS'
twas <- 'AD_TWAS2'

### directories
## HGCC
in_dir <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/', twas, '/')
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/050323/'
plot_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/050323/'
# out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'
# plot_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'

## laptop
# in_dir <- paste0('~/HGCC/YangFSSdata/SR_TWAS/', twas, '/')
# out_dir <- '~/HGCC/YangFSSdata/SR_TWAS/plots/'

### read in data
dat <- read.table(
	paste0(in_dir, 'all_', twas, '_results.txt'), 
	header=TRUE, sep='\t')

dat$GeneName <- trimws(as.character(dat$GeneName))

## SET UP MODELS
# dat <- dat[dat$dataset %in% c('PredDB','GTEx','ROSMAP', 'Naive3_ROSMAP','SR_ROSMAP','Naive_ROSMAP'),]
# dat <- dat[dat$dataset %in% c('PredDB','GTEx','ROSMAP', 'Naive3_ROSMAP','SR_ROSMAP'),]

# dat <- dat[dat$dataset %in% c('PredDB','GTEx','ROSMAP','ROSMAP2','Naive34basemodels_ROSMAP','SR4basemodels_ROSMAP'),]

# dat <- dat[dat$dataset %in% c('PredDB','GTEx','ROSMAP','ROSMAP2','Naive34basemodels_ROSMAP','SR4basemodels_ROSMAP','GTEx_LD_ROSMAP','ROSMAP_LD_ROSMAP','PredDB_LD_ROSMAP','ROSMAP2_LD_ROSMAP','SR4basemodels_ROSMAP_LD_ROSMAP','Naive34basemodels_ROSMAP_LD_ROSMAP','GTEx6tissues_GTEx_LD_GTEx','GTEx6tissues_Naive_GTEx_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP'),]


if (twas == 'PD_TWAS'){
	## PD-TWAS
	dat <- dat[dat$dataset %in% c('GTEx6tissues_GTEx_LD_GTEx', 
		# 'GTEx6tissues_Naive_GTEx_LD_GTEx',
		'GTEx_Brain_Substantia_nigra_LD_GTEx', 
		'GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 
		'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 
		'GTEx_Brain_Cortex_LD_GTEx', 
		'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx',
		'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 
		'GTEx_Whole_Blood_LD_GTEx', 'PredDB_GTEx_Brain_Substantia_nigra_LD_GTEx'),]

} else if (twas == 'AD_TWAS2'){
	## AD-TWAS2
	# ROSMAP was TIGAR DPR ROSMAP job
	# ROSMAP2 was prediXcan ROSMAP job
	dat <- dat[dat$dataset %in% c('GTEx','ROSMAP_LD_ROSMAP','PredDB','ROSMAP2_LD_ROSMAP','SR4basemodels_ROSMAP_LD_ROSMAP', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP','PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP'),]
}

model_order <- c('PredDB','GTEx','ROS','MAP','ROSMAP','ROSMAP2','Naive_MAP','Naive_ROS','Naive3_ROSMAP','SR_MAP','SR_ROS','SR_ROSMAP','Naive_ROS_MAP','SR_ROS_MAP','Naive34basemodels_ROSMAP','SR4basemodels_ROSMAP','GTEx_LD_ROSMAP','ROSMAP_LD_ROSMAP','PredDB_LD_ROSMAP','ROSMAP2_LD_ROSMAP','SR4basemodels_ROSMAP_LD_ROSMAP','Naive34basemodels_ROSMAP_LD_ROSMAP','GTEx6tissues_GTEx_LD_GTEx', 'GTEx6tissues_Naive_GTEx_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP','GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx','PredDB_GTEx_Brain_Substantia_nigra_LD_GTEx','PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP')

names_col <- c('PrediXcan_GTEx', 'PrediXcan_ROSMAP', paste0('TIGAR_', c('GTEx','ROS','MAP','ROSMAP')), 'Naive_MAP','Naive_ROS','Naive_ROSMAP','SR_MAP','SR_ROS','SR_ROSMAP','Naive_ROS_MAP','SR_ROS_MAP', 'Naive_ROSMAP', 'SR_ROSMAP','TIGAR_GTEx_LD_ROSMAP','TIGAR_ROSMAP_LD_ROSMAP','PrediXcan_GTEx_LD_ROSMAP','PrediXcan_ROSMAP_LD_ROSMAP','SR_LD_ROSMAP','Naive_LD_ROSMAP','SR_GTEx6_LD_GTEx', 'Naive_GTEx6_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx','ROSMAP_RNAseq_LD_ROSMAP','GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx','PredDB_GTEx_Brain_Substantia_nigra_LD_GTEx','PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP')
raw_names_col <- c('PredDB','ROSMAP2','GTEx','ROS','MAP','ROSMAP','Naive_MAP','Naive_ROS','Naive3_ROSMAP','SR_MAP','SR_ROS','SR_ROSMAP','Naive_ROS_MAP','SR_ROS_MAP','Naive34basemodels_ROSMAP','SR4basemodels_ROSMAP','GTEx_LD_ROSMAP','ROSMAP_LD_ROSMAP','PredDB_LD_ROSMAP','ROSMAP2_LD_ROSMAP','SR4basemodels_ROSMAP_LD_ROSMAP','Naive34basemodels_ROSMAP_LD_ROSMAP','GTEx6tissues_GTEx_LD_GTEx', 'GTEx6tissues_Naive_GTEx_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP','GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx','PredDB_GTEx_Brain_Substantia_nigra_LD_GTEx','PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP')
cbind(names_col, raw_names_col)

factor_labels_raw <- setNames(names_col, raw_names_col)


factor_levels <- model_order[which(model_order %in% unique(dat$dataset))]
factor_labels <- setNames(factor_labels_raw[factor_levels], NULL)
dat$dataset <- factor(dat$dataset, levels=factor_levels, labels=factor_labels)

## GET ZSCORE COLUMN
method <- 'SPred'
# method <- 'FUSION'
colnames(dat)[which(colnames(dat) == paste0(method,'_Z'))] <- 'Zscore'

## GET PVALUE
dat['Pvalue'] <- exp(pchisq(dat[['Zscore']]^2, 1, lower.tail=FALSE, log.p=TRUE))
dat['Pvalue_raw'] <- dat$Pvalue
dat <- dat[order(dat$dataset, dat$CHROM, dat$Pvalue), ]

## RENAME COLUMNS
colnames(dat)[which(colnames(dat) == 'CHROM')] <- 'Chrom'
colnames(dat)[which(colnames(dat) == 'GeneStart')] <- 'Start'
colnames(dat)[which(colnames(dat) == 'GeneEnd')] <- 'End'
colnames(dat)[which(colnames(dat) == 'GeneName')] <- 'Gene'

dat$Gene <- trimws(dat$Gene)

### function
# geom_label_size = 2; chr_vec = 1:22; ntop = 3;  sig_level = 2.5e-6;  size = 10;  chrGAP = 500;  colors = c('#44B5AD','#3A948E','#36807A','#2f615d');  present = FALSE; theme=theme_grey(); panel.border=element_blank(); sig_col1='#FD923F'; sig_col2='#D92B26'; sig_level_line_col='black'; sig_linetype='dashed'; sig_line_size=1; label_seg_col='black'; label_fill='white'; label_col='black'; segment_size=0.2; point_alpha=0.9; lab_force=2

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


# text = element_text(size = size), 
# axis.text.y=element_text(face='bold', 
# 	size=size+4),
# axis.text.x = element_text(face='bold', 
# 	size=size+4, angle=-90, vjust=0.5, hjust=0,),
# plot.tag = element_text(face='bold', size=size+10, family='Helvetica'),
# panel.grid.major.x = element_blank(),
# axis.title = element_text(face='bold', size=size+5),
# panel.grid.minor.x = element_blank(),
# # panel.border = element_blank(),
# plot.background = element_rect(fill=plot_bg_col)

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


# id_names <- c()
# id_names <- as.character(dat[dat$Pvalue <= 2.5e-6, 'Gene'])

# ## PD-TWAS
# id_names <- c("CD38", "GPNMB", "IDUA", "LA16c-385E7.1", "LRRC37A4P", "MMRN1", 
# "NDUFAF2", "RAB29", "SHROOM3", "SLC30A3", "VKORC1", "ZSWIM7")


output_manplots3 <- function(group1_models, group1_labels, group2_models, group2_labels, group1_name='base', group2_name='naivesr', suffix='', filetype='pdf', group1_combo_labels=NULL, group2_combo_labels=NULL, group1_plot_add=NULL, group2_plot_add=NULL, ...){

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

		# ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group1_name, suffix,'.', filetype), p1 + facet_grid(dataset ~ .) , width=5, height=0.7*base_height)

		ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group1_name, suffix,'.', filetype), p1 + facet_grid(dataset ~ .) , width=4.5, height=6.5)
	
		# ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group1_name, suffix,'_wrap.', filetype), p1 + facet_wrap(dataset ~ ., nrow=2, strip.position='right'), width=5, height=0.7*5)
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
## AD-TWAS2

output_manplots4 <- function(group1_models, group1_labels, group2_models, group2_labels, group1_name='base', group2_name='naivesr', suffix='', filetype='pdf', group1_combo_labels=NULL, group2_combo_labels=NULL, group1_plot_add=NULL, group2_plot_add=NULL, ...){

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

		# ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group1_name, suffix,'.', filetype), p1 + facet_grid(dataset ~ .) , width=5, height=0.7*base_height)

		ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group1_name, suffix,'.', filetype), p1 + facet_grid(dataset ~ .) , width=4.5, height=6.5)
	
		# ggsave(paste0(out_dir, twas, '_', method, '_manplots_', group1_name, suffix,'_wrap.', filetype), p1 + facet_wrap(dataset ~ ., nrow=2, strip.position='right'), width=5, height=0.7*5)
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
		}

		p2 <- p2 + scale_y_break(c(61,145))

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

id_names <- c('AC073842.1', 'ACE', 'CR1', 'DMPK', 'FAM13C', 'GFAP', 'HLA-DRA', 'PPP1R9B', 'SLC15A3', 'TREM2', 'ZSCAN26')


output_manplots4(
	group1_models=c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP'), 
	group1_labels=c('PrediXcan\nGTEx\nFCBA9', 'TIGAR\nGTEx\nFCBA9', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC'),
	group1_name='base',
	group2_models=c('PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP','ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP'), 
	group2_labels=c('PrediXcan\nROSMAP\nSMA','TIGAR\nROSMAP\nSMA','SR-TWAS\nROSMAP\nSMA'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='_SR4tissue',
	filetype='pdf')



dat[dat$dataset == 'PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP' & dat$Pvalue <= 1e-99, ]


dat[dat$Gene %in% c('GPR4','FOSB','BLOC1S3'),]


unique(dat$dataset)

min(dat[dat$dataset == 'PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP', 'Pvalue'])


min(dat[dat$dataset == 'SR_LD_ROSMAP', 'Pvalue'])
min(dat[dat$dataset == 'TIGAR_ROSMAP_LD_ROSMAP', 'Pvalue'])


colnames(dat)



output_manplots3(
	group1_models=c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP'), 
	group1_labels=c('PrediXcan\nGTEx\nFCBA9', 'TIGAR\nGTEx\nFCBA9', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC'),
	group1_name='base',
	group2_models=c('ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP'), 
	group2_labels=c('TIGAR\nROSMAP\nSMA','SR-TWAS\nROSMAP\nSMA'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='_SR4tissue',
	filetype='pdf')


# Brain_Frontal_Cortex_BA9	BRNCTXB

output_manplots3(
	group1_models=c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP'), 
	group1_labels=c('PrediXcan\nGTEx\nBRNCTXB', 'TIGAR\nGTEx\nBRNCTXB', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC'),
	group1_name='base',
	group2_models=c('ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP'), 
	group2_labels=c('TIGAR\nROSMAP\nSMA','SR-TWAS\nROSMAP\nSMA'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='_SR4tissue_abbrv',
	filetype='pdf')



###########
## PD-TWAS
id_names <- c('CD38', 'GPNMB', 'IDUA', 'LA16c-385E7.1', 'LRRC37A4P', 'MMRN1', 'NDUFAF2', 'RAB29', 'SHROOM3', 'SLC30A3', 'VKORC1', 'ZSWIM7')

# tissue abbreviations:
# 'BRNACC', 'GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx'
# 'BRNCDT', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx'
# 'BRNCTXA', 'GTEx_Brain_Cortex_LD_GTEx'
# 'BRNNNCC', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx'
# 'BRNPTM', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx'
# 'BLOOD', 'GTEx_Whole_Blood_LD_GTEx'

# BRNSNG Brain_Substantia_nigra

output_manplots3(
	group1_models=c('GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx'),
	group1_labels=c('BRNACC', 'BRNCDT', 'BRNCTXA', 'BRNNNCC', 'BRNPTM', 'BLOOD'),
	group1_combo_labels=c('TIGAR\nGTEx\nBRNACC', 'TIGAR\nGTEx\nBRNCDT', 'TIGAR\nGTEx\nBRNCTXA', 'TIGAR\nGTEx\nBRNNNCC', 'TIGAR\nGTEx\nBRNPTM', 'TIGAR\nGTEx\nBLOOD'),
	group1_name='base',
	group1_plot_add = scale_y_continuous(sec.axis = sec_axis(~ ., name='TIGAR GTEx Base Models', breaks=NULL, labels=NULL)),
	group2_models=c('PredDB_GTEx_Brain_Substantia_nigra_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx', 'SR_GTEx6_LD_GTEx'), 
	group2_labels=c('PrediXcan\nGTEx\nBRNSNG','TIGAR\nGTEx\nBRNSNG', 'SR-TWAS\nGTEx\nBRNSNG'),
	group2_combo_labels=c('PrediXcan\nGTEx\nBRNSNG','TIGAR\nGTEx\nBRNSNG', 'SR-TWAS\nGTEx\nBRNSNG'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='_GTEx6tissues_LD_GTEx_abbrv',
	filetype='pdf',
	theme=theme_bw())






output_manplots3(
	group1_models=c('GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx'), 
	group1_labels=c('Anterior\ncingulate\ncortex\nBA24', 'Caudate\nbasal\nganglia', 'Cortex', 'Nucleus\naccumbens\nbasal\nganglia', 'Putamen\nbasal\nganglia', 'Whole\nBlood'),
	group1_combo_labels=c('TIGAR\nGTEx\nBRNACC', 'TIGAR\nGTEx\nBRNCDT', 'TIGAR\nGTEx\nBRNCTXA', 'TIGAR\nGTEx\nBRNNNCC', 'TIGAR\nGTEx\nBRNPTM', 'TIGAR\nGTEx\nBLOOD'),
	group1_name='base',
	group1_plot_add = scale_y_continuous(sec.axis = sec_axis(~ ., name='TIGAR GTEx Base Models', breaks=NULL, labels=NULL)),
	group2_models=c('GTEx_Brain_Substantia_nigra_LD_GTEx', 'SR_GTEx6_LD_GTEx'), 
	group2_labels=c('TIGAR\nGTEx\nSubstantia nigra', 'SR-TWAS\nGTEx\nSubstantia nigra'),
	group2_combo_labels=c('TIGAR\nGTEx\nBRNSNG', 'SR-TWAS\nGTEx\nBRNSNG'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='_GTEx6tissues_LD_GTEx_allindsigids',
	filetype='pdf',
	theme=theme_bw())


output_manplots3(
	group1_models=c('GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx'), 
	group1_labels=c('BRNACC', 'BRNCDT', 'BRNCTXA', 'BRNNNCC', 'BRNPTM', 'BLOOD'),
	group1_combo_labels=c('TIGAR\nGTEx\nBRNACC', 'TIGAR\nGTEx\nBRNCDT', 'TIGAR\nGTEx\nBRNCTXA', 'TIGAR\nGTEx\nBRNNNCC', 'TIGAR\nGTEx\nBRNPTM', 'TIGAR\nGTEx\nBLOOD'),
	group1_name='base',
	group1_plot_add = scale_y_continuous(sec.axis = sec_axis(~ ., name='TIGAR GTEx Base Models', breaks=NULL, labels=NULL)),
	group2_models=c('GTEx_Brain_Substantia_nigra_LD_GTEx', 'SR_GTEx6_LD_GTEx'), 
	group2_labels=c('TIGAR\nGTEx\nBRNSNG', 'SR-TWAS\nGTEx\nBRNSNG'),
	group2_combo_labels=c('TIGAR\nGTEx\nBRNSNG', 'SR-TWAS\nGTEx\nBRNSNG'),
	group2_name='sr',
	group2_plot_add=NULL,
	suffix='_GTEx6tissues_LD_GTEx_allindsigids_abbrv',
	filetype='pdf',
	theme=theme_bw())

#####################

## COUNTS

## count nsig
nsig_dat <- cbind(sort(by(dat, dat$dataset, function(x) { sum(x$Pvalue < 2.5e-6, na.rm=TRUE)})))

## count n ind sig
most_sig_no_overlap2 <- function(data = dat){
	table_data <- data[data$Pvalue_raw <= 2.5e-6, c('Gene','Chrom','Start','End', 'Zscore', 'Pvalue')]
	table_data <- na.omit(table_data)

	gr_dat <- table_data[,c('Gene','Chrom','Start','End')]
	colnames(gr_dat) <- c('name','chr','start','end')
	gr <- makeGRangesFromDataFrame(gr_dat, keep.extra.columns=TRUE)
	gr_wide <- resize(gr, width = 2000000 + width(gr), fix = 'center')
	overlaps <- findOverlaps(gr_wide, type = 'any', select = 'first')
	# overlaps <- findOverlaps(gr, maxgap=1000000, type = 'any', select = 'first')
	most_sig_names <- as.data.frame(gr[!duplicated(overlaps)])[['name']]
	return(most_sig_names)
}
n_indsig_dat <- cbind(sort(by(dat, dat$dataset, function(x) {length(most_sig_no_overlap2(x)) } )))

## AD-TWAS2
n_indsig_dat_cols <- c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP','ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP')
n_indsig_dat_colsnames <- c('PrediXcan_GTEx_FCBA9', 'TIGAR_GTEx_FCBA9', 'PrediXcan_ROSMAP_DLPFC', 'TIGAR_ROSMAP_DLPFC','TIGAR_ROSMAP_SMA','SR-TWAS_ROSMAP_SMA')

# n_indsig_dat_cols <- c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP','ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP')
# n_indsig_dat_colsnames <- c('PrediXcan_GTEx_FCBA9', 'TIGAR_GTEx_FCBA9', 'PrediXcan_ROSMAP_DLPFC', 'TIGAR_ROSMAP_DLPFC','TIGAR_ROSMAP_SMA','SR-TWAS_SMA')


## PD-TWAS
n_indsig_dat_cols <- c('GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 
'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 
'GTEx_Brain_Cortex_LD_GTEx', 
'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 
'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 
'GTEx_Whole_Blood_LD_GTEx', 
'Naive_GTEx6_LD_GTEx', 
'GTEx_Brain_Substantia_nigra_LD_GTEx', 
'SR_GTEx6_LD_GTEx')

n_indsig_dat_colsnames <- c('TIGAR-GTEx-Brain_Anterior_cingulate_cortex_BA24', 
'TIGAR-GTEx-Brain_Caudate_basal_ganglia', 
'TIGAR-GTEx-Brain_Cortex', 
'TIGAR-GTEx-Brain_Nucleus_accumbens_basal_ganglia', 
'TIGAR-GTEx-Brain_Putamen_basal_ganglia', 
'TIGAR-GTEx-Whole_Blood', 
'Naive-GTEx-Brain_Substantia_nigra', 
'TIGAR-GTEx-Brain_Substantia_nigra', 
'SRTWAS-GTEx-Brain_Substantia_nigra')

####


count_table <- merge(nsig_dat, n_indsig_dat, by=0)
row.names(count_table) <- count_table$Row.names
colnames(count_table) <- c('Row.names','N_sig', 'N_ind_sig')

count_table <- count_table[n_indsig_dat_cols, c('N_sig', 'N_ind_sig')]
count_table$name <- n_indsig_dat_colsnames
count_table <- count_table[, c('name','N_sig', 'N_ind_sig')]
row.names(count_table) <- count_table$name
count_table <- count_table[, c('N_sig', 'N_ind_sig')]
count_table

##############


targetsz <- by(dat, dat$dataset, function(x) {sort(most_sig_no_overlap2(x)) } )
targetsz
targetsz[['SR-TWAS']]
dput(targetsz[['SR_GTEx6_LD_GTEx']])

# ####

# n_indsig_dat <- cbind(sort(by(dat, dat$dataset, function(x) {length(most_sig_no_overlap2(x)) } )))

# n_indsig_dat_cols <- base_labels

# bleh <- merge(nsig_dat, n_indsig_dat, by=0)
# row.names(bleh) <- bleh$Row.names
# colnames(bleh) <- c('Row.names','N_sig', 'N_ind_sig')
# bleh[n_indsig_dat_cols, c('N_sig', 'N_ind_sig')]



#####################

	# group1_models=c('GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx'), 
	# group1_labels=c('Anterior\ncingulate\ncortex\nBA24', 'Caudate\nbasal\nganglia', 'Cortex', 'Nucleus\naccumbens\nbasal\nganglia', 'Putamen\nbasal\nganglia', 'Whole\nBlood'),
	# group1_combo_labels=c('TIGAR\nGTEx\nAnterior\ncingulate\ncortex\nBA24', 'TIGAR\nGTEx\nCaudate\nbasal\nganglia', 'TIGAR\nGTEx\nCortex', 'TIGAR\nGTEx\nNucleus accumbens\nbasal\nganglia', 'TIGAR\nGTEx\nPutamen\nbasal\nganglia', 'TIGAR\nGTEx\nWhole\nBlood'),

# output_manplots3(
# 	group1_models=c('GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx'), 
# 	group1_labels=c('TIGAR\nGTEx\nAnterior cingulate cortex BA24', 'TIGAR\nGTEx\nCaudate basal ganglia', 'TIGAR\nGTEx\nCortex', 'TIGAR\nGTEx\nNucleus accumbens basal ganglia', 'TIGAR\nGTEx\nPutamen basal ganglia', 'TIGAR\nGTEx\nWhole Blood'),
# 	group1_name='base',
# 	group2_models=c('GTEx_Brain_Substantia_nigra_LD_GTEx', 'SR_GTEx6_LD_GTEx'), 
# 	group2_labels=c('TIGAR\nGTEx\nSubstantia nigra', 'SR-TWAS\nGTEx\nSubstantia nigra'),
# 	group2_name='sr',
# 	suffix='_GTEx6tissues_LD_GTEx_allindsigids',
# 	filetype='jpeg',
# 	theme=theme_bw())



#####################



# id_names <- c("CD38", "GPNMB",  "LRRC37A4P", "MMRN1", "RAB29", "SHROOM3", "VKORC1", "ZSWIM7","SLC30A3")
# # c("IDUA", "LA16c-385E7.1","NDUFAF2","SLC30A3")

# # SHROOM near old ind sig find: CCDC158
# # VKORC1 near old ind sig find: PRSS53
# # LRRC37A4P near old ind sig find: LINC02210 (which was near known TWAS risk gene LRRC37A2)
# # ZSWIM7 near old ind sig find: ADORA2B

# output_manplots3(
# 	group1_models=c('GTEx_Brain_Substantia_nigra_LD_GTEx', 'SR_GTEx6_LD_GTEx'), 
# 	group1_labels=c('GTEx\nSubstantia nigra', 'SR-TWAS\nSubstantia nigra'),
# 	group1_name='basesr',
# 	group2_models=c(), 
# 	group2_labels=c(),
# 	group2_name='',
# 	suffix='_GTEx6tissues_LD_GTEx',
# 	filetype='pdf',
# 	theme=theme_bw())



# # ASHG color slides
# output_manplots3(
# 	group1_models=c('GTEx_Brain_Substantia_nigra_LD_GTEx', 'SR_GTEx6_LD_GTEx'), 
# 	group1_labels=c('GTEx\nSubstantia nigra', 'SR-TWAS\nSubstantia nigra'),
# 	group1_name='basesr',
# 	group2_models=c(), 
# 	group2_labels=c(),
# 	group2_name='',
# 	suffix='_GTEx6tissues_LD_GTEx_ashg',
# 	filetype='pdf',
# 	point_colors=c('#1b683a','#62a323','#09597b'),
# 	theme=theme_bw() + theme(strip.text = element_text(color='white', face = 'bold'), 
# 		strip.background = element_rect(fill = '#09597b', colour = '#09597b')), 
# 	panel_border=element_rect(fill=NA,colour='#333333',size=0.175),
# 	sig_color1='#80c3e9', sig_color2='#dbc741', sig_line_size=0.5, point_alpha=0.75
# 	)




# ###########

# ## AD-TWAS2

# id_names <- c('AC073842.1', 'ACE', 'CR1', 'DMPK', 'FAM13C', 'GFAP', 'HLA-DRA', 'PPP1R9B', 'SLC15A3', 'TREM2', 'ZSCAN26')

# output_manplots3(
# 	group1_models=c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP'), 
# 	group1_labels=c('PrediXcan\nGTEx\nFCBA9', 'TIGAR\nGTEx\nFCBA9', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC'),
# 	group1_name='base',
# 	group2_models=c('ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP'), 
# 	group2_labels=c('TIGAR\nROSMAP\nSMA','SR-TWAS\nSMA'),
# 	group2_name='sr',
# 	suffix='_SR4tissue',
# 	filetype='pdf')


# output_manplots3(
# 	group1_models=c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP'), 
# 	group1_labels=c('PrediXcan\nGTEx\nFCBA9', 'TIGAR\nGTEx\nFCBA9', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC'),
# 	group1_name='base',
# 	group2_models=c('ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP'), 
# 	group2_labels=c('TIGAR\nROSMAP\nSMA','SR-TWAS\nSMA'),
# 	group2_name='sr',
# 	suffix='_SR4tissue_ashg',
# 	filetype='pdf',
# 	point_colors=c('#1b683a','#62a323','#09597b'),
# 	theme=theme_bw() + theme(strip.text = element_text(color='white', face = 'bold'),
# 		strip.background = element_rect(fill = '#09597b', colour = '#09597b')), 
# 	panel_border=element_rect(fill=NA,colour='#333333',size=0.175),
# 	sig_color1='#80c3e9', sig_color2='#dbc741', sig_line_size=0.5, point_alpha=0.75
# 	)


# # # ASHG color slides
# # output_manplots3(
# # 	group1_models=c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP','ROSMAP_RNAseq_LD_ROSMAP'), 
# # 	group1_labels=c('PrediXcan\nGTEx\nFCBA9', 'TIGAR\nGTEx\nFCBA9', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC','TIGAR\nROSMAP\nSMA'),
# # 	group1_name='base',
# # 	group2_models=c('SR_LD_ROSMAP'), 
# # 	group2_labels=c('SR-TWAS\nSMA'),
# # 	group2_name='sr',
# # 	suffix='_SR4tissue_ashg',
# # 	filetype='pdf',
# # 	point_colors=c('#1b683a','#62a323','#09597b'),
# # 	theme=theme_bw() + theme(strip.text = element_text(color='white', face = 'bold'),
# # 		strip.background = element_rect(fill = '#09597b', colour = '#09597b')), 
# # 	panel_border=element_rect(fill=NA,colour='#333333',size=0.175),
# # 	sig_color1='#80c3e9', sig_color2='#dbc741', sig_line_size=0.5, point_alpha=0.75
# # 	)



############
## COUNTS
dat <- temp

## count nsig
nsig_dat <- cbind(sort(by(dat, dat$dataset, function(x) { sum(x$Pvalue < 2.5e-6, na.rm=TRUE)})))

## count n ind sig
most_sig_no_overlap2 <- function(data = dat){
	table_data <- data[data$Pvalue_raw <= 2.5e-6, c('Gene','Chrom','Start','End', 'Zscore', 'Pvalue')]
	table_data <- na.omit(table_data)

	gr_dat <- table_data[,c('Gene','Chrom','Start','End')]
	colnames(gr_dat) <- c('name','chr','start','end')
	gr <- makeGRangesFromDataFrame(gr_dat, keep.extra.columns=TRUE)
	gr_wide <- resize(gr, width = 2000000 + width(gr), fix = 'center')
	overlaps <- findOverlaps(gr_wide, type = 'any', select = 'first')
	# overlaps <- findOverlaps(gr, maxgap=1000000, type = 'any', select = 'first')
	most_sig_names <- as.data.frame(gr[!duplicated(overlaps)])[['name']]
	return(most_sig_names)
}
n_indsig_dat <- cbind(sort(by(dat, dat$dataset, function(x) {length(most_sig_no_overlap2(x)) } )))


##############
# n_indsig_dat_cols <- base_labels


output_manplots3(
	group1_models=c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP'), 
	group1_labels=c('PrediXcan\nGTEx\nFCBA9', 'TIGAR\nGTEx\nFCBA9', 'PrediXcan\nROSMAP\nDLPFC', 'TIGAR\nROSMAP\nDLPFC'),
	group1_name='base',
	group2_models=c('ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP'), 
	group2_labels=c('TIGAR\nROSMAP\nSMA','SR-TWAS\nSMA'),
	group2_name='sr',
	suffix='_SR4tissue',
	filetype='pdf')


## AD-TWAS2
n_indsig_dat_cols <- c('PrediXcan_GTEx','TIGAR_GTEx', 'PrediXcan_ROSMAP_LD_ROSMAP', 'TIGAR_ROSMAP_LD_ROSMAP','ROSMAP_RNAseq_LD_ROSMAP','SR_LD_ROSMAP')
n_indsig_dat_colsnames <- c('PrediXcan_GTEx_FCBA9', 'TIGAR_GTEx_FCBA9', 'PrediXcan_ROSMAP_DLPFC', 'TIGAR_ROSMAP_DLPFC','TIGAR_ROSMAP_SMA','SR-TWAS_SMA')


bleh <- merge(nsig_dat, n_indsig_dat, by=0)
row.names(bleh) <- bleh$Row.names
colnames(bleh) <- c('Row.names','N_sig', 'N_ind_sig')

bleh <- bleh[n_indsig_dat_cols, c('N_sig', 'N_ind_sig')]
bleh$name <- n_indsig_dat_colsnames
bleh <- bleh[, c('name','N_sig', 'N_ind_sig')]
bleh

##############


targetsz <- by(dat, dat$dataset, function(x) {sort(most_sig_no_overlap2(x)) } )
targetsz
targetsz[['SR-TWAS']]
dput(targetsz[['SR_GTEx6_LD_GTEx']])
# PDTWAS most sig no overlap
# dat$dataset: Naive
#  [1] "RAB29"         "SLC30A3"       "MMRN1"         "CD38"         
#  [5] "IDUA"          "SHROOM3"       "NDUFAF2"       "GPNMB"        
#  [9] "VKORC1"        "LA16c-385E7.1" "LRRC37A4P"     "ZSWIM7"       
# ------------------------------------------------------------ 
# dat$dataset: SR-TWAS
#  [1] "RAB29"     "CCNT2-AS1" "MMRN1"     "CD38"      "GAK"       "NDUFAF2"  
#  [7] "GPNMB"     "IGSF9B"    "KAT8"      "SPPL2C"    "ZSWIM7"   

# length(most_sig_no_overlap2(x))


## 

print(dat[dat$dataset == 'GTEx_Brain_Substantia_nigra' & dat$Gene %in% targetsz[['GTEx_Brain_Substantia_nigra']], c(2:4,6,10,11)],row.names = FALSE)
print(dat[dat$dataset == 'ROSMAP' & dat$Gene %in% targetsz[['ROSMAP']], c(2:4,6,10,11)],row.names = FALSE)
print(dat[dat$dataset == 'SR-TWAS' & dat$Gene %in% targetsz[['SR-TWAS']], c(2:4,6,10,11)],row.names = FALSE)  





print(dat[dat$dataset == 'SR-TWAS' & dat$Gene %in% targetsz[['SR-TWAS']], ],row.names = FALSE)  


n_indsig_dat <- cbind(sort(by(dat, dat$dataset, function(x) {length(most_sig_no_overlap2(x)) } )))


# n_indsig_dat_cols <- c('PrediXcan_ROSMAP','PrediXcan_ROSMAP_LD_ROSMAP','PrediXcan_GTEx','PrediXcan_GTEx_LD_ROSMAP','TIGAR_ROSMAP','TIGAR_ROSMAP_LD_ROSMAP','TIGAR_GTEx','TIGAR_GTEx_LD_ROSMAP','Naive_ROSMAP','Naive_LD_ROSMAP','SR_ROSMAP','SR_LD_ROSMAP')

n_indsig_dat_cols <- base_labels

bleh <- merge(nsig_dat, n_indsig_dat, by=0)
row.names(bleh) <- bleh$Row.names
colnames(bleh) <- c('Row.names','N_sig', 'N_ind_sig')
bleh[n_indsig_dat_cols, c('N_sig', 'N_ind_sig')]




# paste0(c('PrediXcan_GTEx','PrediXcan_ROSMAP','TIGAR_GTEx','TIGAR_ROSMAP'), '_LD_ROSMAP')
# output_manplots2(, 
# , naive_sr_models, naive_sr_labels, suffix='_LD_GTEx', filetype='pdf', ...)



# max(abs(dat$Zscore), na.rm=TRUE)

# dat[which(dat$Zscore == Inf),]




# ###

# out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/RIP_121421/'
# # ASHG
# output_manplots('ROSMAP',
# 	colors=c('#1b683a', '#62a323', '#09597b'),
# 	theme=theme_bw() + theme(
# 		strip.text=element_text(color='white', face='bold'), 
# 		strip.background=element_rect(fill='#09597b', colour='#09597b')), 
# 	panel.border=element_rect(fill=NA, colour='#333333', size=0.175), 
# 	sig_col1='#80c3e9', sig_col2='#dbc741', 
# 	sig_linetype='dashed', sig_line_size=0.5, point_alpha=0.75,
# 	suffix='', filetype='pdf',
# 	point_padding=0,
# 	#lab_force=50, #seed=123456, 
# 	lab_force=50,
# 	max_iter=25000)



# ## PD-TWAS
# # includes all ind. sig genes and genes only Id'd by SR-TWAS
# id_names <- unique(c("CCDC158", "CCNT2", "CD38", "CNTN1", "CXCL9", "GAK", "GPNMB", 
# "MCCC1-AS1", "MMRN1", "NDUFAF2", "PRSS53", "PTP4A2P2", "RAB29", 
# "SNORA64", "STX4", "VKORC1", "NUCKS1", "KLHL7-DT", "CNTN1", "RNF40"
# ))

# ## genes Id'd by sr-twas
# id_names <- c("CCDC158", "CCNT2", "CD38", "CNTN1", "FAM13A", "GPNMB", "HSD3B7", 
# "KAT8", "KLHL7-DT", "MIR762HG", "MMRN1", "NUCKS1", "NUP42", "PHKG2", 
# "PRSS53", "PTP4A2P2", "RAB29", "RNF40", "STX4", "VKORC1", "ZNF689"
# )

# ## ind genes id'd by sr-twas, + unique
# id_names <- c("CCDC158", "CCNT2", "CD38", "CNTN1", "CNTN1", "GPNMB", "KLHL7-DT", "MMRN1", "NUCKS1", "PRSS53", "PTP4A2P2", "RAB29", "RNF40", "STX4",  "VKORC1")

# ## ind genes id'd by sr-twas, + unique minus genes id'd in all
# id_names <- c("CCDC158", "CCNT2", "CNTN1", "KLHL7-DT", "MMRN1", "NUCKS1", "PRSS53", "PTP4A2P2", "RAB29", "RNF40")
