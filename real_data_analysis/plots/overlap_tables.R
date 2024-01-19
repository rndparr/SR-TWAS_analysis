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

twas <- 'PD_TWAS'
# twas <- 'AD_TWAS2'

Mb <- 1e+06

# directories
in_dir <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/', twas, '/')
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'
pmr_dir <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/PMR/', twas, '/')
pmr_result_dir <- paste0(pmr_dir, 'results/LD_0.99_lambda_0/')

# load functions
source(paste0(out_dir, 'functions.R'))
	# format_pval(x)
	# wideScreen()
	# strsplits(x, splits, ...)
	# near_gwas_snp(row, assoc_locs=assoc_locs)
	# near_gwas_gene(row, gwas_data=prev_id_gwas_dat)
	# near_gwas(row, gwas_data=prev_id_gwas_dat)
	# near_twas(row, twas_data=prev_id_dat)
	# most_sig_no_overlap(dataset, data=dat)
	# get_overlap_locus_group(data=dat)
	# get_overlap_locus_group_rank(data=dat)
	# annot_name(row)
	# annot_letter(row)

#############################

if (twas == 'PD_TWAS'){
	## PD-TWAS
	factor_levels <- c( 
		'GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 
		'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 
		'GTEx_Brain_Cortex_LD_GTEx', 
		'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx',
		'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 
		'GTEx_Whole_Blood_LD_GTEx',
		'PredDB_GTEx_Brain_Substantia_nigra_LD_GTEx',		
		'GTEx_Brain_Substantia_nigra_LD_GTEx',
		'GTEx6tissues_GTEx_LD_GTEx', 
		'GTEx_BRNSNG_valid_dat_models_all3')

	factor_labels <- c('BRNACC', 'BRNCDT', 'BRNCTXA', 'BRNNNCC', 'BRNPTM', 'BLOOD', 'PrediXcan_valid', 'TIGAR_valid', 'SR', 'Avg_valid_SR')

	SR_dataset_name <- 'Avg_valid_SR'
	other_dataset_names <- c('SR', 'PrediXcan_valid', 'TIGAR_valid')

} else if (twas == 'AD_TWAS2'){
	## AD-TWAS2
	# ROSMAP was TIGAR DPR ROSMAP job
	# ROSMAP2 was prediXcan ROSMAP job
	factor_levels <- c(
		'PredDB', 'ROSMAP2_LD_ROSMAP', # PrediXcan base models
		'GTEx', 'ROSMAP_LD_ROSMAP', # TIGAR base models
		'PrediXcan_ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP', # PrediXcan_valid
		'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP', # TIGAR_valid
		'SR4basemodels_ROSMAP_LD_ROSMAP',
		'ROSMAP_SMA_valid_dat_models_all3')

	factor_labels <- c(
		'PrediXcan_GTEx', 'PrediXcan_ROSMAP',
		'TIGAR_GTEx', 'TIGAR_ROSMAP', 'PrediXcan_valid', 'TIGAR_valid', 'SR', 'Avg_valid_SR')

	SR_dataset_name <- 'Avg_valid_SR'
	other_dataset_names <- c('SR', 'PrediXcan_valid', 'TIGAR_valid')
}
cbind(factor_levels, factor_labels)

# read in data
dat <- read.table(
	paste0(in_dir, 'all_', twas, '_results.txt'), 
	header=TRUE, sep='\t')
dat$GeneName <- trimws(as.character(dat$GeneName))
dat <- dat[dat$dataset %in% factor_levels,]
dat$dataset <- factor(dat$dataset, levels=factor_levels, labels=factor_labels)

# GET ZSCORE COLUMN
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
dat$TargetIDold <- dat$TargetID
dat$TargetID <- do.call(rbind, strsplit(dat$TargetID, '.', fixed=TRUE))[,1]

dat_raw <- dat

###############
## Output data to conduct PMR

# # gene_info <-  na.omit(dat[dat$Pvalue <= 2.5e-6 & dat$dataset=='Avg_valid_SR', c('Chrom', 'Start', 'End', 'TargetIDold', 'Gene', 'n_snps', 'Zscore')])
# # rownames(gene_info) <- NULL
# # colnames(gene_info) <- c('Chrom', 'Start', 'End', 'TargetID', 'Gene', 'n_snps', 'Zscore')

# ## compare with SR sig genes
# gene_info1 <-  na.omit(dat[dat$Pvalue <= 2.5e-6 & dat$dataset=='Avg_valid_SR', c('Chrom', 'Start', 'End', 'TargetIDold', 'Gene', 'n_snps', 'Zscore')])
# rownames(gene_info1) <- NULL
# colnames(gene_info1) <- c('Chrom', 'Start', 'End', 'TargetID', 'Gene', 'n_snps', 'Zscore')
# nrow(gene_info1)

# gene_info <-  na.omit(dat[dat$Pvalue <= 2.5e-6 & dat$dataset=='SR', c('Chrom', 'Start', 'End', 'TargetIDold', 'Gene', 'n_snps', 'Zscore')])
# rownames(gene_info) <- NULL
# colnames(gene_info) <- c('Chrom', 'Start', 'End', 'TargetID', 'Gene', 'n_snps', 'Zscore')
# nrow(gene_info)

# gene_info <- gene_info[gene_info$Gene %in% setdiff(gene_info$Gene, gene_info1$Gene), ]
# rownames(gene_info) <- NULL
# nrow(gene_info)

# save(gene_info, file=paste0(pmr_dir, 'gene_info.RData'))
# write.table(gene_info$TargetID, paste0(pmr_dir, 'gene_list.txt'), row.names=FALSE, sep='\t', quote=FALSE, col.names=FALSE)

###############
## GWAS CATALOG DATA
assoc_locs_genes <- read.csv(
	paste0(out_dir, list('AD_TWAS'='AD_TWAS', 'PD_TWAS'='PD_TWAS', 'AD_TWAS2'='AD_TWAS')[[twas]], '_GWAS_catalog.csv'))
assoc_locs_genes <- assoc_locs_genes[order(assoc_locs_genes$Location), ]
rownames(assoc_locs_genes) <- NULL


# by SNP location
assoc_locs <- as.data.frame(do.call(rbind, strsplit(
	unlist(strsplit(assoc_locs_genes[!(assoc_locs_genes$Location %in% c('', 'Mapping not available')), 'Location'], '|', fixed=TRUE)), 
	':', 
	fixed=TRUE)))
colnames(assoc_locs) <- c('Chrom', 'loc')

assoc_locs <- assoc_locs[assoc_locs[['Chrom']] != 'X', ]
assoc_locs$Chrom <- as.integer(assoc_locs$Chrom)
assoc_locs$loc <- as.integer(assoc_locs$loc)

# by gene
assoc_genes <- strsplits(unique(assoc_locs_genes$Mapped.gene), c('[[:space:]]', ',', ';'))

# list of genes with locations
dat$prev_id_gwas <- as.integer(dat$Gene %in% assoc_genes)
prev_id_gwas_dat <- dat[dat$prev_id_gwas > 0, c('Chrom', 'Start', 'End', 'prev_id_gwas')]
prev_id_gwas_dat <- prev_id_gwas_dat[order(prev_id_gwas_dat$Chrom, prev_id_gwas_dat$Start), ]
prev_id_gwas_dat <- prev_id_gwas_dat[!duplicated(prev_id_gwas_dat), ]


###############
## TWAS data
load(paste0(out_dir, 'prev_id_twas_', twas, '.RData'))
prev_id_dat$Gene <- trimws(as.character(prev_id_dat$Gene))

##############################
## data set up
# windows for each gene
dat$w_start <- as.integer(dat$Start - Mb)
dat$w_end <- as.integer(dat$End + Mb)

# significant
for(dataset in levels(dat$dataset)){
	sig_str <- paste0('sig_', dataset)
	dat[[sig_str]] <- as.integer(dat$Gene %in% dat[dat$dataset == dataset & dat$Pvalue <= 2.5e-6, 'Gene'])
}

# most sig in region
for(dataset in levels(dat$dataset)){
	most_sig_str <- paste0('most_sig_', dataset)
	dat[[most_sig_str]] <- as.integer(dat$Gene %in% most_sig_no_overlap(dataset, data=dat))
}

dat[['near_prev_id_gwas_snp']] <- apply(dat, 1, near_gwas_snp, gwas_assoc_locs=assoc_locs)
dat[['near_prev_id_gwas_gene']] <- apply(dat, 1, near_gwas_gene, gwas_data=prev_id_gwas_dat)
dat[['near_prev_id_gwas']] <- ifelse(dat$near_prev_id_gwas_snp > 0 | dat$near_prev_id_gwas_gene > 0, 1, 0)


###############
## SIG_DAT DATAFRAME FOR MAKING TABLES; updated with new dat columns
sig_dat <- na.omit(dat[dat$Pvalue <= 2.5e-6, ])
rownames(sig_dat) <- NULL

sig_dat$most_sig_n <- rowSums(sig_dat[,colnames(sig_dat)[which(startsWith(colnames(sig_dat), 'most_sig'))]])
sig_dat$most_sig_any <- as.integer(rowSums(sig_dat[,colnames(sig_dat)[which(startsWith(colnames(sig_dat), 'most_sig'))]]) > 0)

sig_dat$prev_id_twas <- as.integer(sig_dat$TargetID %in% prev_id_dat$TargetID)
sig_dat$near_prev_id_twas <- apply(sig_dat, 1, near_twas, twas_data=prev_id_dat)
sig_dat$GeneName <- apply(sig_dat, 1, annot_name)

################################################################################################
pmr_results <- read.table(paste0(pmr_result_dir, 'results_all.txt'), sep='\t', header=TRUE)
pmr_results$pmr_causal_pvalue <- as.numeric(pmr_results$causal_pvalue)
pmr_results$pmr_causal_pval_bonf <- p.adjust(pmr_results$pmr_causal_pvalue, 'bonferroni')
pmr_results$pmr_sig <- ifelse(pmr_results$pmr_causal_pval_bonf < 0.05, 1, 0)

pmr_results$pmr_causal_pval_fdr <- p.adjust(pmr_results$pmr_causal_pvalue, 'fdr')
pmr_results$pmr_fdr_sig <- ifelse(pmr_results$pmr_causal_pval_fdr < 0.05, 1, 0)
save(pmr_results, file=paste0(pmr_dir, 'pmr_results.RData'))

sig_dat <- merge(sig_dat, pmr_results[, c('i','causal_effect', 'pleiotropy_effect', 'pleiotropy_pvalue', 'sigma_cisSNP', 'sigma_error_1', 'sigma_error_2','TargetID', 'pmr_causal_pvalue', 'pmr_causal_pval_bonf', 'pmr_sig')], by.x='TargetIDold', by.y='TargetID', all.x=TRUE, all.y=TRUE)


save(list=ls(all.names=TRUE), file=paste0(out_dir, 'overlap_tables_data_', twas, '.RData'))

####################################################################################################################

print_res_info <- function(res_dat_annots, res_dat, print_genes=TRUE, prop_tab=FALSE, prop_merge_tab=FALSE, ...){
	## setup for count results tables
	# name_convert <- c('a'='GWAS', 'b'='near GWAS', 'ad'='GWAS, near TWAS', 'bd'='near GWAS, near TWAS', 'ac'='GWAS, TWAS', 'd'='near TWAS', 'c'='TWAS', 'bc'=''near GWAS, TWAS')
	names_order <- c('0', 'a', 'b', 'ad', 'bd', 'ac', 'd', 'c', 'bc')
	name_convert <- c('0'='novel','a'='GWAS', 'b'='~GWAS', 'ad'='GWAS,~TWAS', 'bd'='~GWAS,~TWAS', 'ac'='GWAS,TWAS', 'd'='~TWAS', 'c'='TWAS', 'bc'='~GWAS,TWAS')
	name_convert_long <- c('0'='novel:\t\t\t','a'='GWAS:\t\t\t', 'b'='near GWAS:\t\t', 'ad'='GWAS, near TWAS:\t', 'bd'='near GWAS, near TWAS:\t', 'ac'='GWAS, TWAS:\t\t', 'd'='near TWAS:\t\t', 'c'='TWAS:\t\t\t', 'bc'='near GWAS, TWAS:\t')
	name_convert <- name_convert[names_order]

	format_res_row <- function(search) {
		x <- paste(sort(res_dat[names(res_dat_annots[res_dat_annots==search]), 'Gene']), collapse=', ', sep='')
		return(paste0(name_convert_long[search], ifelse(x != '', x, '-'), '\n'))
	}

	annot_tab_dat <- data.frame('str'=res_dat_annots)
	annot_tab_dat$gwas <- 0
	annot_tab_dat[grepl('a', annot_tab_dat$str), 'gwas'] <- 'GWAS'
	annot_tab_dat[grepl('b', annot_tab_dat$str), 'gwas'] <- 'near_GWAS'
	annot_tab_dat$gwas <- factor(annot_tab_dat$gwas, levels=c('GWAS','near_GWAS', 0), labels=c('GWAS','near_GWAS', 'NA'))
	annot_tab_dat$twas <- 0
	annot_tab_dat[grepl('c', annot_tab_dat$str), 'twas'] <- 'TWAS'
	annot_tab_dat[grepl('d', annot_tab_dat$str), 'twas'] <- 'near_TWAS'
	annot_tab_dat$twas <- factor(annot_tab_dat$twas, levels=c('TWAS','near_TWAS', 0), labels=c('TWAS','near_TWAS','NA'))
	# addmargins(table(annot_tab_dat$gwas, annot_tab_dat$twas))

	## prep formatted
	res_vec <- table(res_dat_annots)
	not_in <- setdiff(names(name_convert), names(res_vec))
	res_vec_letters <- c(setNames(rep(0, length(not_in)), not_in), res_vec)[names(name_convert)]
	res_vec <- res_vec_letters
	# res_vec
	names(res_vec) <- name_convert[names(res_vec)]


	tab_dat <- table(annot_tab_dat$gwas, annot_tab_dat$twas)

	cat('-------------------------------------------------------')
	if (prop_tab){
		print(round(addmargins(prop.table(table(annot_tab_dat$gwas, annot_tab_dat$twas))), 2))
	} else if (prop_merge_tab) {
		#
		ntab <- as.data.frame.matrix(addmargins(tab_dat))
		ptab <- as.data.frame.matrix(addmargins(prop.table(tab_dat)), check.names=FALSE)

		tab <- ntab
		for (tcol in c('TWAS', 'near_TWAS', 'NA', 'Sum')){
			for (grow in c('GWAS', 'near_GWAS', 'NA', 'Sum')){
			tab[grow,tcol] <- paste0(ntab[grow,tcol], ' (', sprintf("%.2f",ptab[grow,tcol]), ')')
			}
		}
		colnames(tab) <- c('TWAS', 'near_TWAS', 'na_TWAS', 'Total')
		rownames(tab) <- c('GWAS', 'near_GWAS', 'na_GWAS', 'Total')

		# if (all(ntab[,'NA']==0) & all(ntab['NA',]==0)) { 
		# 	tab <- tab[c('GWAS', 'near_GWAS', 'Total'), c('TWAS', 'near_TWAS', 'Total')]
		# }
		if (all(tab_dat[,'NA']==0)) { 
			tab <- tab[, c('TWAS', 'near_TWAS', 'Total')]
		}
		if (all(tab_dat['NA',]==0)) { 
			tab <- tab[c('GWAS', 'near_GWAS', 'Total'), ]
		}

		cat('\n')
		print(tab)

	} else {
		print(addmargins(table(annot_tab_dat$gwas, annot_tab_dat$twas)))
	}

	if(print_genes){
		cat('-------------------------\n')
		cat(sapply(names_order, format_res_row, USE.NAMES=FALSE), sep='')
	}
	cat('-------------------------------------------------------')
}

# table_data_sr <- res_dat

print_res_table <- function(table_data_sr, start_end_cols=FALSE, print_locus=FALSE, ...){
	### merge other data

	table_data_raw <- table_data_sr

	n_other <- length(other_dataset_names)
	for (i in 1:n_other) {
		other_dataset_name <- other_dataset_names[i]
		table_data_other <- sig_dat[sig_dat$dataset == other_dataset_name , c('GeneName', 'Zscore', 'Pvalue')]
		# colnames(table_data_other) <- c('GeneName', paste0('Zscore_', LETTERS[i]), paste0('Pvalue_', LETTERS[i]))
			colnames(table_data_other) <- c('GeneName', paste0('Zscore_', other_dataset_name), paste0('Pvalue_', other_dataset_name))

		table_data_raw <- merge(table_data_raw, table_data_other, all.x=TRUE, by='GeneName')
		table_data_raw <- table_data_raw[order(table_data_raw$Chrom, table_data_raw$Start, table_data_raw$Pvalue), ]
	}


	############
	## FORMAT TABLE
	# pmr_causal_pval_bonf
	other_cols <- c()
	for (i in 1:n_other) {
		# other_cols <- c(other_cols, paste0('Zscore_', LETTERS[i]), paste0('Pvalue_', LETTERS[i]))
		other_dataset_name <- other_dataset_names[i]
		other_cols <- c(other_cols, paste0('Zscore_', other_dataset_name), paste0('Pvalue_', other_dataset_name))
	}

	if (start_end_cols) {
		table_data <- table_data_raw[, c('GeneName', 'Chrom', 'Start', 'End', 'pmr_causal_pval_bonf', 'Zscore', 'Pvalue', other_cols)]
		colnames(table_data) <- c('Gene', 'CHR', 'Start', 'End', 'PMRCausalPvalue','Zscore', 'Pvalue', other_cols)
		display_vec <- c('s','s','d','d','d','e','f','e', rep(c('f','e'), n_other))
	} else if (print_locus) {
		table_data <- table_data_raw[, c('GeneName', 'Chrom', 'Start', 'End', 'pmr_causal_pval_bonf', 'Zscore', 'Pvalue', other_cols)]
		table_data$locus <- paste(table_data$Chrom, ':', table_data$Start, '-', table_data$End, sep='')

		table_data <- table_data[, c('GeneName', 'locus', 'pmr_causal_pval_bonf', 'Zscore', 'Pvalue', other_cols)]
		colnames(table_data) <- c('Gene', 'Locus', 'PMRCausalPvalue','Zscore', 'Pvalue', other_cols)
		display_vec <- c('s','s','s','e','f','e', rep(c('f','e'), n_other))	

	} else {
		table_data <- table_data_raw[, c('GeneName', 'Chrom', 'pmr_causal_pval_bonf', 'Zscore', 'Pvalue', other_cols)]
		colnames(table_data) <- c('Gene', 'CHR', 'PMRCausalPvalue','Zscore', 'Pvalue', other_cols)
		display_vec <- c('s','s','d','e','f','e', rep(c('f','e'), n_other))
	}


	xtable_data <- xtable(table_data,
		math.style.exponent=TRUE, booktabs=TRUE, display=display_vec)

	print(xtable_data, booktabs=TRUE, include.rownames=FALSE, 
		sanitize.text.function=identity, 
		table.placement='H', comment=FALSE,
		NA.string='-')
}


print_all_res <- function(print_all_sig_table=FALSE, ...){
	# all sig
	# non-ind sig table
	cat(paste('* ', twas, '\n# all sig genes:\n'))
	res_dat_annots <- apply(sig_dat[sig_dat[['dataset']] == SR_dataset_name, ], 1, annot_letter)
	res_dat <- sig_dat[sig_dat[['dataset']] == SR_dataset_name, ]

	# info about annotations
	print_res_info(res_dat_annots, res_dat, ...)

	# non-ind sig table
	if (print_all_sig_table) {
		print_res_table(res_dat, ...)
		cat('-------------------------------------------------------\n')
	}

	cat('\n# ind sig genes:\n')

	#$ no PMR:
	ind_sig_which <- sig_dat[['dataset']] == SR_dataset_name & sig_dat[['final_ind_sig']] == 1

	res_dat_annots <- apply(sig_dat[ind_sig_which, ], 1, annot_letter)
	res_dat <- sig_dat[ind_sig_which, ]

	# info about annotations
	print_res_info(res_dat_annots, res_dat, ...)
	# ind_sig table
	cat('\n')
	print_res_table(res_dat, ...)
	cat('-------------------------------------------------------\n')

}

sig_dat$final_ind_sig <- as.integer(sig_dat$Gene %in% most_sig_no_overlap(SR_dataset_name, sig_dat[sig_dat$pmr_sig == 1, ]))

print(SR_dataset_name)
print_all_res(print_locus=TRUE, print_genes=TRUE, prop_merge_tab=TRUE)
print_all_res(print_locus=FALSE, print_genes=TRUE, prop_merge_tab=TRUE)

