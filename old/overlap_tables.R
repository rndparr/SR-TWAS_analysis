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

####
### TWAS to do
# twas <- 'AD_TWAS'
twas <- 'PD_TWAS'
# twas <- 'AD_TWAS2'

### directories
## HGCC
in_dir <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/', twas, '/')
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'
plot_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'

## laptop
# in_dir <- paste0('~/HGCC/YangFSSdata/SR_TWAS/', twas, '/')
# out_dir <- '~/HGCC/YangFSSdata/SR_TWAS/plots/'


#############################
## load data if not already made
load(paste0(out_dir, 'overlap_tables_data_', twas, '.RData'))


#############################
## MAKE DATA

format_pval <- function(x){
	return(format(x, scientific = TRUE, digits = 3))
}

wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
  options(width=as.integer(howWide))
}

wideScreen()

### read in data
dat <- read.table(
	paste0(in_dir, 'all_', twas, '_results.txt'), 
	header=TRUE, sep='\t')

dat$GeneName <- trimws(as.character(dat$GeneName))

if (twas == 'PD_TWAS'){
	## PD-TWAS
	# dat <- dat[dat$dataset %in% c('GTEx6tissues_GTEx_LD_GTEx','GTEx6tissues_Naive_GTEx_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP'),]
	dat <- dat[dat$dataset %in% c('GTEx6tissues_GTEx_LD_GTEx', 
		# 'GTEx6tissues_Naive_GTEx_LD_GTEx',
		'GTEx_Brain_Substantia_nigra_LD_GTEx', 
		'GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 
		'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 
		'GTEx_Brain_Cortex_LD_GTEx', 
		'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx',
		'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 
		'GTEx_Whole_Blood_LD_GTEx'),]

} else if (twas == 'AD_TWAS2'){
	## AD-TWAS2
	# ROSMAP was TIGAR DPR ROSMAP job
	# ROSMAP2 was prediXcan ROSMAP job
	dat <- dat[dat$dataset %in% c('GTEx','ROSMAP_LD_ROSMAP','PredDB','ROSMAP2_LD_ROSMAP','SR4basemodels_ROSMAP_LD_ROSMAP', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP'),]
}

model_order <- c('PredDB','GTEx','ROS','MAP','ROSMAP','ROSMAP2','Naive_MAP','Naive_ROS','Naive3_ROSMAP','SR_MAP','SR_ROS','SR_ROSMAP','Naive_ROS_MAP','SR_ROS_MAP','Naive34basemodels_ROSMAP','SR4basemodels_ROSMAP','GTEx_LD_ROSMAP','ROSMAP_LD_ROSMAP','PredDB_LD_ROSMAP','ROSMAP2_LD_ROSMAP','SR4basemodels_ROSMAP_LD_ROSMAP','Naive34basemodels_ROSMAP_LD_ROSMAP','GTEx6tissues_GTEx_LD_GTEx', 'GTEx6tissues_Naive_GTEx_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP','GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx')

names_col <- c('PrediXcan_GTEx', 'PrediXcan_ROSMAP', paste0('TIGAR_', c('GTEx','ROS','MAP','ROSMAP')), 'Naive_MAP','Naive_ROS','Naive_ROSMAP','SR_MAP','SR_ROS','SR_ROSMAP','Naive_ROS_MAP','SR_ROS_MAP', 'Naive_ROSMAP', 'SR_ROSMAP','TIGAR_GTEx_LD_ROSMAP','TIGAR_ROSMAP_LD_ROSMAP','PrediXcan_GTEx_LD_ROSMAP','PrediXcan_ROSMAP_LD_ROSMAP','SR_LD_ROSMAP','Naive_LD_ROSMAP','SR_GTEx6_LD_GTEx', 'Naive_GTEx6_LD_GTEx','GTEx_BrainSubNigra_LD_GTEx','ROSMAP_RNAseq_LD_ROSMAP','GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx')
raw_names_col <- c('PredDB','ROSMAP2','GTEx','ROS','MAP','ROSMAP','Naive_MAP','Naive_ROS','Naive3_ROSMAP','SR_MAP','SR_ROS','SR_ROSMAP','Naive_ROS_MAP','SR_ROS_MAP','Naive34basemodels_ROSMAP','SR4basemodels_ROSMAP','GTEx_LD_ROSMAP','ROSMAP_LD_ROSMAP','PredDB_LD_ROSMAP','ROSMAP2_LD_ROSMAP','SR4basemodels_ROSMAP_LD_ROSMAP','Naive34basemodels_ROSMAP_LD_ROSMAP','GTEx6tissues_GTEx_LD_GTEx', 'GTEx6tissues_Naive_GTEx_LD_GTEx','GTEx_Brain_Substantia_nigra_LD_GTEx', 'ROSMAP_RNAseq_brain_expr_WGS_b38_LD_ROSMAP','GTEx_Brain_Anterior_cingulate_cortex_BA24_LD_GTEx', 'GTEx_Brain_Caudate_basal_ganglia_LD_GTEx', 'GTEx_Brain_Cortex_LD_GTEx', 'GTEx_Brain_Nucleus_accumbens_basal_ganglia_LD_GTEx', 'GTEx_Brain_Putamen_basal_ganglia_LD_GTEx', 'GTEx_Whole_Blood_LD_GTEx')
cbind(names_col, raw_names_col)

factor_labels_raw <- setNames(names_col, raw_names_col)


# may be missing some models due to jobs not being finished
factor_levels <- model_order[which(model_order %in% unique(dat$dataset))]
factor_labels <- setNames(factor_labels_raw[factor_levels], NULL)
dat$dataset <- factor(dat$dataset, levels=factor_levels, labels=factor_labels)

# GET ZSCORE COLUMN
method <- 'SPred'
# method <- 'FUSION'
colnames(dat)[which(colnames(dat) == paste0(method,'_Z'))] <- 'Zscore'
# endsWith(colnames(dat), '_Z')

## GET PVALUE
dat['Pvalue'] <- exp(pchisq(dat[['Zscore']]^2, 1, lower.tail=FALSE, log.p=TRUE))
dat['Pvalue_raw'] <- dat$Pvalue
dat <- dat[order(dat$dataset, dat$CHROM, dat$Pvalue), ]

## RENAME COLUMNS
# colnames(dat)[c(2:4,6)] <- c('Chrom', 'Start', 'End', 'Gene')
colnames(dat)[which(colnames(dat) == 'CHROM')] <- 'Chrom'
colnames(dat)[which(colnames(dat) == 'GeneStart')] <- 'Start'
colnames(dat)[which(colnames(dat) == 'GeneEnd')] <- 'End'
colnames(dat)[which(colnames(dat) == 'GeneName')] <- 'Gene'

dat$Gene <- trimws(dat$Gene)
dat$TargetIDold <- dat$TargetID
dat$TargetID <- do.call(rbind, strsplit(dat$TargetID, '.', fixed=TRUE))[,1]

dat_raw <- dat


###############
## GWAS CATALOG DATA
# assoc_locs_genes <- read.csv(
# 	paste0(plot_dir, twas, '_GWAS_catalog.csv'))

assoc_locs_genes <- read.csv(
	paste0(plot_dir, list('AD_TWAS'='AD_TWAS', 'PD_TWAS'='PD_TWAS', 'AD_TWAS2'='AD_TWAS')[[twas]], '_GWAS_catalog.csv'))


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
strsplits <- function(x, splits, ...){
	for (split in splits) {
		x <- unlist(strsplit(x, split, ...))
	}
	return(x[!x == ""]) # Remove empty values
}
assoc_genes <- strsplits(unique(assoc_locs_genes$Mapped.gene), c('[[:space:]]', ',', ';'))


# list of genes with locations
dat$prev_id_gwas <- as.integer(dat$Gene %in% assoc_genes)
prev_id_gwas_dat <- dat[dat$prev_id_gwas > 0, c('Chrom', 'Start', 'End', 'prev_id_gwas')]
prev_id_gwas_dat <- prev_id_gwas_dat[order(prev_id_gwas_dat$Chrom, prev_id_gwas_dat$Start), ]
prev_id_gwas_dat <- prev_id_gwas_dat[!duplicated(prev_id_gwas_dat), ]

##########
# function for GWAS_CATALOG comparison

# near snp
near_gwas_snp <- function(row){
	chrom <- as.integer(row[['Chrom']])
	wstart <- as.integer(row[['w_start']])
	wend <- as.integer(row[['w_end']])

	chrom_assoc_locs <- assoc_locs[as.integer(assoc_locs[['Chrom']]) == chrom, 'loc']

	return(as.integer(any(chrom_assoc_locs >= wstart & chrom_assoc_locs <= wend)))
}

# near gene; (search whole dat dataframe, but only for rows in sig_data)
# near_gwas_gene <- function(row, data = prev_id_gwas_dat){
# 	chrom <- as.integer(row[['Chrom']])
# 	wstart <- as.integer(row[['w_start']])
# 	wend <- as.integer(row[['w_end']])

# 	locs <- data[data[['Chrom']] == chrom & data[['prev_id_gwas']] == 1, ]
# 	locs <- locs[!duplicated(locs),]
# 	locs <- locs[,c('Start', 'End')]

# 	return(as.integer(any(locs >= wstart & locs <= wend)))
# }

near_gwas_gene <- function(row, data = prev_id_gwas_dat){
	chrom <- as.integer(row[['Chrom']])
	wstart <- as.integer(row[['w_start']])
	wend <- as.integer(row[['w_end']])

	locs <- data[data[['Chrom']] == chrom & data[['prev_id_gwas']] == 1, ]
	if (nrow(locs) == 0) {
		return(0)
	}
	locs <- locs[!duplicated(locs),]
	locs <- locs[,c('Start', 'End')]

	return(as.integer(any(locs >= wstart & locs <= wend)))
}

# near either
near_gwas <- function(row, data = prev_id_gwas_dat){
	out <- (near_gwas_gene(row, data) + near_gwas_snp(row)) > 0
	return(as.integer(out))
}

##########
## FUNCTIONS BASED ON GENE LOCATION AND SIGNIFICANCE

## SET UP SIG_DAT DATAFRAME FOR MAKING TABLES
sig_dat <- na.omit(dat[dat$Pvalue <= 2.5e-6, ])
rownames(sig_dat) <- NULL


# function to find most significant genes in 1MB window
most_sig_no_overlap <- function(dataset, data = sig_dat){
	table_data <- data[data$dataset==dataset & data$Pvalue_raw <= 2.5e-6, c('Gene','Chrom','Start','End', 'Zscore', 'Pvalue')]
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


# select = 'all' to get all matches
prev_id_gwas_overlap <- function(dataset, data = sig_dat, select = 'first'){

	filter_dat <- data[data$dataset==dataset & data$Pvalue_raw <= 2.5e-6, c('Gene','Chrom','Start','End', 'Zscore','Pvalue', 'prev_id_gwas')]
	filter_dat <- na.omit(filter_dat)

	# previously id'd genes
	gr_gwas_dat <- filter_dat[filter_dat[['prev_id_gwas']] == 1, c('Gene','Chrom','Start','End')]
	colnames(gr_gwas_dat) <- c('name','chr','start','end')
	gr_gwas_dat <- gr_gwas_dat[!duplicated(gr_gwas_dat),]

	# other genes
	gr_dat <- filter_dat[filter_dat[['prev_id_gwas']] != 1, c('Gene','Chrom','Start','End')]
	colnames(gr_dat) <- c('name','chr','start','end')

	# intersect overlapping columns
	chr_overlap <- intersect(unique(gr_gwas_dat$chr), unique(gr_dat$chr))
	gr_gwas_dat <- gr_gwas_dat[gr_gwas_dat$chr %in% chr_overlap, ]

	row.names(gr_gwas_dat) <- NULL
	row.names(gr_dat) <- NULL

	# make ranges
	gr_gwas <- makeGRangesFromDataFrame(gr_gwas_dat, keep.extra.columns=TRUE)
	gr <- makeGRangesFromDataFrame(gr_dat, keep.extra.columns=TRUE)

	### conservative only
	gr_gwas_wide <- resize(gr_gwas, width = 2000000 + width(gr_gwas), fix = 'center') 
	gr_wide <- resize(gr, width = 2000000 + width(gr), fix = 'center')
	overlaps <- findOverlaps(gr_wide, gr_gwas_wide, type = 'any', select = select) 

	### non conservative
	# overlaps <- findOverlaps(gr, gr_gwas, maxgap=1000000, type = 'any', select = 'first') # non-conservative
	if (select == 'first'){
		gr_dat$previd_overlap <- gr_gwas_dat$name[overlaps]
		gr_dat_overlapped_previd <- gr_dat
		row.names(gr_dat_overlapped_previd) <- NULL
	} else if (select == 'all'){
		# get ALL matches as comma separated string
		odf <- as.data.frame(overlaps)
		odf$subjectHits <- gr_gwas_dat$name[as.integer(odf$subjectHits)]
		odf$queryHits <- gr_dat$name[as.integer(odf$queryHits)]
		prev_id_gwas_matches <- as.data.frame(cbind(by(odf$subjectHits, odf$queryHits, paste, collapse=',')))
		prev_id_gwas_matches$name <- row.names(prev_id_gwas_matches)
		colnames(prev_id_gwas_matches) <- c('prev_id_gwas_overlap', 'name')
		row.names(prev_id_gwas_matches) <- NULL

		gr_dat_overlapped_previd <- merge(gr_dat, prev_id_gwas_matches, all = TRUE)
		row.names(gr_dat_overlapped_previd) <- NULL
	}

	# return genes that are near previd'd, with column denoting which gene they're near
	return(gr_dat_overlapped_previd)
}


#### SEEING IF THERE'S A PVAL, EVEN NON SIG, FOR ALL METHODS FOR GENES SIGNIFICANT IN SR_ROSMAP
# windows for each gene
dat$w_start <- as.integer(dat$Start - 1000000)
dat$w_end <- as.integer(dat$End + 1000000)

# significant
for(dataset in levels(dat$dataset)){
	sig_str <- paste0('sig_', dataset)
	dat[[sig_str]] <- as.integer(dat$Gene %in% dat[dat$dataset == dataset & dat$Pvalue <= 2.5e-6, 'Gene'])
}

# most sig in region
for(dataset in levels(dat$dataset)){
	most_sig_str <- paste0('most_sig_', dataset)
	dat[[most_sig_str]] <- as.integer(dat$Gene %in% most_sig_no_overlap(dataset))
}

dat[['near_prev_id_gwas_snp']] <- apply(dat, 1, near_gwas_snp)
dat[['near_prev_id_gwas_gene']] <- apply(dat, 1, near_gwas_gene)
# dat[['near_prev_id_gwas']] <- apply(dat, 1, near_gwas)
dat[['near_prev_id_gwas']] <- ifelse(dat$near_prev_id_gwas_snp > 0 | dat$near_prev_id_gwas_gene > 0, 1, 0)


# near_gwas_gene_trycatch <- function(row, data = prev_id_gwas_dat) {
# 	out <- tryCatch({ near_gwas_gene(row, data) },
# 		error=function(cond){ message(row); return(NA) },
# 		warning=function(cond){return(NA)}
# 		)
# 	return(out)
# }


# dat$most_sig_n <- rowSums(dat[,colnames(dat)[which(startsWith(colnames(dat), 'most_sig'))]])
# dat$most_sig_any <- as.integer(rowSums(dat[,colnames(dat)[which(startsWith(colnames(dat), 'most_sig'))]]) > 0)


################
## UPDATE SIG_DAT DATAFRAME FOR MAKING TABLES
sig_dat <- na.omit(dat[dat$Pvalue <= 2.5e-6, ])
rownames(sig_dat) <- NULL

sig_dat$most_sig_n <- rowSums(sig_dat[,colnames(sig_dat)[which(startsWith(colnames(sig_dat), 'most_sig'))]])
sig_dat$most_sig_any <- as.integer(rowSums(sig_dat[,colnames(sig_dat)[which(startsWith(colnames(sig_dat), 'most_sig'))]]) > 0)

# head(sig_dat)

# gene names if most sig in region for any dataset
ind_genes <- sort(unique(sig_dat[sig_dat$most_sig_any > 0, 'Gene']))

paste(sort(unique(sig_dat[sig_dat$most_sig_any > 0, 'Gene'])), sep="','", collapse="','")
## 'APOE','CEACAM19','CLASRP','CR1','DDAH2','DMPK','GIPR','HARBI1','HLA-DRA','HLA-DRB1','MS4A7','PICALM','PROB1','SLC15A3','SLC39A13','SPINDOC','TMEM173','TRAPPC6A','VWA7','ZBTB12'

################################################################################################

load(paste0(plot_dir, 'prev_id_twas_', twas, '.RData'))
prev_id_dat$Gene <- trimws(as.character(prev_id_dat$Gene))
# loads dataset called prev_id_dat

# near gene; (search whole dat dataframe, but only for rows in sig_data)
near_twas <- function(row, data = prev_id_dat){
	chrom <- as.integer(row[['Chrom']])
	wstart <- as.integer(row[['w_start']])
	wend <- as.integer(row[['w_end']])

	locs <- data[data[['Chrom']] == chrom, ]
	if (length(locs) == 0) {
		return(0)
	}
	locs <- locs[!duplicated(locs),]
	locs <- locs[,c('Start', 'End')]

	return(as.integer(any(locs >= wstart & locs <= wend)))
}

annot_name <- function(row){
	search <- c('prev_id_gwas', 'near_prev_id_gwas', 'prev_id_twas', 'near_prev_id_twas')

	# gwas
	if (row[[search[1]]] == 1) {
		# prev id
		note <- letters[1]
	} else if (row[[search[2]]] == 1) {
		# near
		note <- letters[2]
	} else {
		note <- ''
	}

	# twas
	if (row[[search[3]]] == 1) {
		# prev id
		note <- paste0(note, letters[3])
	} else if (row[[search[4]]] == 1) {
		# near
		note <- paste0(note, letters[4])
	}

	# add math notation
	if (note != '') {
		note <- paste0('$^{', note, '}$')
	}

	return(paste0(row[['Gene']], note))
}

annot_letter <- function(row){
	search <- c('prev_id_gwas', 'near_prev_id_gwas', 'prev_id_twas', 'near_prev_id_twas')
	# gwas
	if (row[[search[1]]] == 1) {
		# prev id
		note <- letters[1]
	} else if (row[[search[2]]] == 1) {
		# near
		note <- letters[2]
	} else {
		note <- ''
	}
	# twas
	if (row[[search[3]]] == 1) {
		# prev id
		note <- paste0(note, letters[3])
	} else if (row[[search[4]]] == 1) {
		# near
		note <- paste0(note, letters[4])
	}
	if (note == '') {
		note <- 0
	}
	return(note)
}

# levels(sig_dat$dataset)

sig_dat$prev_id_twas <- as.integer(sig_dat$TargetID %in% prev_id_dat$TargetID)
sig_dat$near_prev_id_twas <- apply(sig_dat, 1, near_twas)
sig_dat$GeneName <- apply(sig_dat, 1, annot_name)


###############

# save(file=paste0(out_dir, 'overlap_tables_data_', twas, '.RData'))

# save.image(file=paste0(out_dir, 'overlap_tables_data_', twas, '.RData'))

## END MAKE DATA
####################################################################################################################

###############



if (twas == 'PD_TWAS'){
	## PD-TWAS
	SR_dataset_name <- 'SR_GTEx6_LD_GTEx'
	other_dataset_name <- 'GTEx_BrainSubNigra_LD_GTEx'

} else if (twas == 'AD_TWAS2'){
	## AD-TWAS2
	SR_dataset_name <- 'SR_LD_ROSMAP'
	other_dataset_name <- 'ROSMAP_RNAseq_LD_ROSMAP'
}


print_res_info <- function(res_dat_annots, res_dat){
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

	## prep formatted
	res_vec <- table(res_dat_annots)
	not_in <- setdiff(names(name_convert), names(res_vec))
	res_vec_letters <- c(setNames(rep(0, length(not_in)), not_in), res_vec)[names(name_convert)]
	res_vec <- res_vec_letters
	# res_vec
	names(res_vec) <- name_convert[names(res_vec)]

	cat('-------------------------------------------------------\n')
	cat(paste0('total: ', nrow(res_dat), '\n'))
	cat('-------------------------\n')
	print(table(unlist(strsplit(res_dat_annots, '')), dnn=NULL))
	cat('-------------------------\n')
	print(res_vec_letters, dnn=NULL)
	cat('-------------------------\n')
	cat(paste(capture.output(print(cbind(res_vec)))[-1], collapse='\n'), '\n')
	cat('-------------------------\n')
	cat(sapply(names_order, format_res_row, USE.NAMES=FALSE), sep='')
	cat('-------------------------------------------------------\n')
}


print_res_table <- function(table_data_sr){
	### merge other data
	table_data_other <- sig_dat[sig_dat$dataset == other_dataset_name , c('GeneName', 'Zscore', 'Pvalue')]
	colnames(table_data_other) <- c('GeneName', 'Zscore_V', 'Pvalue_V')

	table_data_raw <- merge(table_data_sr, table_data_other, all.x=TRUE, by='GeneName')
	table_data_raw <- table_data_raw[order(table_data_raw$Chrom, table_data_raw$Start, table_data_raw$Pvalue), ]

	############
	## FORMAT TABLE
	table_data <- table_data_raw[, c('GeneName', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')]
	colnames(table_data) <- c('Gene', 'CHR', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')

	display_vec <- c('s','s','d','d','d','f','e','f','e')

	xtable_data <- xtable(table_data,
		math.style.exponent=TRUE, booktabs=TRUE, display=display_vec)

	print(xtable_data, booktabs=TRUE, include.rownames=FALSE, 
		sanitize.text.function=identity, 
		table.placement='H', comment=FALSE,
		NA.string='-')
}


print_all_res <- function(print_all_sig_table=FALSE){
	# all sig
	# non-ind sig table
	cat(paste('* ', twas, '\n# all sig genes:\n'))
	res_dat_annots <- apply(sig_dat[sig_dat[['dataset']] == SR_dataset_name, ], 1, annot_letter)
	res_dat <- sig_dat[sig_dat[['dataset']] == SR_dataset_name, ]

	# info about annotations
	print_res_info(res_dat_annots, res_dat)

	# non-ind sig table
	if (print_all_sig_table) {
		print_res_table(res_dat)
		cat('-------------------------------------------------------\n')
	}

	# cat(paste('\n* ', twas, 'ind sig genes:\n'))

	cat('\n# ind sig genes:\n')

	# ind sig
	res_dat_annots <- apply(sig_dat[sig_dat[['dataset']] == SR_dataset_name & sig_dat[[paste0('most_sig_', SR_dataset_name)]] == 1, ], 1, annot_letter)
	res_dat <- sig_dat[sig_dat[['dataset']] == SR_dataset_name & sig_dat[[paste0('most_sig_', SR_dataset_name)]] == 1, ]

	# info about annotations
	print_res_info(res_dat_annots, res_dat)
	# ind_sig table
	print_res_table(res_dat)
	cat('-------------------------------------------------------\n')

}

print_all_res()

################################

bleh <- summary(dat$dataset)
bleh2 <- cbind(names(bleh), bleh)
colnames(bleh2) <- c('', 'N_results')
bleh2 <- as.data.frame(bleh2)
bleh2$N_results <- as.integer(bleh2$N_results)
row.names(bleh2) <- NULL
bleh2[,c('V1','N_results')]

print(bleh2, row.names=FALSE)

################################

## see if figure out way to merge all base model results too
table_data_sr <- res_dat


display_vec <- c('s','s','d','d','d','f','e')
table_data_cols <- c('Gene', 'CHR', 'Start', 'End', 'Zscore', 'Pvalue')


	table_data <- table_data_raw[, c('GeneName', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')]
	colnames(table_data) <- c('Gene', 'CHR', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')

for (other_dataset_name in ) {

	### merge other data
	table_data_other <- sig_dat[sig_dat$dataset == other_dataset_name , c('GeneName', 'Zscore', 'Pvalue')]

	colnames(table_data_other) <- c('GeneName', 'Zscore_V', 'Pvalue_V')

	table_data_raw <- merge(table_data_sr, table_data_other, all.x=TRUE, by='GeneName')
	table_data_raw <- table_data_raw[order(table_data_raw$Chrom, table_data_raw$Start, table_data_raw$Pvalue), ]

	############
	## FORMAT TABLE
	table_data <- table_data_raw[, c('GeneName', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')]
	colnames(table_data) <- c('Gene', 'CHR', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')

	display_vec <- c(display_vec,'f','e')

}


tail(sig_dat)

# # all sig
# res_dat_annots <- apply(sig_dat[sig_dat[['dataset']] == SR_dataset_name,], 1, annot_letter)
# res_dat <- sig_dat[sig_dat[['dataset']] == SR_dataset_name,]

# print_res(res_dat_annots, res_dat)

# # ind sig
# res_dat_annots <- apply(sig_dat[sig_dat[['dataset']] == SR_dataset_name & sig_dat[[paste0('most_sig_', SR_dataset_name)]] == 1,], 1, annot_letter)
# res_dat <- sig_dat[sig_dat[['dataset']] == SR_dataset_name & sig_dat[[paste0('most_sig_', SR_dataset_name)]] == 1,]


# bleh1 <- sig_dat[sig_dat[['dataset']] == SR_dataset_name & sig_dat[[paste0('sig_', SR_dataset_name)]] == 1, ]
# bleh2 <- sig_dat[sig_dat[['dataset']] == SR_dataset_name, ]
# all(bleh1==bleh2)

# ############
# ## ASSEMBLE LATEX TABLE WITHOUT COMPARISON
# # non-ind sig table
# table_data_sr <- sig_dat[sig_dat[['dataset']] == SR_dataset_name & sig_dat[[paste0('sig_', SR_dataset_name)]] == 1, ]

# # ind_sig table
# table_data_sr <- sig_dat[sig_dat[['dataset']] == SR_dataset_name & sig_dat[[paste0('most_sig_', SR_dataset_name)]] == 1, ]


# ############
# # ## FORMAT TABLE
# # table_data <- table_data_sr[, c('GeneName', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue')]
# # colnames(table_data) <- c('Gene', 'CHR', 'Start', 'End', 'Zscore', 'Pvalue')

# # display_vec <- c('s','s','d','d','d','f','e')

# # xtable_data <- xtable(table_data,
# # 	math.style.exponent=TRUE, booktabs=TRUE, display=display_vec)

# # print(xtable_data, booktabs=TRUE, include.rownames=FALSE, 
# # 	sanitize.text.function=identity, 
# # 	table.placement='H', comment=FALSE,
# # 	NA.string='-')

# ############
# ### merge other data

# # table_data_other <- sig_dat[sig_dat$dataset == other_dataset_name & sig_dat[[paste0('most_sig_', SR_dataset_name)]] == 1, c('GeneName', 'Zscore', 'Pvalue')]
# table_data_other <- sig_dat[sig_dat$dataset == other_dataset_name , c('GeneName', 'Zscore', 'Pvalue')]
# colnames(table_data_other) <- c('GeneName', 'Zscore_V', 'Pvalue_V')

# table_data_raw <- merge(table_data_sr, table_data_other, all.x=TRUE, by='GeneName')
# table_data_raw <- table_data_raw[order(table_data_raw$Chrom, table_data_raw$Start, table_data_raw$Pvalue), ]

# ############
# ## FORMAT TABLE
# table_data <- table_data_raw[, c('GeneName', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')]
# colnames(table_data) <- c('Gene', 'CHR', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_V', 'Pvalue_V')

# display_vec <- c('s','s','d','d','d','f','e','f','e')

# xtable_data <- xtable(table_data,
# 	math.style.exponent=TRUE, booktabs=TRUE, display=display_vec)

# print(xtable_data, booktabs=TRUE, include.rownames=FALSE, 
# 	sanitize.text.function=identity, 
# 	table.placement='H', comment=FALSE,
# 	NA.string='-')


# ############
# ### merge naive data
# table_data_naive <- sig_dat[sig_dat$dataset == 'Naive_ROSMAP' & sig_dat$most_sig_SR_ROSMAP == 1, c('GeneName', 'Zscore', 'Pvalue')]
# colnames(table_data_naive) <- c('GeneName', 'Zscore_N', 'Pvalue_N')

# table_data_raw <- merge(table_data_sr, table_data_naive, all.x=TRUE, by='GeneName')
# table_data_raw <- table_data_raw[order(table_data_raw$Chrom, table_data_raw$Start, table_data_raw$Pvalue), ]

# ############
# ## FORMAT TABLE
# table_data <- table_data_raw[, c('GeneName', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_N', 'Pvalue_N')]
# colnames(table_data) <- c('Gene', 'CHR', 'Start', 'End', 'Zscore', 'Pvalue', 'Zscore_N', 'Pvalue_N')

# display_vec <- c('s','s','d','d','d','f','e','f','e')

# xtable_data <- xtable(table_data,
# 	math.style.exponent=TRUE, booktabs=TRUE, display=display_vec)

# print(xtable_data, booktabs=TRUE, include.rownames=FALSE, 
# 	sanitize.text.function=identity, 
# 	table.placement='H', comment=FALSE,
# 	NA.string='-')




################################################################################################

################
# N significant genes by method
n_genes_all <- data.frame('N_sig_genes'=t(t(as.matrix(table(sig_dat$dataset)))))
n_genes_all$dataset <- rownames(n_genes_all)
n_genes_all <- n_genes_all[order(n_genes_all$N_sig_genes), c('dataset','N_sig_genes')]
rownames(n_genes_all) <- NULL

# print(n_genes_all, row.names=FALSE, justify='left')

################
### NUMBER OF INDEPENDENT
n_ind_genes <- data.frame()

for (dataset in levels(dat$dataset)) {
	most_sig_str <- paste0('most_sig_', dataset)
	n_ind_genes <- rbind(n_ind_genes, data.frame(dataset=dataset, N_ind_genes=sum(dat[dat$dataset == dataset, most_sig_str])))
}
n_ind_genes <- n_ind_genes[order(n_ind_genes$N), c('dataset','N_ind_genes')]


# print(n_ind_genes, row.names=FALSE)

################
## COMBINE SIG GENES, NUMBER INDEPENDENT
ngenes <- merge(n_ind_genes, n_genes_all)
ngenes <- ngenes[order(ngenes$N_ind_genes, ngenes$N_sig_genes),]

# levels_models <- c('PrediXcan_GTEx','TIGAR_GTEx','TIGAR_ROSMAP', 'Naive_ROSMAP', 'SR_ROSMAP')
levels_models <- c('PrediXcan_GTEx','PrediXcan_ROSMAP','TIGAR_GTEx','TIGAR_ROSMAP', 'Naive_ROSMAP', 'SR_ROSMAP')
ngenes$dataset <- factor(ngenes$dataset, levels=levels_models)
ngenes <- ngenes[order(ngenes$dataset),]

print(ngenes, row.names=FALSE)

# ngenes[c('PrediXcan_GTEx','TIGAR_GTEx','TIGAR_ROSMAP', 'Naive_ROSMAP', 'SR_ROSMAP'),]

# ngenes[ngenes$dataset %in% c('PrediXcan_GTEx','TIGAR_GTEx','TIGAR_ROSMAP', paste0(c('Naive','SR'), '_', 'ROSMAP')),]



#######################################################################################################
suffix <- '_4basemodels'
# suffix <- ''
save(sig_dat, dat, ind_genes, ngenes, file=paste0(in_dir, 'overlap_table_dat_', twas, suffix, '.Rdata'))




#######################################################################################################
################
# ROSMAP DATA
# sr_rosmap_sig_genes <- sig_dat[sig_dat$dataset == 'SR_ROSMAP', 'Gene']
sr_rosmap_sig_genes <- sig_dat[sig_dat$dataset == 'SR_ROSMAP' & sig_dat$most_sig_SR_ROSMAP==1, 'Gene']

# get gene position info from SR_ROSMAP rows
sr_rosmap_sig_gene_info <- dat[dat$dataset == 'SR_ROSMAP' & dat$Gene %in% sr_rosmap_sig_genes, c('Chrom','Start','End','Gene','most_sig_SR_ROSMAP')]
rownames(sr_rosmap_sig_gene_info) <- NULL
sr_rosmap_sig_gene_info

## get dataframe
sr_rosmap_sig_gene_dat <- merge(
	sr_rosmap_sig_gene_info,
	dcast(Chrom + Gene ~ dataset, value.var='Pvalue', data=dat[dat$Gene %in% sr_rosmap_sig_genes, ]))
sr_rosmap_sig_gene_dat <- with(sr_rosmap_sig_gene_dat, sr_rosmap_sig_gene_dat[order(Chrom, SR_ROSMAP, Start), -5])


## formatting non-info cols
for (col in colnames(sr_rosmap_sig_gene_dat)[-c(1:4)]){
	sr_rosmap_sig_gene_dat[[col]] <- trimws(format_pval(sr_rosmap_sig_gene_dat[[col]]))
}
sr_rosmap_sig_gene_dat[sr_rosmap_sig_gene_dat=='NA'] <- '-'
rownames(sr_rosmap_sig_gene_dat) <- NULL

sr_rosmap_sig_gene_dat

## for printing
final_models <- c('PrediXcan_GTEx','TIGAR_GTEx','TIGAR_ROSMAP', paste0(c('Naive','SR'), '_', 'ROSMAP'))
sr_rosmap_sig_gene_dat_clean <- sr_rosmap_sig_gene_dat[, c('Chrom','Start','End','Gene', final_models)]
colnames(sr_rosmap_sig_gene_dat_clean)[8:9] <- c('Naive','SR-TWAS')

print(sr_rosmap_sig_gene_dat_clean, row.names=FALSE)

print(sr_rosmap_sig_gene_dat_clean, digits=3, row.names=FALSE)

print(sr_rosmap_sig_gene_dat_clean[,c(1,4:9)], digits=3, row.names=FALSE)

##
sr_rosmap_sig_gene_dat_clean$prev_id <- as.integer(sr_rosmap_sig_gene_dat_clean$Gene %in% assoc_genes)
colnames(dat)

sr_rosmap_sig_gene_dat_clean$w_start <- as.integer(sr_rosmap_sig_gene_dat_clean$Start - 1000000)
sr_rosmap_sig_gene_dat_clean$w_end <- as.integer(sr_rosmap_sig_gene_dat_clean$End + 1000000)

sr_rosmap_sig_gene_dat_clean[['near_prev_id_snp']] <- apply(sr_rosmap_sig_gene_dat_clean, 1, near_gwas_snp)
sr_rosmap_sig_gene_dat_clean[['near_prev_id_gwas']] <- apply(sr_rosmap_sig_gene_dat_clean, 1, near_gwas_gene)
sr_rosmap_sig_gene_dat_clean[['near_prev_id']] <- apply(sr_rosmap_sig_gene_dat_clean, 1, near_gwas)

sr_rosmap_sig_gene_dat_clean


### PD-TWAS, UNIQUELY ID'D BY NAIVE, SR-TWAS
sr_rosmap_sig_gene_dat_clean[sr_rosmap_sig_gene_dat_clean$Gene %in% c('NUCKS1','KLHL7-DT','CNTN1','RNF40'),]
## RNF40 only one not previously ID'd

################
# ROS_MAP DATA
sr_ros_map_sig_genes <- sig_dat[sig_dat$dataset == 'SR_ROS_MAP', 'Gene']

# get gene position info from SR_ROS_MAP rows
sr_ros_map_sig_gene_info <- dat[dat$dataset == 'SR_ROS_MAP' & dat$Gene %in% sr_ros_map_sig_genes, c('Chrom','Start','End','Gene','most_sig_SR_ROS_MAP')]
rownames(sr_ros_map_sig_gene_info) <- NULL

## get dataframe
sr_ros_map_sig_gene_dat <- merge(
	sr_ros_map_sig_gene_info,
	dcast(Chrom + Gene ~ dataset, value.var='Pvalue', data=dat[dat$Gene %in% sr_ros_map_sig_genes, ]))
sr_ros_map_sig_gene_dat <- with(sr_ros_map_sig_gene_dat, sr_ros_map_sig_gene_dat[order(Chrom, SR_ROS_MAP, Start), -5])

## formatting non-info cols
for (col in colnames(sr_ros_map_sig_gene_dat)[-c(1:4)]){
	sr_ros_map_sig_gene_dat[[col]] <- trimws(format_pval(sr_ros_map_sig_gene_dat[[col]]))
}
sr_ros_map_sig_gene_dat[sr_ros_map_sig_gene_dat=='NA'] <- '-'
rownames(sr_ros_map_sig_gene_dat) <- NULL

sr_ros_map_sig_gene_dat

## for printing
final_models <- c('PrediXcan_GTEx','TIGAR_GTEx',paste0('TIGAR_', unlist(strsplit('ROS_MAP','_'))), paste0(c('Naive','SR'), '_', 'ROS_MAP'))
sr_ros_map_sig_gene_dat_clean <- sr_ros_map_sig_gene_dat[, c('Chrom','Start','End','Gene', final_models)]
colnames(sr_ros_map_sig_gene_dat_clean)[9:10] <- c('Naive','SR-TWAS')

print(sr_ros_map_sig_gene_dat_clean, row.names=FALSE)



# ####################################################################################################

################
# ALL DATA
sig_genes <- unique(sig_dat$Gene)

# get gene position info from SR_ROSMAP rows
sig_gene_info <- dat_raw[dat_raw$Gene %in% sig_genes, c('Chrom','Start','End','Gene')]
sig_gene_info <- sig_gene_info[nrow(sig_gene_info):1,]
sig_gene_info <- sig_gene_info[!duplicated(sig_gene_info$Gene),]
rownames(sig_gene_info) <- NULL
sig_gene_info


## get dataframe
sig_gene_dat <- merge(
	sig_gene_info,
	dcast(Chrom + Gene ~ dataset, value.var='Pvalue', data=dat_raw[dat_raw$Gene %in% sig_genes, ]))
sig_gene_dat <- with(sig_gene_dat, sig_gene_dat[order(Chrom, SR_ROSMAP, Start), ])


## formatting non-info cols
for (col in colnames(sr_rosmap_sig_gene_dat)[-c(1:4)]){
	sig_gene_dat[[col]] <- trimws(format_pval(sig_gene_dat[[col]]))
}
sig_gene_dat[sig_gene_dat=='NA'] <- '-'
rownames(sig_gene_dat) <- NULL

print(sig_gene_dat, row.names=FALSE)



########
sig_gene_info <- dat_raw[dat_raw$Gene %in% sig_genes, c('Chrom','Start','End','Gene')]

rownames(sig_gene_info) <- NULL
sig_gene_info

setdiff(sig_genes, sig_gene_info$Gene)



sig_gene_dat$w_start <- as.integer(sig_gene_dat$Start - 1000000)
sig_gene_dat$w_end <- as.integer(sig_gene_dat$End + 1000000)

sig_gene_dat[['near_prev_id_snp']] <- apply(sig_gene_dat, 1, near_gwas_snp)
sig_gene_dat[['near_prev_id_gwas']] <- apply(sig_gene_dat, 1, near_gwas_gene)
sig_gene_dat[['near_prev_id']] <- apply(sig_gene_dat, 1, near_gwas)

sig_gene_dat$prev_id <- as.integer(sig_gene_dat$Gene %in% assoc_genes)

sig_gene_dat$ind_sig <- as.integer(sig_gene_dat$Gene %in% assoc_genes)

sort(unique(unlist(vennlist2)))



sig_gene_dat$no_overlap <- as.integer(sig_gene_dat$Gene %in% no_overlap(sig_gene_dat))
sig_gene_dat[sig_gene_dat$Gene %in% ]

no_overlap(sig_gene_dat)

no_overlap <- function(data = sig_dat){
	table_data <- data[,c('Gene','Chrom','Start','End')]
	table_data <- na.omit(table_data)

	gr_dat <- table_data[,c('Gene','Chrom','Start','End')]
	gr_dat <- table_data[!duplicated(table_data),]
	colnames(gr_dat) <- c('name','chr','start','end')
	gr <- makeGRangesFromDataFrame(gr_dat, keep.extra.columns=TRUE)
	gr_wide <- resize(gr, width = 2000000 + width(gr), fix = 'center')
	overlaps <- findOverlaps(gr_wide, type = 'any', select = 'first')
	# overlaps <- findOverlaps(gr, maxgap=1000000, type = 'any', select = 'first')
	most_sig_names <- as.data.frame(gr[!duplicated(overlaps)])[['name']]
	return(most_sig_names)
}




################################################################################################
bleh <- function(x){
	return(sig_dat[sig_dat$dataset == x , 'Gene'])
}

bleh2 <- function(x){
	return(sig_dat[sig_dat$dataset == x & sig_dat[[paste0('most_sig_', x)]] == 1, 'Gene'])
}

# set directory, tissue
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'

vennlist <- sapply(levels(sig_dat$dataset), bleh)[c('PrediXcan_GTEx','TIGAR_GTEx','TIGAR_ROSMAP', paste0(c('Naive','SR'), '_', 'ROSMAP'))]
dput(vennlist)

sort(unique(unlist(vennlist)))

intersect(vennlist[['TIGAR_GTEx']], vennlist[['SR_ROSMAP']])
setdiff(vennlist[['TIGAR_GTEx']], vennlist[['SR_ROSMAP']])
setdiff( vennlist[['SR_ROSMAP']],vennlist[['TIGAR_GTEx']])

setdiff( vennlist[['SR_ROSMAP']],vennlist[['TIGAR_ROSMAP']])
setdiff( vennlist[['SR_ROSMAP']],vennlist[['PrediXcan_GTEx']])


vennlist2 <- sapply(levels(sig_dat$dataset), bleh2)[c('PrediXcan_GTEx','TIGAR_GTEx','TIGAR_ROSMAP', paste0(c('Naive','SR'), '_', 'ROSMAP'))]
dput(vennlist2)

dput(c(sort(unique(unlist(vennlist2))), c('NUCKS1','KLHL7-DT','CNTN1','RNF40')))
dput(sort(intersect(sr_rosmap_sig_gene_dat_clean$Gene, unique(unlist(vennlist)))))
dput(sort(c(intersect(sr_rosmap_sig_gene_dat_clean$Gene, unique(unlist(vennlist2))), c('NUCKS1','KLHL7-DT','CNTN1','RNF40'))))

c(intersect(sr_rosmap_sig_gene_dat_clean$Gene, unique(unlist(vennlist2))), c('NUCKS1','KLHL7-DT','CNTN1','RNF40'))

dput(sort(setdiff(c(intersect(sr_rosmap_sig_gene_dat_clean$Gene, unique(unlist(vennlist2))), c('NUCKS1','KLHL7-DT','CNTN1','RNF40'))
, c("CD38" ,"GPNMB","STX4", "VKORC1"))))

## genes I'd by all
Reduce(intersect, vennlist)

## ind genes I'd by all
Reduce(intersect, vennlist2)

print(sr_rosmap_sig_gene_dat_clean[sr_rosmap_sig_gene_dat_clean$Gene %in% sort(unique(unlist(vennlist2))), ], row.names=FALSE)


print(sr_rosmap_sig_gene_dat_clean, row.names=FALSE)


sr_rosmap_sig_gene_dat_clean[sr_rosmap_sig_gene_dat_clean$Gene %in% sort(unique(unlist(vennlist2))), ]

########
library(VennDiagram)
venn.diagram(
	x=vennlist,
	filename=paste0(out_dir, twas, '_venn_diagram.png')
	)

venn.diagram(
	x=vennlist2,
	filename=paste0(out_dir, twas, '_venn_diagram2.png')
	)



# # ####################################################################################################