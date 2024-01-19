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

format_pval <- function(x){ return(format(x, scientific=TRUE, digits=3)) }

wideScreen <- function(howWide=Sys.getenv("COLUMNS")) { options(width=as.integer(howWide)) }

strsplits <- function(x, splits, ...){
	for (split in splits) { x <- unlist(strsplit(x, split, ...)) }
	return(x[!x == ""]) # Remove empty values
}

## functions for GWAS_CATALOG & previous TWAS comparison

# near snp
near_gwas_snp <- function(row, gwas_assoc_locs){
	chrom <- as.integer(row[['Chrom']]); wstart <- as.integer(row[['w_start']]); wend <- as.integer(row[['w_end']])
	chrom_gwas_assoc_locs <- gwas_assoc_locs[as.integer(gwas_assoc_locs[['Chrom']]) == chrom, 'loc']
	return(as.integer(any(chrom_gwas_assoc_locs >= wstart & chrom_gwas_assoc_locs <= wend)))
}

near_gwas_gene <- function(row, gwas_data){
	chrom <- as.integer(row[['Chrom']]); wstart <- as.integer(row[['w_start']]); wend <- as.integer(row[['w_end']])
	locs <- gwas_data[gwas_data[['Chrom']] == chrom & gwas_data[['prev_id_gwas']] == 1, ]
	if (nrow(locs) == 0) { return(0) }
	locs <- locs[!duplicated(locs),]
	locs <- locs[,c('Start', 'End')]
	return(as.integer(any(locs >= wstart & locs <= wend)))
}

# near either
near_gwas <- function(row, gwas_data){
	out <- (near_gwas_gene(row, gwas_data) + near_gwas_snp(row)) > 0
	return(as.integer(out))
}

# near gene; (search whole dat dataframe, but only for rows in sig_data)
near_twas <- function(row, twas_data){
	chrom <- as.integer(row[['Chrom']]); wstart <- as.integer(row[['w_start']]);  wend <- as.integer(row[['w_end']])
	locs <- twas_data[twas_data[['Chrom']] == chrom, ]
	if (length(locs) == 0) { return(0) }
	locs <- locs[!duplicated(locs),]
	locs <- locs[,c('Start', 'End')]

	return(as.integer(any(locs >= wstart & locs <= wend)))
}



## FUNCTIONS BASED ON GENE LOCATION AND SIGNIFICANCE
# function to find most significant genes in 1MB window
most_sig_no_overlap <- function(dataset, data){
	table_data <- data[data$dataset==dataset & data$Pvalue_raw <= 2.5e-6, c('Gene', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue')]
	table_data <- na.omit(table_data)
	table_data <- table_data[order(table_data$Pvalue), ]

	gr_dat <- table_data[,c('Gene','Chrom','Start','End')]
	colnames(gr_dat) <- c('name','chr','start','end')
	gr <- makeGRangesFromDataFrame(gr_dat, keep.extra.columns=TRUE)
	gr_wide <- resize(gr, width=2000000+width(gr), fix='center')
	overlaps <- findOverlaps(gr_wide, type='any', select='first')
	# overlaps <- findOverlaps(gr, maxgap=1000000, type = 'any', select = 'first')
	most_sig_names <- as.data.frame(gr[!duplicated(overlaps)])[['name']]
	return(most_sig_names)
}

get_overlap_locus_group <- function(data){
	table_data <- data[data$Pvalue_raw <= 2.5e-6, c('Gene', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue')]
	table_data <- na.omit(table_data)
	table_data <- table_data[order(table_data$Pvalue),]

	gr_dat <- table_data[,c('Gene','Chrom','Start','End')]
	colnames(gr_dat) <- c('name','chr','start','end')
	gr <- makeGRangesFromDataFrame(gr_dat, keep.extra.columns=TRUE)
	gr_wide <- resize(gr, width = 2000000 + width(gr), fix = 'center')
	overlaps <- findOverlaps(gr_wide, type = 'any', select = 'first')
	# overlaps <- findOverlaps(gr, maxgap=1000000, type = 'any', select = 'first')
	table_data$overlaps <- overlaps
	table_data$overlaps_rank <- with(table_data, ave(seq(overlaps), overlaps, FUN = function(x) seq(length(x))))
	return(overlaps)
}

get_overlap_locus_group_rank <- function(data){
	table_data <- data[data$Pvalue_raw <= 2.5e-6, c('Gene', 'Chrom', 'Start', 'End', 'Zscore', 'Pvalue')]
	table_data <- na.omit(table_data)
	table_data <- table_data[order(table_data$Pvalue),]

	gr_dat <- table_data[,c('Gene','Chrom','Start','End')]
	colnames(gr_dat) <- c('name','chr','start','end')
	gr <- makeGRangesFromDataFrame(gr_dat, keep.extra.columns=TRUE)
	gr_wide <- resize(gr, width = 2000000 + width(gr), fix = 'center')
	overlaps <- findOverlaps(gr_wide, type = 'any', select = 'first')
	# overlaps <- findOverlaps(gr, maxgap=1000000, type = 'any', select = 'first')
	table_data$overlaps <- overlaps
	table_data$overlaps_rank <- with(table_data, ave(seq(overlaps), overlaps, FUN = function(x) seq(length(x))))
	return(table_data$overlaps_rank)
}

## annotation functions
annot_name <- function(row){
	search <- c('prev_id_gwas', 'near_prev_id_gwas', 'prev_id_twas', 'near_prev_id_twas')
	# gwas
	if (row[[search[1]]] == 1) { note <- letters[1] # prev id
		} else if (row[[search[2]]] == 1) { note <- letters[2] # near
		} else { note <- '' }
	# twas
	if (row[[search[3]]] == 1) { note <- paste0(note, letters[3]) # prev id
		} else if (row[[search[4]]] == 1) { note <- paste0(note, letters[4]) # near
		}
	if (note != '') { note <- paste0('$^{', note, '}$') } # add math notation
	return(paste0(row[['Gene']], note)) }

annot_letter <- function(row){
	search <- c('prev_id_gwas', 'near_prev_id_gwas', 'prev_id_twas', 'near_prev_id_twas')
	# gwas
	if (row[[search[1]]] == 1) { note <- letters[1] # prev id
		} else if (row[[search[2]]] == 1) { note <- letters[2] # near
		} else { note <- '' }
	# twas
	if (row[[search[3]]] == 1) { note <- paste0(note, letters[3]) # prev id
		} else if (row[[search[4]]] == 1) { note <- paste0(note, letters[4]) # near
		}
	if (note == '') { note <- 0 }
	return(note)}

########################
# unused
prev_id_gwas_overlap <- function(dataset, data, select='first'){

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
		prev_id_gwas_matches <- as.data.frame(cbind(by(odf$subjectHits, odf$queryHits, paste, collapse=', ')))
		prev_id_gwas_matches$name <- row.names(prev_id_gwas_matches)
		colnames(prev_id_gwas_matches) <- c('prev_id_gwas_overlap', 'name')
		row.names(prev_id_gwas_matches) <- NULL

		gr_dat_overlapped_previd <- merge(gr_dat, prev_id_gwas_matches, all = TRUE)
		row.names(gr_dat_overlapped_previd) <- NULL
	}

	# return genes that are near previd'd, with column denoting which gene they're near
	return(gr_dat_overlapped_previd)
}

################################################################



