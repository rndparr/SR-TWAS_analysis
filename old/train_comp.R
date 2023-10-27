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


train_dir='/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/GTEx_6tissue/'
out_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/plots/'

# naive_path='GTEx6tissue_Naive_train_GTEx_LD_GTEx_train_results.txt'
# sr_path='GTEx6tissue_SR_train_GTEx_LD_GTEx_train_results.txt'


naive_train_dat <- read.table(
	paste0(train_dir, 'GTEx6tissue_Naive_train_GTEx_LD_GTEx_train_results.txt'), 
	header=TRUE, sep='\t')


sr_train_dat <- read.table(
	paste0(train_dir, 'GTEx6tissue_SR_train_GTEx_LD_GTEx_train_results.txt'), 
	header=TRUE, sep='\t')


## RENAME COLUMNS
colnames(naive_train_dat)[which(colnames(naive_train_dat) == 'CHROM')] <- 'Chrom'
colnames(naive_train_dat)[which(colnames(naive_train_dat) == 'GeneStart')] <- 'Start'
colnames(naive_train_dat)[which(colnames(naive_train_dat) == 'GeneEnd')] <- 'End'
colnames(naive_train_dat)[which(colnames(naive_train_dat) == 'GeneName')] <- 'Gene'


# read in dat from manplots, filter down etc
naive_twas_dat <- dat[dat$dataset == 'Naive', -1]
naive_dat <- merge(naive_twas_dat, naive_train_dat)


# figure out which rows have only 1 base model
cols1 <- colnames(naive_dat)[c(1:5, which(sapply(colnames(naive_dat), endsWith, 'NSNP')))]
naive_base_model_comp <- bleh
naive_base_model_comp <- naive_dat[, cols1]
naive_base_model_comp[,6:11] <- (naive_base_model_comp[, 6:11] > 0)*1
naive_base_model_comp$N_base_models <- rowSums(naive_base_model_comp[,6:11])

head(naive_base_model_comp,20)

sum((naive_base_model_comp$N_base_models == 1) * 1)
# 228 of the output genes of 22825
# ~1%

nrow(dat[dat$dataset == 'SR-TWAS', -1])
# PDTWAS, AD-TWAS2
# 22921

# check if any only 1 model were significant 
naive_base_model_comp[which(naive_base_model_comp$N_base_models==1),]

sig_naive_targets <- dat[(dat$dataset == 'Naive') & (dat$Pvalue < 2.5e-6), 'TargetID']
sig_sr_targets <- dat[(dat$dataset == 'SR-TWAS') & (dat$Pvalue < 2.5e-6), 'TargetID']

# any overlap between naive sig targets and single base model?
naive_base_model_comp[ (naive_base_model_comp$TargetID %in% sig_naive_targets) & (naive_base_model_comp$N_base_models == 1), ]
# PDTWAS:
# nope
# ADTWAS2
# 1 overlap: ENSG00000112195.8 TREML2

# any overlap between SR sig targets and single base model?
naive_base_model_comp[ (naive_base_model_comp$TargetID %in% sig_sr_targets) & (naive_base_model_comp$N_base_models == 1), ]
# PDTWAS:
# nope
# ADTWAS2
# 1 overlap: ENSG00000112195.8 TREML2

naive_base_model_comp[ (naive_base_model_comp$TargetID %in% sig_naive_targets), ]





####### SR-TWAS


sr_train_dat <- read.table(
	paste0(train_dir, 'GTEx6tissue_SR_train_GTEx_LD_GTEx_train_results.txt'), 
	header=TRUE, sep='\t')


## RENAME COLUMNS
colnames(sr_train_dat)[which(colnames(sr_train_dat) == 'CHROM')] <- 'Chrom'
colnames(sr_train_dat)[which(colnames(sr_train_dat) == 'GeneStart')] <- 'Start'
colnames(sr_train_dat)[which(colnames(sr_train_dat) == 'GeneEnd')] <- 'End'
colnames(sr_train_dat)[which(colnames(sr_train_dat) == 'GeneName')] <- 'Gene'


# read in dat from manplots, filter down etc
sr_twas_dat <- dat[dat$dataset == 'SR-TWAS', -1]
sr_dat <- merge(sr_twas_dat, sr_train_dat)


# figure out which rows have only 1 base model
cols2 <- colnames(sr_dat)[c(1:5, which(sapply(colnames(sr_dat), startsWith, 'Z_')))]
sr_base_model_comp <- sr_dat[, cols2]
# exact equal to 0
sr_base_model_comp[,6:11] <- (sr_base_model_comp[, 6:11] > 0)*1
sr_base_model_comp$N_sr_base_models <- rowSums(sr_base_model_comp[,6:11])
head(sr_base_model_comp)

head(merge(naive_base_model_comp[,-c(6:11)], sr_base_model_comp[,-c(6:11)]))
naive_sr_comp <- merge(naive_base_model_comp[,-c(6:11)], sr_base_model_comp[,-c(6:11)])

# all models actually used the same number of basemodels; ie, sr-twas didn't exclude any base models except those that did not have a model, ie did not set any Z to 0 unless that base model didn't have a GReX prediction model for that gene
all(naive_sr_comp$N_sr_base_models == naive_sr_comp$N_base_models)



sr_base_model_comp <- sr_dat[, cols2]
head(sr_base_model_comp)
sr_base_model_comp <- format(sr_base_model_comp, digits = 3)


################
# hmm. looks like should look at Z_tissue values that are approx equal to 0
approxnonzero <- function(x, tol = 1e-12) {
	abs(x) >= tol
}

valapproxnonzero <- function(x, tol = 1e-12) {
	ifelse(abs(x) >= tol, x, 0)
}

# figure out which rows have only 1 base model
cols2 <- colnames(sr_dat)[c(1:5, which(sapply(colnames(sr_dat), startsWith, 'Z_')))]
sr_base_model_comp1 <- sr_dat[, cols2]
# approx nonequal to 0
sr_base_model_comp <- sr_base_model_comp1
sr_base_model_comp[,6:11] <- (approxnonzero(sr_base_model_comp1[, 6:11]))*1

sr_base_model_comp$N_sr_base_models <- rowSums(sr_base_model_comp[,6:11])
head(sr_base_model_comp)

head(merge(naive_base_model_comp[,-c(6:11)], sr_base_model_comp[,-c(6:11)]))
naive_sr_comp <- merge(naive_base_model_comp[,-c(6:11)], sr_base_model_comp[,-c(6:11)])

# all models actually used the same number of basemodels; ie, sr-twas didn't exclude any base models except those that did not have a model, ie did not set any Z to 0 unless that base model didn't have a GReX prediction model for that gene
all(naive_sr_comp$N_sr_base_models == naive_sr_comp$N_base_models)
## FALSE

naive_sr_comp$N_base_model_diff <- naive_sr_comp$N_base_models - naive_sr_comp$N_sr_base_models
head(naive_sr_comp)

# naive_sr_comp[, -c(6,7)]

# number of genes where same number of base models actually used:
sum(naive_sr_comp$N_base_model_diff == 0)

nrow(naive_sr_comp)
# total genes: 22831

# with tol of 1e-16
# 4528
# 4528/22831 =  0.1983268
# about 1/5

# when tol set to 1e-12
# 2608
# 2608/22831 =  0.1142307
# about 1/10



#### compare with sig dat
sig_sr_targets <- dat[(dat$dataset == 'SR-TWAS') & (dat$Pvalue < 2.5e-6), 'TargetID']

# number of sig targets where naive, sr have same number of base models
naive_sr_comp[ (naive_sr_comp$TargetID %in% sig_sr_targets) & (naive_sr_comp$N_base_model_diff == 0), ]
x <- nrow(naive_sr_comp[ (naive_sr_comp$TargetID %in% sig_sr_targets) & (naive_sr_comp$N_base_model_diff == 0), ])
x
targets <- naive_sr_comp[ (naive_sr_comp$TargetID %in% sig_sr_targets) & (naive_sr_comp$N_base_model_diff == 0), 'TargetID']
# PDTWAS:
# tol 1e-12
# 6
# ADTWAS2
# tol 1e-16
# 15
# tol 1e-12
# 6


# sr_base_model_comp1[sr_base_model_comp1$TargetID %in% targets, ]


# # any overlap between SR sig targets and single base model?
# naive_base_model_comp[ (naive_base_model_comp$TargetID %in% sig_sr_targets) & (naive_base_model_comp$N_base_models == 1), ]
# # PDTWAS:
# # 
# # ADTWAS2
# # 

# see if can find number of genes where the value of each non-zero Z is equal to approximately 1/(# of used base models)
cols2 <- colnames(sr_dat)[c(1:5, which(sapply(colnames(sr_dat), startsWith, 'Z_')))]
sr_base_model_comp1 <- sr_dat[, cols2]
# approx nonequal to 0
sr_base_model_comp <- sr_base_model_comp1
sr_base_model_comp$N_sr_base_models <- rowSums((approxnonzero(sr_base_model_comp1[, 6:11]))*1)

bleh <- merge(sr_base_model_comp, naive_base_model_comp[, -c(6:11)])


head(bleh)
bleh$mean_Z <- 1 / bleh$N_base_models
bleh$mean_Z_sr <- 1 / bleh$N_sr_base_models


# GENES WHERE Zs ARE EQUAL TO 1/NUMBER BASE MODELS TOTAL
bleh2 <- bleh
# replace approx 0 intries with NA
bleh2[,c(6:11)] <- replace(bleh2[,c(6:11)], !approxnonzero(bleh2[,c(6:11)]), NA)
bleh2[,c(6:11)] <- (bleh2[,c(6:11)] - bleh2$mean_Z)
bleh2$total_mean_diff <- rowSums(abs(bleh2[,c(6:11)]), na.rm=TRUE)
bleh2$total_mean_diff_approx_0 <- (!approxnonzero(bleh2$total_mean_diff))*1

head(bleh2[-c(1:3,12,14)], 20)

sum(bleh2$total_mean_diff_approx_0)
# training:
# 871

# genes
bleh2_1 <- bleh2[which(bleh2$total_mean_diff_approx_0 == 1), -c(1:3,6:11,12,14)]
head(bleh2_1, 20)

bleh2_1[bleh2_1$TargetID %in% sig_sr_targets, ]
targets <- bleh2_1$TargetID[bleh2_1$TargetID %in% sig_sr_targets]
bleh[bleh$TargetID %in% targets,]
# ADTWAS2 approx 4 genes where SR result was approximately the mean of the used base models
#       Chrom    Start      End           TargetID              Gene
# 3948     11 60455752 60470760 ENSG00000156738.17             MS4A1
# 11731    19 44846175 44889228  ENSG00000130202.9           NECTIN2
# 19441     6 32255711 32265838  ENSG00000225914.1 XXbac-BPG154L12.4
# 19590     6 41190277 41201194  ENSG00000112195.8            TREML2
# PDTWAS approx 2 genes where SR result was approximately the mean of the used base models
#      Chrom    Start      End           TargetID        Gene
# 9656    17 45249534 45268630  ENSG00000267278.5 MAP3K14-AS1
# 9673    17 45894382 46028334 ENSG00000186868.15        MAPT

# GENES WHERE Zs ARE EQUAL TO 1/NUMBER BASE MODELS USED IN SR
bleh3 <- bleh

# replace approx 0 intries with NA
bleh3[,c(6:11)] <- replace(bleh3[,c(6:11)], !approxnonzero(bleh3[,c(6:11)]), NA)
bleh3[,c(6:11)] <- bleh3[,c(6:11)] - bleh3$mean_Z_sr
bleh3$total_mean_diff <- rowSums(abs(bleh3[,c(6:11)]), na.rm=TRUE)
bleh3$total_mean_diff_approx_0 <- (!approxnonzero(bleh3$total_mean_diff))*1

head(bleh3[,-c(1:3,13,15)], 20)

sum(bleh3$total_mean_diff_approx_0)
# training
# 1526

bleh3_2 <- bleh3[bleh3$total_mean_diff > 0,]
head(bleh3_2[order(bleh3_2$total_mean_diff),])

bleh3_3 <- bleh3[approxnonzero(bleh3$total_mean_diff),]
head(bleh3_3[order(bleh3_3$total_mean_diff),], 50)

bleh3_4 <- merge(bleh3_3[,-c(6:11)],bleh[bleh$TargetID %in% bleh3_3$TargetID,])
head(bleh3_4)


# genes
bleh3_1 <- bleh3[which(bleh3$total_mean_diff_approx_0 == 1), -c(1:3,6:11,13,15)]
head(bleh3_1, 20)

bleh3_1[bleh3_1$TargetID %in% sig_sr_targets, ]
targets <- bleh3_1$TargetID[bleh3_1$TargetID %in% sig_sr_targets]
bleh[bleh$TargetID %in% targets,]
# ADTWAS2 approx 4 genes where SR result was approximately the mean of the used base models
#       Chrom    Start      End           TargetID              Gene
# 3948     11 60455752 60470760 ENSG00000156738.17             MS4A1
# 11731    19 44846175 44889228  ENSG00000130202.9           NECTIN2
# 19441     6 32255711 32265838  ENSG00000225914.1 XXbac-BPG154L12.4
# 19590     6 41190277 41201194  ENSG00000112195.8            TREML2
# PDTWAS approx 5 genes where SR result was approximately the mean of the used base models
#      Chrom    Start      End           TargetID        Gene
# 8293    16 30761198 30762710 ENSG00000196118.11     CCDC189
# 8318    16 31090842 31095980 ENSG00000167397.14      VKORC1
# 9653    17 45222223 45247266 ENSG00000184922.13       FMNL1
# 9656    17 45249534 45268630  ENSG00000267278.5 MAP3K14-AS1
# 9673    17 45894382 46028334 ENSG00000186868.15        MAPT


targetsz <- by(dat, dat$dataset, function(x) {most_sig_no_overlap2(x) } )
targetsz
# naive
sort(targetsz[[1]])
# sr
sort(targetsz[[2]])

## sr most sig no overlap
sr_most_sig_no_overlap <- dat[(dat$Gene %in% targetsz[[2]]) & (dat$dataset == 'SR-TWAS'), c(2:6,12)]
sr_most_sig_no_overlap
nrow(sr_most_sig_no_overlap)
# 11

## naive most sig no overlap
naive_most_sig_no_overlap <- dat[(dat$Gene %in% targetsz[[1]]) & (dat$dataset == 'Naive'), c(2:6,12)]
naive_most_sig_no_overlap
nrow(naive_most_sig_no_overlap)
# 12


## sr most sig no overlap unique
targetsz_overlap <- intersect(targetsz[[1]], targetsz[[2]])
sr_most_sig_no_overlap_unique <- dat[(dat$Gene %in% setdiff(targetsz[[2]], targetsz_overlap)) & (dat$dataset == 'SR-TWAS'), c(2:6,12)]
sr_most_sig_no_overlap_unique
nrow(sr_most_sig_no_overlap_unique)
# 5

## naive most sig no overlap unique
naive_most_sig_no_overlap_unique <- dat[(dat$Gene %in% setdiff(targetsz[[1]], targetsz_overlap)) & (dat$dataset == 'Naive'), c(2:6,12)]
naive_most_sig_no_overlap_unique
nrow(naive_most_sig_no_overlap_unique)
# 6



sr_sig_z_comp <- merge(dat[(dat$Gene %in% targetsz[[2]]) & (dat$dataset == 'SR-TWAS'),], bleh)
