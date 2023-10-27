
options(stringsAsFactors=FALSE)


### TWAS to do
twas <- 'AD_TWAS'
# twas <- 'PD_TWAS'
twas <- 'AD_TWAS2'

### directories
## HGCC
in_dir <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/', twas, '/')


# suffix <- ''
suffix <- '_4basemodels'
load(paste0(in_dir, 'overlap_table_dat_', twas, suffix, '.Rdata'))


sig_dat[which(sig_dat$most_sig_SR_ROSMAP == 1),]

srgenes <- unique(sig_dat[which(sig_dat$most_sig_SR_ROSMAP == 1),'Gene'])

naivegenes <- unique(sig_dat[which(sig_dat$most_sig_Naive_ROSMAP == 1),'Gene'])

setdiff(srgenes, naivegenes)

setdiff(naivegenes,srgenes)

sig_dat[sig_dat$Gene %in% setdiff(naivegenes,srgenes),]





sig_dat[which(sig_dat$most_sig_SR_ROSMAP == 1 & sig_dat$dataset %in% c('Naive_ROSMAP','SR_ROSMAP')),]



dput(sort(sig_dat[which(sig_dat$most_sig_SR_ROSMAP == 1 & sig_dat$dataset =='SR_ROSMAP'),'Gene']))