#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE)


###############
# load libraries
library(ggplot2)
library(reshape2)


######################################################
# directory for simulation files
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/'

dat_raw <- read.table(paste0(dir, 'results/pred_sig_results.txt'), header=FALSE, colClasses=c('character', rep('numeric',4)))
colnames(dat_raw) <- c('i', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')

dat <- melt(dat_raw, measure.vars=c('ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN'),
	variable.name='model', value.name='pvalue', na.rm=TRUE)

##############


# dat$o = -log10(sort(dat$pvalue, decreasing=FALSE))
# dat$e = -log10( 1:length(dat$o)/length(dat$o) )
# max_lim <- max(c(dat$e, dat$o[dat$o<Inf]))

# pdat <- data.frame(e, o)


p <- ggplot(dat, aes(sample=pvalue)) + 
	geom_qq(distribution=stats::qchisq, dparams=list(df=1)) +
	geom_qq_line(color='red') +
	labs(x='Expected', y='Observed') +
	facet_wrap(model ~ .)


ggsave(paste0(dir, 'qqplot', '.pdf'), p)
ggsave(paste0(dir, 'qqplot', '.png'), p)


# pdf(file=dpr_man_path)
# grid.arrange(
# 	do_man_plot('DPR','Breast', drop_most_sig=TRUE, id_names=id_names) + ggtitle('A'), 
# 	do_man_plot('DPR','Ovary', id_names=id_names) + ggtitle('B'),
# 	nrow = 2)
# dev.off()


# ##############

# make_qq_plot <- function(pvector, size=8) {
#     pvector = pvector[!is.na(pvector)]
#     o = -log10(sort(pvector, decreasing=FALSE))
#     e = -log10( 1:length(o)/length(o) )
#     max_lim <- max(c(e, o[o<Inf]))
#     qplot(e, o) + xlim(c(0, max_lim)) + ylim(c(0, max_lim)) + 
#         labs(x = "Expected", y = "Observed") + 
#         geom_abline(intercept = 0, slope = 1, colour = "red") + 
#     theme(text = element_text(size = size))
# }

# make_qq_plot <- function(pvector, size=8) {


# 	p <- ggplot( ) +
# 		geom_point()

# }

