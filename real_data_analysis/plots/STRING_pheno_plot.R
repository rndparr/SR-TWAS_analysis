#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE)

###############
# load libraries
library(ggplot2)

wkdir <- '~/pCloud Drive/Documents/Important/YangLab/SR_TWAS_manuscript/STRING/'


# import data
get_dat <- function(twas, twidth=15){
		dat <- read.table(paste0(wkdir, twas, '_enrichment.HPO.tsv'), header=TRUE, sep='\t', comment.char='')

	if(twas=='AD'){
		dat[dat$term.description=='Family history of Alzheimer\342\200\231s disease','term.description'] <- "Family history of Alzheimer's disease"


		dat <- dat[ !(dat$term.description %in% c('Nervous system disease', 'Serum alanine aminotransferase measurement', 'Bone measurement', 'Disease')), ]
	}

	dat <- dat[order(-dat$false.discovery.rate), ]
	row.names(dat) <- NULL

	dat$term <- factor(dat$term.description, levels=dat$term.description, labels=stringr::str_wrap(dat$term.description, twidth) )

	pal <- c('#3A948E','#36807A','#2f615d','#44B5AD')
	text_size <- 11

	plot_colors <- rep(pal, (nrow(dat) %/% 4) + as.integer((nrow(dat) %% 4) != 0))[1:nrow(dat)]

	dat$twas <- list('AD'='AD-TWAS', 'PD'='PD-TWAS')[[twas]]

	dat$neglogFDR <- -log10(dat$false.discovery.rate)

	return(dat)
}


pdat <- get_dat('PD')
adat <- get_dat('AD')

get_plot <- function(twas){

	dat <- get_dat(twas, twidth=25)
	pal <- c('#3A948E','#36807A','#2f615d','#44B5AD')

	if (twas == 'AD'){
		text_size <- 16
	} else {
		text_size <- 14
	}

	plot_colors <- rep(pal, (nrow(dat) %/% 4) + as.integer((nrow(dat) %% 4) != 0))[1:nrow(dat)]

	p <- ggplot(data=dat, aes(x=term, y=neglogFDR, color=term, fill=term)) + 
		geom_col(alpha=0.95) +
		labs(y=bquote(bold(-"log"["10"]("FDR"))), x='') +
		scale_color_manual(values=plot_colors) +
		scale_fill_manual(values=plot_colors) +
		theme_bw() +
		guides(color='none', fill='none') +
		theme(text=element_text(size=text_size+1),
			axis.text.y=element_text(size=text_size+1, color='black', face='bold'),
			axis.text.x=element_text(size=text_size+1, color='black', face='bold'),
			axis.title.x=element_text(size=text_size+3, color='black', face='bold'),
			# panel.grid.major.x=element_blank(),
			# panel.grid.minor.x=element_blank(),
			panel.grid.major.y=element_blank(),
			panel.grid.minor.y=element_blank(),
			plot.tag=element_text(face='bold', size=text_size+5, family='Helvetica')) +
		coord_flip()

	return(p)
}


adp <- get_plot('AD')
pdp <- get_plot('PD')



file.remove(paste0(wkdir, 'AD_phenotype_bar.pdf'))
ggsave(paste0(wkdir, 'AD_phenotype_bar.pdf'), adp, width=11.5, height=11)


file.remove(paste0(wkdir, 'PD_phenotype_bar.pdf'))
ggsave(paste0(wkdir, 'PD_phenotype_bar.pdf'), pdp, width=10, height=7)


