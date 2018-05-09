#### running new CN pipeline on examples from dogs

library(ggplot2)
library(readr)
library(GenomicRanges)

setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/")

readcounts_file <- "../data/matched_pair_559_readcounts.tsv.gz"

readcounts <- read_tsv(file = readcounts_file, col_names = T, col_types = "cicciiiiii", progress = T)
baflogr <- as.data.frame(readcounts[, -c(7,10)])
colnames(baflogr) <- c("chr", "cumpos", "REF", "ALT", "h_ref", "h_alt", "t_ref", "t_alt")

## some filtering
mincounts <- 10
has_good_depth <- baflogr$h_ref + baflogr$h_alt >= mincounts
has_single_read_tumour <- baflogr$t_ref + baflogr$t_alt >= 1

baflogr <- baflogr[has_good_depth & has_single_read_tumour, ]

### quick check of total CN reconstruction
rho <- .76
psit <- 2



nminonepeak <- 1
nminzeropeak <- 0
