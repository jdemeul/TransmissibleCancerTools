## reading allelecounts
# check 100% DFT and 100% fibroblast first
library(ggplot2)
library(GenomicRanges)
library(readr)

devilinfo_file <- "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/reference/Devil7.1.fa.fai"

tumour_allelecounts_file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/Fibroblast_91H_DFT2_202T2_50-50MIX_DNA_alleleCount.txt"
host_allelecounts_file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/Fibroblast_91H_cell_line_DNA_alleleCount.txt"
snploci_file <- "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/reference/Devil_SNPList.txt"


# dft2_dna <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/DFT2_202T2_cell_line_DNA_alleleCount.txt", as.is = T)
# fibro_dna <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/Fibroblast_91H_cell_line_DNA_alleleCount.txt", as.is = T)
# 
# # fake 50-50 mix tumour normal
# tumour_dna <- cbind(dft2_dna[ , 1:2], round(32*(dft2_dna[ , 3:7]/median(dft2_dna$Good_depth) + fibro_dna[ , 3:7]/median(fibro_dna$Good_depth)) ))
# colnames(tumour_dna) <- c('#CHR',"POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")
# write.table(x = tumour_dna, file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/Fibroblast_91H_DFT2_202T2_50-50MIX_DNA_alleleCount.txt", sep = "\t", quote = F, row.names = F)

setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/")

source("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/TransmissibleCancerTools/CN_calling/transmissibleASCAT.R")

devilinfo_gr <- get_devil_genome(devilinfofile = devilinfo_file)

# debug(get_baf_logr)
baflogr <- get_baf_logr(tumour_allelecounts_file = tumour_allelecounts_file,
                        host_allelecounts_file = host_allelecounts_file,
                        reference_alleles_file = snploci_file,
                        genomeinfo = devilinfo_gr,
                        min_contig_size = 1e6,
                        min_depth = 10)


## adding the "position along the chromosome"
# snploci$REF <- factor(snploci$REF, levels = c("A", "C", "G", "T"))
# snploci$ALT <- factor(snploci$ALT, levels = c("A", "C", "G", "T"))

### quick check of total CN reconstruction
check_rho_psi_estimates(sampleid = "Fibro_DFT2_50-50MIX", baflogr = baflogr, rho = .5, psit = 1.85)

# debug(genotype_loci)
genotype_loci(baflogr = baflogr, rho = .5, psit = 2)


findpeaks <- function (x, thresh = 0.00001, from = -.75, to = 1.75, mindens = .25) {
  xdens <- stats::density(x, kernel = "gaussian", from = from, to = to)
  pks <- which(diff(sign(diff(xdens$y, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  pks <- pks[xdens$y[pks - 1] - xdens$y[pks] > thresh & xdens$y[pks] > mindens]
  return(xdens$x[pks])
}

nminpeaks <- findpeaks(x = baflogr$nmin, thresh = .00001)
nminonepeak <- nminpeaks[which.min(abs(nminpeaks - 1))]
nminzeropeak <- nminpeaks[which.min(abs(nminpeaks - 0))]
# ntotpeaks <- findpeaks(x = baflogr$ntotsm, thresh = .00001, from = .25, to =1.75)

cnonesegs <- rle(x = baflogr$ntotsm <= nminonepeak + 1/3)
cnonesegs$values[cnonesegs$lengths < 100] <- F
cnonesegs <- inverse.rle(cnonesegs)

segstarts <- which(c(cnonesegs, F) & !c(F, cnonesegs))
segends <- which(!c(cnonesegs, F) & c(F, cnonesegs))


nminzerosegs <- rle(x = baflogr$nminsm <= nminzeropeak + nminonepeak/50)
nminzerosegs$values[nminzerosegs$lengths < 1000] <- F
nminzerosegs <- inverse.rle(nminzerosegs)

segstarts2 <- which(c(nminzerosegs, F) & !c(F, nminzerosegs))
segends2 <- which(!c(nminzerosegs, F) & c(F, nminzerosegs))


### selecting SNPs to run on:
# get all hetloci (nmin ≥ 1) which are het in host as well (so htgeno = AB/*)
hetidxs <- which(baflogr$htgeno == "AB/*")
# nminoneidxs <- which(baflogr$nmin >= (.75*nminonepeak + .25*nminzeropeak)/2)
nminoneidxs <- which(baflogr$nmin >= (.66*nminonepeak + .33*nminzeropeak))
cnoneidxs <- unlist(mapply(segstarts, segends, FUN = ':'))
nminzeroidxs <- unlist(mapply(segstarts2, segends2, FUN = ':'))
hostAAidxs <- which(baflogr$htgeno == 'AA/B*')
hostBBidxs <- which(baflogr$htgeno == 'BB/A*')

# finalidxs <- intersect(hetidxs, c(nminoneidxs, cnoneidxs, nminzeroidxs))

baflogr$type <- NA
baflogr[intersect(hetidxs, nminoneidxs), "type"] <- "AB/AB"
baflogr[intersect(hetidxs, c(cnoneidxs, nminzeroidxs)), "type"] <- 'AB/A'
baflogr[intersect(hostAAidxs, nminoneidxs), "type"] <- 'AA/AB'
baflogr[intersect(hostBBidxs, nminoneidxs), "type"] <- 'BB/AB'

# finalbaflogr <- baflogr[finalidxs, ]

# p1 <- ggplot(data = baflogr[baflogr$type %in% c("AB/AB", "AB/A") & baflogr$chr != "Chrx", ], mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = baf), colour = "goldenrod2")
p1 <- ggplot(data = baflogr[baflogr$type %in% c("AB/AB", "AB/A") & baflogr$chr != "X", ], mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = baf), colour = "goldenrod2")
p1 <- p1 + facet_wrap(~chr)
# p1 <- p1+ylim(c(-.2, 6))
p1
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180426_input_BAF_first_ASCAT_run.png", plot = p1, width = 16, height = 6)


# write.table(x = baflogr, file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180425_selectedSNPs_counts.txt", quote = F, sep = "\t")
write_tsv(x = baflogr, path = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180508_selectedSNPs_counts.txt", col_names = T)
# write.table(x = baflogr, file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180508_selectedSNPs_counts.txt", quote = F, sep = "\t", row.names = F)
