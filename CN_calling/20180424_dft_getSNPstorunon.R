## reading allelecounts
# check 100% DFT and 100% fibroblast first
library(ggplot2)
library(GenomicRanges)

setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/")
devilinfo_file <- "../reference/Devil7.1.fa.fai"
snploci_file <- "../reference/Devil_SNPList.txt"

dft2_dna_file <- "../results/DFT2_202T2_cell_line_DNA_alleleCount.txt"
fibro_dna_file <- "../results/Fibroblast_91H_cell_line_DNA_alleleCount.txt"


devilinfo <- read.delim(file = devilinfo_file, as.is = T, header = F)
colnames(devilinfo) <- c("contig", "width")
devilinfo$chr <- substr(x = devilinfo$contig, start = 1, stop = 4)
devilinfo$start <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = function(x) cumsum(c(1,x))[-length(x)]))
devilinfo$end <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = cumsum))
devilinfo_gr <- GRanges(seqnames = devilinfo$chr, ranges = IRanges(start = devilinfo$start, end = devilinfo$end), mcols = DataFrame(contig = devilinfo$contig))

snploci <- read.delim(file = snploci_file, as.is = T)

dft2_dna <- read.delim(file = dft2_dna_file, as.is = T)
fibro_dna <- read.delim(file = fibro_dna_file, as.is = T)

# fake 50-50 mix tumour normal
tumour_dna <- cbind(dft2_dna[ , c("X.CHR", "POS")], dft2_dna[ , 3:6] + fibro_dna[ , 3:6])


## adding the "position along the chromosome"
# snploci$REF <- factor(snploci$REF, levels = c("A", "C", "G", "T"))
# snploci$ALT <- factor(snploci$ALT, levels = c("A", "C", "G", "T"))

baflogr <- snploci[, -5]
baflogr$chr <- substr(baflogr$Scaffold, 1, 4)
baflogr$cumpos <- start(devilinfo_gr[match(x = baflogr$Scaffold, table = mcols(devilinfo_gr)$mcols.contig)]) + baflogr$Position
baflogr$h_ref <- ifelse(snploci$REF == "A", fibro_dna$Count_A,
                        ifelse(snploci$REF == "C", fibro_dna$Count_C, 
                               ifelse(snploci$REF == "G", fibro_dna$Count_G, fibro_dna$Count_T)))
baflogr$h_alt <- ifelse(snploci$ALT == "A", fibro_dna$Count_A,
                        ifelse(snploci$ALT == "C", fibro_dna$Count_C, 
                               ifelse(snploci$ALT == "G", fibro_dna$Count_G, fibro_dna$Count_T)))

baflogr$t_ref <- ifelse(snploci$REF == "A", tumour_dna$Count_A,
                        ifelse(snploci$REF == "C", tumour_dna$Count_C, 
                               ifelse(snploci$REF == "G", tumour_dna$Count_G, tumour_dna$Count_T)))
baflogr$t_alt <- ifelse(snploci$ALT == "A", tumour_dna$Count_A,
                        ifelse(snploci$ALT == "C", tumour_dna$Count_C, 
                               ifelse(snploci$ALT == "G", tumour_dna$Count_G, tumour_dna$Count_T)))


## some filtering
mincounts <- 10
mincontiglength <- 1e6

allowed_contigs <- devilinfo[devilinfo$width >= 1e6, "contig"]
has_good_depth <- baflogr$h_ref + baflogr$h_alt >= mincounts
# on_sel_chrom <- dft2_dna$chr == "Chr4"

baflogr <- baflogr[has_good_depth & baflogr$Scaffold %in% allowed_contigs, ]
# dft2_dna <- dft2_dna[has_good_depth & dft2_dna$X.CHR %in% allowed_contigs, ]
# fibro_dna <- fibro_dna[has_good_depth & fibro_dna$X.CHR %in% allowed_contigs, ]


### quick check of total CN reconstruction
rho <-  .5
psit <- 2.1

totalhost <- baflogr$h_ref + baflogr$h_alt
totaltum <- baflogr$t_ref + baflogr$t_alt

baflogr$logr <- log2( (totaltum/totalhost) / median((totaltum/totalhost), na.rm = T) )
rm(totalhost,totaltum)

baflogr$ntot <- ((2*(1-rho) + rho*psit)*2^baflogr$logr - 2*(1-rho))/rho

p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos, y = ntot)) + geom_point(shape = ".", alpha = .2) + facet_wrap(~chr) + ylim(c(-.2, 6))
p1


### continue with reconstruction of BAF and minor allele CN
### get "host" het and homo SNP loci

### host = hom if (coverage ≥ 10, as above) and 0 reads ref | alt
### host = het if ≥ 3 reads for both ALT and REF & .25 ≤ BAF ≤ .75
baflogr$htgeno <- ifelse(baflogr$h_ref == 0 & baflogr$t_ref >= 3, "BB/A*", ifelse(baflogr$h_alt  == 0 & baflogr$t_alt >= 3, "AA/B*",
                                                                                  ifelse(baflogr$h_ref >= 3 & baflogr$h_alt >= 3 & ( (baflogr$h_alt/(baflogr$h_alt+baflogr$h_ref)) >= .25 | (baflogr$h_alt/(baflogr$h_alt+baflogr$h_ref)) <= .75 ), "AB/*", NA)))

baflogr$baf <- baflogr$t_alt/(baflogr$t_alt + baflogr$t_ref)

baflogr$na <- ifelse(baflogr$htgeno == "AB/*", ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr - (1-rho))/rho,
                     ifelse(baflogr$htgeno == "AA/B*", ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr - 2*(1-rho))/rho, 
                            ifelse(baflogr$htgeno == "BB/A*", ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr)/rho, NA)))
baflogr$nb <- ifelse(baflogr$htgeno == "AB/*", ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr - (1-rho))/rho,
                     ifelse(baflogr$htgeno == "AA/B*", ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr)/rho, 
                            ifelse(baflogr$htgeno == "BB/A*", ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr - 2*(1-rho))/rho, NA)))

baflogr <- baflogr[!is.na(baflogr$htgeno), ]

baflogr$nmin <- ifelse(baflogr$na <= baflogr$nb, baflogr$na, baflogr$nb)

# source("/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/kataegis/20180309_dpclust3p_fastPCF.R")
# 
# debug(run_pcf_chr)
# cnsegs <- by(data = baflogr[baflogr$chr == "Chr2", ], INDICES = baflogr[baflogr$chr == "Chr2", "chr"], FUN = run_pcf_chr, kmin = 1001, gamma = 25)
# 
# # run pcf, omitting first mut in dinucs
# run_pcf_chr <- function(cndf, kmin, gamma = 25) {
#   goodcalls <- which(!is.na(cndf$nmin))
#   sdev <- getMad(cndf[goodcalls, "nmin"], k = 101)
#   res <- selectFastPcf(cndf[goodcalls, "nmin"], kmin, gamma*sdev, T)
#   return(res)
# }

# baflogr$rmed <- runmed(x = baflogr$nmin, k = 501, endrule = "constant")
baflogr$nminsm <- unlist(by(data = baflogr$nmin, INDICES = baflogr$chr, FUN = function(x) as.numeric(runmean(x = Rle(x), k = 1001, endrule = "constant"))))
baflogr$ntotsm <- unlist(by(data = baflogr$ntot, INDICES = baflogr$chr, FUN = function(x) runmed(x = x, k = 251, endrule = "constant")))


# checking
# testdf <- baflogr[sample(x = 1:nrow(baflogr), size = 100000, replace = F), ]
p1 <- ggplot(data = baflogr[baflogr$chr != "Chrx",], mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = ntot), colour = "goldenrod2")
p1 <- p1 + geom_point(shape = ".", alpha = .2, mapping = aes(y = nmin), colour = "darkslategrey") + geom_point(mapping = aes(y = nminsm), colour = "red", alpha = .5, shape = ".")
p1 <- p1 + geom_point(mapping = aes(y = ntotsm), colour = "black", alpha = .5, shape = ".")
p1 <- p1 + facet_wrap(~chr)
p1 <- p1+ylim(c(-.2, 6))
p1

# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180426_CN-transformed_BAFLogR.png", plot = p1, width = 16, height = 6)

# plot(1:length(cnsegs$Chr2$yhat), cnsegs$Chr2$yhat)

# p1 <- ggplot(data = baflogr, mapping = aes(x = ntotsm, fill = htgeno)) + geom_histogram(binwidth = .01)
p1 <- ggplot(data = baflogr, mapping = aes(x = ntot, fill = htgeno)) + geom_histogram(binwidth = .01)
p1 <- p1 + geom_vline(xintercept = .1)
p1 <- p1 + facet_wrap(~chr)
p1 <- p1 + xlim(c(-.2, 3))
p1


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
nminoneidxs <- which(baflogr$nmin >= (nminonepeak + nminzeropeak)/2)
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

p1 <- ggplot(data = baflogr[baflogr$type %in% c("AB/AB", "AB/A") & baflogr$chr != "Chrx", ], mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = baf), colour = "goldenrod2")
p1 <- p1 + facet_wrap(~chr)
# p1 <- p1+ylim(c(-.2, 6))
p1
# ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180426_input_BAF_first_ASCAT_run.png", plot = p1, width = 16, height = 6)


write.table(x = baflogr, file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180425_selectedSNPs_counts.txt", quote = F, sep = "\t")
