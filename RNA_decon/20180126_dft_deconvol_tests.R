## reading allelecounts
# check 100% DFT and 100% fibroblast first
library(ggplot2)
library(GenomicRanges)
library(nnls)

devilinfo_file <- "../reference/Devil7.1.fa.fai"

dft2_dna_file <- "../results/DFT2_202T2_cell_line_DNA_alleleCount.txt"
fibro_dna_file <- "../results/Fibroblast_91H_cell_line_DNA_alleleCount.txt"
# dft2_cn_file <- "../results/DFT2_CN_profiles.txt"
dft2_cn_file <- "../results/DFT2.ASCATprofile.segments_clean_noLOH.txt"
fibro_cn_file <- "../results/Fibro.ASCATprofile.segments_clean_noLOH.txt"

devilinfo <- read.delim(file = devilinfo_file, as.is = T, header = F)
colnames(devilinfo) <- c("contig", "width")
devilinfo$chr <- substr(x = devilinfo$contig, start = 1, stop = 4)
devilinfo$start <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = function(x) cumsum(c(1,x))[-length(x)]))
devilinfo$end <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = cumsum))
devilinfo_gr <- GRanges(seqnames = devilinfo$chr, ranges = IRanges(start = devilinfo$start, end = devilinfo$end), mcols = DataFrame(contig = devilinfo$contig))

dft2_dna <- read.delim(file = dft2_dna_file, as.is = T)
fibro_dna <- read.delim(file = fibro_dna_file, as.is = T)

## adding the "position along the chromosome" to one of the dataframes
dft2_dna$chr <- substr(dft2_dna$X.CHR, 1, 4)
dft2_dna$cumpos <- start(devilinfo_gr[match(x = dft2_dna$X.CHR, table = mcols(devilinfo_gr)$mcols.contig)]) + dft2_dna$POS

fibro_dna$chr <- substr(fibro_dna$X.CHR, 1, 4)
fibro_dna$cumpos <- start(devilinfo_gr[match(x = fibro_dna$X.CHR, table = mcols(devilinfo_gr)$mcols.contig)]) + fibro_dna$POS

# load and "fix" CN
# dft2_cn <- read.delim(file = dft2_cn_file, as.is = T)
# dft2_cn$chr <- substr(dft2_cn$Scaffold.Start, 1, 4)
# 
# dft2_cn$end <- unlist(by(data = as.numeric(dft2_cn$Width), INDICES = dft2_cn$chr, FUN = cumsum))
# dft2_cn$start <- unlist(by(data = as.numeric(dft2_cn$Width), INDICES = dft2_cn$chr, FUN = function(x) cumsum(c(1,x))[-length(x)]))
# 
# dict <- c(CN0 = 0, CN1 = 1, CN2 = 2, GAIN = 3)
# dft2_cn$CN <- dict[dft2_cn$X202T2]
# dft2_cn_gr <- GRanges(seqnames = dft2_cn$chr, ranges = IRanges(start = dft2_cn$start, end = dft2_cn$end), mcols = DataFrame(cn = dft2_cn$CN))

### new bit, loading ASCAT generated + LOH removed CN
dft2_cn <- read.delim(file = dft2_cn_file, as.is = T)
dft2_cn_gr <- GRanges(seqnames = dft2_cn$seqnames, ranges = IRanges(start = dft2_cn$start, end = dft2_cn$end), mcols = dft2_cn[ , 4:5])
colnames(mcols(dft2_cn_gr)) <- c("nMajor", "nMinor")
fibro_cn <- read.delim(file = fibro_cn_file, as.is = T)
fibro_cn_gr <- GRanges(seqnames = fibro_cn$seqnames, ranges = IRanges(start = fibro_cn$start, end = fibro_cn$end), mcols = fibro_cn[ , 4:5])
colnames(mcols(fibro_cn_gr)) <- c("nMajor", "nMinor")

# check coverage
# p1 <- ggplot(data = dft2_dna[dft2_dna$Good_depth > 0, ], mapping = aes(x = Good_depth)) + geom_histogram() + scale_y_log10() + scale_x_log10() 
# p1


# make genotype calls in "tumour" and "host"
# for tumour
mintotdepth <- 6
dft2_dna_cov <- dft2_dna[rowSums(dft2_dna[, paste0("Count_", c("A", "C", "G", "T"))]) >= mintotdepth, ]

# bafcheck <- data.frame(chr = dft2_dna_cov$chr, pos = dft2_dna_cov$chrpos, baf = apply(X = dft2_dna_cov[, paste0("Count_", c("A", "C", "G", "T"))], MARGIN = 1, FUN = max)/dft2_dna_cov$Good_Depth)
# p1 <- ggplot(bafcheck[bafcheck$chr == "Chr2", ], mapping = aes(x = pos, y = baf)) + geom_point(shape = ".", alpha = .5)
# p1
# 
# bafcheck_fibro <- data.frame(chr = fibro_dna$chr, pos = fibro_dna$cumpos, baf = apply(X = fibro_dna[, paste0("Count_", c("A", "C", "G", "T"))], MARGIN = 1, FUN = max)/fibro_dna$Good_depth)
# p1 <- ggplot(bafcheck_fibro[bafcheck_fibro$chr == "Chr1", ], mapping = aes(x = pos, y = baf)) + geom_point(shape = ".", alpha = .5)
# p1


get_genotype <- function(allelecount, bases = c("A", "C", "G", "T"), mindepth = 3, minbaf = .1) {
  counts_sorted <- sort(setNames(object = allelecount, nm = bases), decreasing = T)
  if (counts_sorted[2] >= mindepth && counts_sorted[2]/sum(counts_sorted) >= minbaf) {
    genotype <- names(counts_sorted)[1:2]
  } else {
    genotype <- rep(names(counts_sorted)[1], 2)
  }
  return(genotype)
}

dft2_geno <- as.data.frame(t(apply(X = dft2_dna_cov[, paste0("Count_", c("A", "C", "G", "T"))], MARGIN = 1, FUN = get_genotype)))
dft2_dna_cov <- cbind(dft2_dna_cov, dft2_geno)
colnames(dft2_dna_cov) <- c("contig", "contigpos", paste0("Count_", c("A", "C", "G", "T")), "Good Depth", "chr", "chrpos", "aMaj", "aMin")

write.table(x = dft2_dna_cov[, c("contig", "contigpos", "aMaj", "aMin")], file = "../results/DFT2_202T2_cell_line_DNA_genotype.txt", sep = "\t", quote = F, row.names = F)


# for "normal"
fibro_dna_cov <- fibro_dna[rowSums(fibro_dna[, paste0("Count_", c("A", "C", "G", "T"))]) >= mintotdepth, ]

fibro_geno <- as.data.frame(t(apply(X = fibro_dna_cov[, paste0("Count_", c("A", "C", "G", "T"))], MARGIN = 1, FUN = get_genotype)))
fibro_dna_cov <- cbind(fibro_dna_cov, fibro_geno)
colnames(fibro_dna_cov) <- c("contig", "contigpos", paste0("Count_", c("A", "C", "G", "T")), "Good Depth", "chr", "chrpos", "aMaj", "aMin")

write.table(x = fibro_dna_cov[, c("contig", "contigpos", "aMaj", "aMin")], file = "../results/Fibroblast_91H_cell_line_DNA_genotype.txt", sep = "\t", quote = F, row.names = F)



####### combine genotype calls and subset informative sites
informative_loci <- merge(x = dft2_dna_cov[, c("contig", "contigpos", "chr", "chrpos", "aMaj", "aMin")], y = fibro_dna_cov[, c("contig", "contigpos", "aMaj", "aMin")], by = c("contig", "contigpos"))
informative_loci <- informative_loci[!(informative_loci$aMaj.x == informative_loci$aMin.x & informative_loci$aMaj.x == informative_loci$aMaj.y & informative_loci$aMaj.x == informative_loci$aMin.y), ]

classify_locus <- function(genotypes) {
  if (length(unique(genotypes)) == 1) {
    return(c(geno = "AA/AA", A = genotypes[1], B = genotypes[1]))
  } else if (length(unique(genotypes)) >= 3) {
    return(c(geno = NA, A = NA, B = NA))
  } else if (setequal(genotypes[1:2], genotypes[3:4])) {
    return(c(geno = "AB/AB", A = genotypes[1], B = genotypes[2]))
  } else if (genotypes[1] == genotypes[2]) {
    if (genotypes[3] == genotypes[4]) {
      return(c(geno = "AA/BB", A = genotypes[1], B = genotypes[3]))
    } else {
      return(c(geno = "AA/AB", A = genotypes[1], B = setdiff(genotypes[3:4], genotypes[1:2])))
    }
  } else if (genotypes[3] == genotypes[4]) {
    if (genotypes[1] == genotypes[3]) {
      return(c(geno = "AB/AA", A = genotypes[1], B = genotypes[2]))
    } else {
      return(c(geno = "AB/BB", A = genotypes[1], B = genotypes[2]))
    }
  } else {
    return(c(geno = NA, A = NA, B = NA))
  }
}


informative_loci <- cbind(informative_loci, as.data.frame(t(apply(X = informative_loci[, c("aMaj.x", "aMin.x", "aMaj.y", "aMin.y")], MARGIN = 1, FUN = classify_locus))))
colnames(informative_loci)[9:11] <- c("geno", "A", "B")


## create contig-disjoined CN segments (basically have every contig or CN segment annotated with CN)
matchsegs <- disjoin(x = c(granges(dft2_cn_gr), granges(fibro_cn_gr), granges(devilinfo_gr)))
matchseg_backmap_dft2 <- findOverlaps(query = matchsegs, subject = dft2_cn_gr)
matchseg_backmap_fibro <- findOverlaps(query = matchsegs, subject = fibro_cn_gr)
matchsegs_cndict <- matchsegs
mcols(matchsegs_cndict)[, c("nMajT", "nMinT", "nMajH", "nMinH")] <- DataFrame(nMajT = integer(length = length(matchsegs_cndict)),
                                                                              nMinT = integer(length = length(matchsegs_cndict)),
                                                                              nMajH = integer(length = length(matchsegs_cndict)),
                                                                              nMinH = integer(length = length(matchsegs_cndict)))
mcols(matchsegs_cndict[queryHits(matchseg_backmap_dft2)])[, c("nMajT", "nMinT")] <- mcols(dft2_cn_gr[subjectHits(matchseg_backmap_dft2)])
mcols(matchsegs_cndict[queryHits(matchseg_backmap_fibro)])[, c("nMajH", "nMinH")] <- mcols(fibro_cn_gr[subjectHits(matchseg_backmap_fibro)])


## CN matching for SNP loci
snphits <- findOverlaps(query = GRanges(seqnames = informative_loci$chr, IRanges(start = informative_loci$chrpos, end = informative_loci$chrpos)),
                        subject = matchsegs_cndict)

informative_loci <- informative_loci[queryHits(snphits), ]
informative_loci[ , c("nMajT", "nMinT", "nMajH", "nMinH")] <- as.data.frame(mcols(matchsegs_cndict[subjectHits(snphits)]))
informative_loci <- informative_loci[rowSums(informative_loci[, c("nMajT", "nMinT", "nMajH", "nMinH")]) > 0, ]

## in this case, with ASCAT running on tumour only, do not trust LOH calls
informative_loci <- informative_loci[informative_loci$nMinH > 0 & informative_loci$nMinT > 0, ]

# for "host"
# informative_loci$nHtot <- 2
# informative_loci[informative_loci$chr == "Chr6", "nHtot"] <- 3

## below is replaced by ASCAT CN
# ## create allele-specific CN
# informative_loci[ , c("nMajT", "nMinT", "nMajH", "nMinH")] <- NA
# informative_loci$nMajT <- ifelse(informative_loci$nTtot == 3, 2, ifelse(informative_loci$nTtot == 2, 1, ifelse(informative_loci$nTtot == 1, 1, 0)))
# informative_loci$nMinT <- ifelse(informative_loci$nTtot == 3, 1, ifelse(informative_loci$nTtot == 2, 1, 0))
# informative_loci$nMajH <- ifelse(informative_loci$nHtot == 3, 2, ifelse(informative_loci$nHtot == 2, 1, ifelse(informative_loci$nTtot == 1, 1, 0)))
# informative_loci$nMinH <- ifelse(informative_loci$nHtot == 3, 1, ifelse(informative_loci$nHtot == 2, 1, 0))
# 
# ## fix cnLOH on Chr2 in T and chr1 in H
# informative_loci$nMajT <- ifelse(informative_loci$chr == "Chr2" & ((informative_loci$chrpos > 0.6E8 & informative_loci$chrpos < 1.6E8) | (informative_loci$chrpos > 1.8E8 & informative_loci$chrpos < 2.35E8)), 2, informative_loci$nMajT )
# informative_loci$nMinT <- ifelse(informative_loci$chr == "Chr2" & ((informative_loci$chrpos > 0.6E8 & informative_loci$chrpos < 1.6E8) | (informative_loci$chrpos > 1.8E8 & informative_loci$chrpos < 2.35E8)), 0, informative_loci$nMinT )
# 
# informative_loci$nMajH <- ifelse(informative_loci$chr == "Chr1" & (informative_loci$chrpos > 3.9E8 & informative_loci$chrpos < 5.1E8), 2, informative_loci$nMajH )
# informative_loci$nMajH <- ifelse(informative_loci$chr == "Chr1" & (informative_loci$chrpos > 3.9E8 & informative_loci$chrpos < 5.1E8), 0, informative_loci$nMinH )

## create allelic counts homo-heterozygous
informative_loci$nAT <- ifelse(grepl(pattern = "AA/", x = informative_loci$geno), informative_loci$nMajT + informative_loci$nMinT, informative_loci$nMajT)
informative_loci$nBT <- ifelse(grepl(pattern = "AA/", x = informative_loci$geno), 0, informative_loci$nMinT)
informative_loci$nAH <- ifelse(grepl(pattern = "/AA", x = informative_loci$geno), informative_loci$nMajH + informative_loci$nMinH,
                               ifelse(grepl(pattern = "/AB", x = informative_loci$geno), informative_loci$nMajH, 0))
informative_loci$nBH <- ifelse(grepl(pattern = "/AA", x = informative_loci$geno), 0,
                               ifelse(grepl(pattern = "/AB", x = informative_loci$geno), informative_loci$nMinH, informative_loci$nMajH + informative_loci$nMinH))

informative_loci <- informative_loci[informative_loci$nAT != informative_loci$nBT | informative_loci$nAH != informative_loci$nBH, ]

## write informative loci file
write.table(x = informative_loci, file = "../results/DFT2_202T2_Fibro_91H_informative_loci.txt", sep = "\t", quote = F, row.names = F)





###### start RNA-Seq data
## load mixture RNA-Seq
dft2_50pc_file <- "../results/DFT2_50perc_RNA_alleleCount.txt"
dft2_50pc <- read.delim(file = dft2_50pc_file, as.is = T)
dft2_50pc <- dft2_50pc[dft2_50pc$Good_depth >= 6, ]

dft2_50pc <- merge(x = dft2_50pc, y = informative_loci[, c("contig", "contigpos", "chr", "chrpos", "A", "B", "nAT", "nBT", "nAH", "nBH")], by.x = c("X.CHR", "POS"), by.y = c("contig", "contigpos"), sort = F)
dft2_50pc$eA <- ifelse(dft2_50pc$A == "A", dft2_50pc$Count_A,
                       ifelse(dft2_50pc$A == "C", dft2_50pc$Count_C, 
                              ifelse(dft2_50pc$A == "G", dft2_50pc$Count_G, dft2_50pc$Count_T)))
dft2_50pc$eB <- ifelse(dft2_50pc$B == "A", dft2_50pc$Count_A,
                        ifelse(dft2_50pc$B == "C", dft2_50pc$Count_C, 
                               ifelse(dft2_50pc$B == "G", dft2_50pc$Count_G, dft2_50pc$Count_T)))
dft2_50pc <- dft2_50pc[!(is.na(dft2_50pc$eA) | is.na(dft2_50pc$eB)), ]

head(dft2_50pc)

get_et_en <- function(sysdata, rho = .5) {
  A <- matrix(data = c((1-rho)*sysdata[["nAH"]], rho*sysdata[["nAT"]], (1-rho)*sysdata[["nBH"]], rho*sysdata[["nBT"]]), nrow = 2, byrow = T)
  b <- sysdata[c("eA", "eB")]
  x <- nnls(A = A, b = b)
  # x <- nnls(A = A, b = b)$x
  return(x$x)
  # return(c(eN = x$x[1], eT = x$x[2], res = x$residuals[1], dev = x$deviance))
}

dft2_50pc[, c("eN", "eT", "res", "dev")] <- t(apply(X = dft2_50pc[, c("nAT", "nBT", "nAH", "nBH", "eA", "eB")], MARGIN = 1, FUN = get_et_en))

p1 <- ggplot(data = dft2_50pc, mapping = aes(x = eN + 1, y = eT +1, colour = interaction(nAT, nBT))) + geom_point(alpha = .5) + scale_x_log10(limits = c(1E-1, 1E3)) + scale_y_log10(limits = c(1E-1, 1E3))
p1

#### out of pipeline: add in results for 25 and 75% DFT2
dft2_50pc_file <- "../results/DFT2_25perc_RNA_alleleCount.txt"
dft2_50pc_file <- "../results/DFT2_75perc_RNA_alleleCount.txt"

dft2_50pc[, c("eN25", "eT25")] <- t(apply(X = dft2_50pc[, c("nAT", "nBT", "nAH", "nBH", "eA", "eB")], MARGIN = 1, FUN = get_et_en, rho = .25))
comparisondf <- merge(x = comparisondf, y = dft2_50pc[, c("X.CHR", "POS", "eN25", "eT25")], by.x = c("contig", "contigpos"), by.y = c("X.CHR", "POS"), all.x = T)

dft2_50pc[, c("eN75", "eT75")] <- t(apply(X = dft2_50pc[, c("nAT", "nBT", "nAH", "nBH", "eA", "eB")], MARGIN = 1, FUN = get_et_en, rho = .75))
comparisondf <- merge(x = comparisondf, y = dft2_50pc[, c("X.CHR", "POS", "eN75", "eT75")], by.x = c("contig", "contigpos"), by.y = c("X.CHR", "POS"), all.x = T)
#### end of this out of pipeline test

### cross-check with 100% DFT2 / 100% Fibro expression
dft2_100pc_file <- "../results/DFT2_100perc_RNA_alleleCount.txt"
dft2_0pc_file <- "../results/DFT2_0perc_RNA_alleleCount.txt"

dft2_100pc <- read.delim(file = dft2_100pc_file, as.is = T)
dft2_0pc <- read.delim(file = dft2_0pc_file, as.is = T)

comparisondf <- merge(x = dft2_50pc, y = dft2_100pc[, c("X.CHR", "POS", "Good_depth")], by = c("X.CHR", "POS"))
comparisondf <- merge(x = comparisondf, y = dft2_0pc[, c("X.CHR", "POS", "Good_depth")], by = c("X.CHR", "POS"))


colnames(comparisondf) <- c("contig", "contigpos", paste0("Count_", c("A", "C", "G", "T")), "Good_depth_50", "chr", "chrpos", "A", "B", "nAT", "nBT", "nAH", "nBH", "eA", "eB", "eN", "eT", "res", "dev", "eDFT2", "eFibro")
write.table(x = comparisondf, file = "../results/DFT2_202T2_Fibro_91H_fit_eT-eN_vs_pure_25and75.txt", sep = "\t", quote = F, row.names = F)


View(comparisondf[abs(comparisondf$res) >= .5, ])

p1 <- ggplot(data = comparisondf[comparisondf$Good_depth_50 < 1000,], mapping = aes(x = log10(Good_depth_50), y = res)) + geom_point(mapping = aes(colour = chr), alpha = .5)
p1


p1 <- ggplot(data = comparisondf, mapping = aes(x = eN+1, y = (eFibro+1)/(nAH+nBH), colour = interaction(nAT, nBT))) + geom_point(alpha = .5) + scale_x_log10() + scale_y_log10()
p1
ggsave(filename = "../results/20180401_eN_vs_Fibro.png", plot = p1)

p1 <- ggplot(data = comparisondf, mapping = aes(x = eT+1, y = (eDFT2+1)/(nAT+nBT), colour = interaction(nAT, nBT))) + geom_point(alpha = .5) + scale_x_log10() + scale_y_log10()
p1
ggsave(filename = "../results/20180401_eT_vs_DFT2.png", plot = p1)

p1 <- ggplot(data = comparisondf, mapping = aes(x = log2((eT+1)/(eN+1)), y = log2(((eDFT2+1)/(nAT+nBT))/((eFibro+1)/(nAH+nBH))))) + geom_abline(slope = 1) + geom_point(alpha = .5, mapping = aes(colour = interaction(nAT, nBT)))
p1 <- p1 + geom_density2d(alpha = .5)
p1
ggsave(filename = "../results/20180401_eT-eN_vs_DFT2-Fibro.png", plot = p1)

p1 <- ggplot(data = comparisondf, mapping = aes(x = log2((eFibro + 1)/(nAH+nBH)), y = log2((eDFT2+1)/(nAT+nBT)))) + geom_point(alpha = .5, mapping = aes(colour = interaction(nAT, nBT)))
p1
ggsave(filename = "../results/20180401_DFT2_vs_Fibro_sanitycheck.png", plot = p1)

p1 <- ggplot(data = comparisondf, mapping = aes(x = chrpos, y = log2((eT+1)/((eDFT2+1)/(nAT+nBT))) )) + geom_point(alpha = .5) + facet_wrap(~chr)
p1


p1 <- ggplot(data = comparisondf[comparisondf$eT < 1 & comparisondf$eDFT2 >= 6, ], mapping = aes(x = log2(eN+1), y = log2((eFibro+1)/(nAT+nBT)), colour = interaction(nAH, nBH))) + geom_point()
p1

p1 <- ggplot(data = comparisondf[, ], mapping = aes(x = chrpos, y = log2((eT+1)/((eDFT2+1)/(nAT+nBT))), colour = interaction(nAT, nBT, nAH, nBH))) + geom_point(alpha = .5) + facet_wrap(~chr)
p1
ggsave(filename = "../results/20180401_error_on_eT_across_genome.png", plot = p1)

p1 <- ggplot(data = comparisondf, mapping = aes(x = chrpos)) + geom_point(mapping = aes(y = nAT), colour = "red", alpha = .5) + geom_point(mapping = aes(y = nBT), colour = "blue", alpha = .5) + facet_wrap(~chr)
p1

p1 <- ggplot(data = comparisondf, mapping = aes(x = chrpos)) + geom_point(mapping = aes(y = log10(eFibro)), colour = "red", alpha = .5) + facet_wrap(~chr)
p1

View(comparisondf[comparisondf$eFibro >= 10^3.5, ])

p1 <- ggplot(data = comparisondf, mapping = aes(x = log2((eT+1)/((eDFT2+1)/(nAT+nBT))) )) + geom_histogram(bins = 100) + lims(x = c(-10, 10))
p1
ggsave(filename = "../results/20180401_error_on_eT_histo.png", plot = p1)

p1 <- ggplot(data = fibro_dna[grepl(pattern = "Chr4_supercontig_00000018", x = fibro_dna[fibro_dna$X.CHR]), ] , mapping = aes(x = cumpos)) + geom_point(mapping = aes(y = Good_depth), alpha = .5) + geom_point(mapping = aes(y = nBT), colour = "blue", alpha = .5) + facet_wrap(~chr)
p1


p1 <- ggplot(data = comparisondf, mapping = aes(x = log10(Good_depth_50+1), y = log2((eT+1)/((eDFT2+1)/(nAT+nBT)))) ) + geom_point(mapping = aes(colour = interaction(nAT, nBT)))
p1
ggsave(filename = "../results/20180401_error_on_eT_vs_bulkcoverage.png", plot = p1)

