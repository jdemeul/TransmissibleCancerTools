## redo CN calls

devilinfo_file <- "../reference/Devil7.1.fa.fai"
dft2_dna_file <- "../results/DFT2_202T2_cell_line_DNA_alleleCount.txt"
fibro_dna_file <- "../results/Fibroblast_91H_cell_line_DNA_alleleCount.txt"

devilinfo <- read.delim(file = devilinfo_file, as.is = T, header = F)
colnames(devilinfo) <- c("contig", "width")
devilinfo$chr <- substr(x = devilinfo$contig, start = 1, stop = 4)
devilinfo$start <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = function(x) cumsum(c(1,x))[-length(x)]))
devilinfo$end <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = cumsum))
devilinfo_gr <- GRanges(seqnames = devilinfo$chr, ranges = IRanges(start = devilinfo$start, end = devilinfo$end), mcols = DataFrame(contig = devilinfo$contig))

dft2_dna <- read.delim(file = dft2_dna_file, as.is = T)
fibro_dna <- read.delim(file = fibro_dna_file, as.is = T)

## create "goodSNPs file"
### subset contigs?
contiglengths <- sort(width(devilinfo_gr), decreasing = T)
# plot((1:length(contiglengths))/length(contiglengths), cumsum(as.numeric(contiglengths)))
# plot(log10(contiglengths), cumsum(as.numeric(contiglengths)))

# allowed_contigs <- c(paste0("Chr1_supercontig_000000", sprintf("%03d", 0:398)),
#                      paste0("Chr2_supercontig_000000", sprintf("%03d", 0:498)),
#                      paste0("Chr3_supercontig_000000", sprintf("%03d", 0:416)),
#                      paste0("Chr4_supercontig_000000", sprintf("%03d", 0:316)),
#                      paste0("Chr5_supercontig_000000", sprintf("%03d", 0:217)),
#                      paste0("Chr6_supercontig_000000", sprintf("%03d", 0:193)))
allowed_contigs_length <- devilinfo[devilinfo$width >= 1e6, "contig"]
mincoverage <- 6

dft2_dna <- dft2_dna[dft2_dna$Good_depth >= mincoverage & dft2_dna$X.CHR %in% allowed_contigs_length, ]
fibro_dna <- fibro_dna[fibro_dna$Good_depth >= mincoverage & fibro_dna$X.CHR %in% allowed_contigs_length, ]

dft2_dna$chrs <- substr(dft2_dna$X.CHR, 4, 4)
# dft2_dna$chrs <- ifelse(dft2_dna$chrs == "x", 7, dft2_dna$chrs)
dft2_dna$pos <- start(devilinfo_gr[match(x = dft2_dna$X.CHR, table = mcols(devilinfo_gr)$mcols.contig)]) + dft2_dna$POS

fibro_dna$chrs <- substr(fibro_dna$X.CHR, 4, 4)
# fibro_dna$chrs <- ifelse(fibro_dna$chrs == "x", 7, fibro_dna$chrs)
fibro_dna$pos <- start(devilinfo_gr[match(x = fibro_dna$X.CHR, table = mcols(devilinfo_gr)$mcols.contig)]) + fibro_dna$POS

## create logR data
dft2_dna$DFT2 <- log2(dft2_dna$Good_depth/median(dft2_dna$Good_depth))
fibro_dna$Fibro <- log2(fibro_dna$Good_depth/median(fibro_dna$Good_depth))
rownames(dft2_dna) <- paste0("SNP", 1:nrow(dft2_dna))
rownames(fibro_dna) <- paste0("SNP", 1:nrow(fibro_dna))

write.table(x = dft2_dna[, c("chrs", "pos", "DFT2")], file = "../results/DFT2_202T2_cell_line_DNA_LogR.txt", sep = "\t", row.names = T, quote = F)
write.table(x = fibro_dna[, c("chrs", "pos", "Fibro")], file = "../results/Fibroblast_91H_cell_line_DNA_LogR.txt", sep = "\t", row.names = T, quote = F)

## create BAF data
get_baf <- function(counts) {
  baf <- max(counts)/sum(counts)
  baf <- sample(x = c(baf, 1 - baf), size = 1)
  return(baf)
}

dft2_dna$DFT2 <- apply(X = dft2_dna[, paste0("Count_", c("A", "C", "G", "T"))], MARGIN = 1, FUN = get_baf)
fibro_dna$Fibro <- apply(X = fibro_dna[, paste0("Count_", c("A", "C", "G", "T"))], MARGIN = 1, FUN = get_baf)


write.table(x = dft2_dna[, c("chrs", "pos", "DFT2")], file = "../results/DFT2_202T2_cell_line_DNA_BAF.txt", sep = "\t", row.names = T, quote = F)
write.table(x = fibro_dna[, c("chrs", "pos", "Fibro")], file = "../results/Fibroblast_91H_cell_line_DNA_BAF.txt", sep = "\t", row.names = T, quote = F)


# find homozygous stretches/contigs and exclude these from CN calls

homfracdft <- c(by(data = dft2_dna$DFT2, INDICES = dft2_dna$X.CHR, FUN = function(x) sum(x >= .95 | x <= 0.05)/length(x)))
homfracfibro <- c(by(data = fibro_dna$Fibro, INDICES = fibro_dna$X.CHR, FUN = function(x) sum(x >= .95 | x <= 0.05)/length(x)))
plot(1:length(homfracdft), homfracdft)
plot(1:length(homfracfibro), homfracfibro)
lohcontigs <- names(which(homfracdft >= .95 | homfracfibro >= .95))

### running ascat do it twice for Fibro and DFT2
library(ASCAT)
ascat.bc <- ascat.loadData("../results/Fibroblast_91H_cell_line_DNA_LogR.txt","../results/Fibroblast_91H_cell_line_DNA_BAF.txt")
setwd("../results/")
ascat.plotRawData(ascat.bc)

source("../code/20180126_dft_modifiedASCAT.R")

ascat.gg <- ascat.predictGermlineGenotypes.mod(ascat.bc, "DevilNGS")
ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg=ascat.gg)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma = 1, rho_manual = 1, psi_manual = 2.18)
# View(ascat.output$segments)
write.table(x = ascat.output$segments, file = "../results/Fibro.ASCATprofile.segments.txt", sep = "\t", quote = F, row.names = F)



cnsegments <- GRanges(seqnames = paste0("Chr", ascat.output$segments$chr), ranges = IRanges(start = ascat.output$segments$startpos, end = ascat.output$segments$endpos), mcols = ascat.output$segments[, c("nMajor", "nMinor")])
lohcontigs_gr <- reduce(devilinfo_gr[devilinfo$contig %in% lohcontigs])
cnsegments_disjoint <- disjoin(c(granges(cnsegments), granges(lohcontigs_gr)))
cnsegments_disjoint <- cnsegments_disjoint[!overlapsAny(query = cnsegments_disjoint, subject = lohcontigs_gr) & width(cnsegments_disjoint) >= 1e6]
cnsegment_overlaps <- findOverlaps(query = cnsegments_disjoint, subject = cnsegments)
mcols(cnsegments_disjoint) <- mcols(cnsegments[subjectHits(cnsegment_overlaps)])
write.table(x = as.data.frame(cnsegments_disjoint)[, c(1:3, 6:7)], file = "../results/Fibro.ASCATprofile.segments_clean_noLOH.txt", sep = "\t", quote = F, row.names = F)
