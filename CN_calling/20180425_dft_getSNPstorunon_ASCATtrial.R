## read baf/logr data
# finalbaflogr <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180425_selectedSNPs_counts.txt", as.is = T)
# finalbaflogr <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180508_selectedSNPs_counts.txt", as.is = T)
finalbaflogr <- as.data.frame(read_tsv(file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180508_selectedSNPs_counts.txt", col_names = T, col_types = "cicciiiinncnnnnnnc"))


# finalbaflogr$chr <- toupper(sub(pattern = "Chr", replacement = "", x = finalbaflogr$chr))

## create BAF/logR data ASCAT run 1
baflogr_run1 <- finalbaflogr[finalbaflogr$type %in% c("AB/AB", "AB/A"), ]
setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/ascat_snpsubset/")
write.table(x = data.frame(chr = baflogr_run1$chr, pos = baflogr_run1$cumpos, '559T' = baflogr_run1$logr), file = "20180508_559T_tumour_LogR.txt", sep = "\t", row.names = T, quote = F)
write.table(x = data.frame(chr = baflogr_run1$chr, pos = baflogr_run1$cumpos, '559N' = 0), file = "20180508_559T_host_LogR.txt", sep = "\t", row.names = T, quote = F)

write.table(x = data.frame(chr = baflogr_run1$chr, pos = baflogr_run1$cumpos, '559T' = abs(baflogr_run1$baf - sample(x = c(0,1), size = nrow(baflogr_run1), replace = T))), file = "20180508_559T_tumour_BAF.txt", sep = "\t", row.names = T, quote = F)
write.table(x = data.frame(chr = baflogr_run1$chr, pos = baflogr_run1$cumpos, '559N' = abs(baflogr_run1$h_alt/(baflogr_run1$h_ref+baflogr_run1$h_alt) - sample(x = c(0,1), size = nrow(baflogr_run1), replace = T))), file = "20180508_559T_host_BAF.txt", sep = "\t", row.names = T, quote = F)


### running ascat do it twice for Fibro and DFT2
library(ASCAT)
ascat.bc <- ascat.loadData("20180508_559T_tumour_LogR.txt", "20180508_559T_tumour_BAF.txt","20180508_559T_host_LogR.txt","20180508_559T_host_BAF.txt")
ascat.plotRawData(ascat.bc)

ascat.bc <- ascat.aspcf(ascat.bc, penalty = 250)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma = 1)
# View(ascat.output$segments)
write.table(x = ascat.output$segments, file = "559T.ASCATprofile.segments.txt", sep = "\t", quote = F, row.names = F)

# str(ascat.output)
psi <- ascat.output$psi
rho <- ascat.output$aberrantcellfraction


### ASCAT run 2
setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/ascat_snpsubset2/")
### use ASCAT-called hom segments to redo the AB/A selection
baflogr_run2 <- finalbaflogr
ascat_lohsegs <- ascat.output$segments[ascat.output$segments$nMinor == 0, ]
ascat_lohsegs_gr <- GRanges(seqnames = ascat_lohsegs$chr, ranges = IRanges(start = ascat_lohsegs$startpos, end = ascat_lohsegs$endpos))

baflogr_run2_gr <- GRanges(seqnames = baflogr_run2$chr, ranges = IRanges(start = baflogr_run2$cumpos, end = baflogr_run2$cumpos))
lohidxs <- which(baflogr_run2_gr %over% ascat_lohsegs_gr)
hetidxs <- which(baflogr_run2$htgeno == "AB/*")

# reset typing for LOH regions
baflogr_run2[baflogr_run2$type %in% c('AB/A'), "type"] <- NA
baflogr_run2[intersect(hetidxs, lohidxs), "type"] <- 'AB/A'


### create BAF/logR data ASCAT run 2

baflogr_run2 <- baflogr_run2[!is.na(finalbaflogr$type), ]
# baflogr_run2$baf <- ifelse(baflogr_run2$type %in% c("AA/AB", "BB/AB"), (baflogr_run2$nb*rho + (1-rho))/(2*(1-rho) + rho*baflogr_run2$ntot), baflogr_run2$baf)
baflogr_run2$baf <- ifelse(baflogr_run2$type %in% c("AA/AB"), (psi*2^baflogr_run2$logr*baflogr_run2$baf+(1-rho))/(psi*2^baflogr_run2$logr),
                           ifelse(baflogr_run2$type %in% c("BB/AB"), (psi*2^baflogr_run2$logr*baflogr_run2$baf-1*(1-rho))/(psi*2^baflogr_run2$logr),
                                  baflogr_run2$baf))
baflogr_run2$normalbaf <- ifelse(baflogr_run2$type %in% c("AA/AB", "BB/AB"), .5, abs(baflogr_run2$h_alt/(baflogr_run2$h_ref+baflogr_run2$h_alt) - sample(x = c(0,1), size = nrow(baflogr_run2), replace = T)))

write.table(x = data.frame(chr = baflogr_run2$chr, pos = baflogr_run2$cumpos, '559T' = baflogr_run2$logr), file = "20180508_559T_tumour_LogR.txt", sep = "\t", row.names = T, quote = F)
write.table(x = data.frame(chr = baflogr_run2$chr, pos = baflogr_run2$cumpos, '559N' = 0), file = "20180508_559T_host_LogR.txt", sep = "\t", row.names = T, quote = F)

write.table(x = data.frame(chr = baflogr_run2$chr, pos = baflogr_run2$cumpos, '559T' = abs(baflogr_run2$baf - sample(x = c(0,1), size = nrow(baflogr_run2), replace = T))), file = "20180508_559T_tumour_BAF.txt", sep = "\t", row.names = T, quote = F)
write.table(x = data.frame(chr = baflogr_run2$chr, pos = baflogr_run2$cumpos, '559N' = abs(baflogr_run2$normalbaf - sample(x = c(0,1), size = nrow(baflogr_run2), replace = T))), file = "20180508_559T_host_BAF.txt", sep = "\t", row.names = T, quote = F)


### running ascat
ascat.bc <- ascat.loadData("20180508_559T_tumour_LogR.txt", "20180508_559T_tumour_BAF.txt","20180508_559T_host_LogR.txt","20180508_559T_host_BAF.txt")
ascat.plotRawData(ascat.bc)

ascat.bc <- ascat.aspcf(ascat.bc, penalty = 250)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma = 1)
# View(ascat.output$segments)
write.table(x = ascat.output$segments, file = "559T.ASCATprofile.segments.txt", sep = "\t", quote = F, row.names = F)

str(ascat.output)
ascat.output$psi
ascat.output$aberrantcellfraction

setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/")
