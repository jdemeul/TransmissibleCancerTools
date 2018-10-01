## read baf/logr data
# finalbaflogr <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180425_selectedSNPs_counts.txt", as.is = T)
# finalbaflogr <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180508_selectedSNPs_counts.txt", as.is = T)
allbaflogrfile <- "Fibro_DFT2_50-50MIX_genotyped_BAF-Logr_all.txt"
finalbaflogr <- as.data.frame(read_tsv(file = allbaflogrfile))

# # str(ascat.output)
# psi <- ascat.output$psi
# rho <- ascat.output$aberrantcellfraction



### ASCAT run 2
setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/Fibro_DFT2_50-50MIX/ascatrun2")
### use ASCAT-called hom segments to redo the AB/A selection
baflogr_run2 <- finalbaflogr
# ascat_lohsegs <- ascat.output$segments[ascat.output$segments$nMinor == 0, ]
# ascat_lohsegs_gr <- GRanges(seqnames = ascat_lohsegs$chr, ranges = IRanges(start = ascat_lohsegs$startpos, end = ascat_lohsegs$endpos))
# ascat_lohsegs_gr <- ascat_lohsegs_gr[width(ascat_lohsegs_gr) >= 1e5]
# 
# baflogr_run2_gr <- GRanges(seqnames = baflogr_run2$chr, ranges = IRanges(start = baflogr_run2$cumpos, end = baflogr_run2$cumpos))
# 
# p1 <- ggplot(data = baflogr_run2, mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = ntot), colour = "goldenrod2")
# p1 <- p1 + geom_point(shape = ".", alpha = .2, mapping = aes(y = nmin), colour = "darkslategrey")
# p1 <- p1 + geom_segment(data = ascat_lohsegs, mapping = aes(x = startpos, xend = endpos, y = nMajor, yend = nMajor), colour = "red", size = 5, alpha = .75)
# p1 <- p1 + ylim(c(-.2, 6)) + theme_minimal() + facet_wrap(~chr)
# p1
# ggsave(filename = paste0(sampleid, "_identifiedLOHsegments_ascatrun1.png"), plot = p1)
# 
# lohidxs <- which(baflogr_run2_gr %over% ascat_lohsegs_gr)
# hetidxs <- which(baflogr_run2$htgeno == "AB/*")
# 
# # reset typing for LOH regions
# baflogr_run2[baflogr_run2$type %in% c('AB/A'), "type"] <- NA
# 
# baflogr_run2[intersect(hetidxs, lohidxs), "type"] <- 'AB/A'


str(ascat.output)
ascat.output$psi
ascat.output$aberrantcellfraction

setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/")
