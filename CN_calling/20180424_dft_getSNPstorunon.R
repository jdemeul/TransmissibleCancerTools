## reading allelecounts
# check 100% DFT and 100% fibroblast first
library(ggplot2)
library(GenomicRanges)
library(readr)
library(ASCAT)

devilinfo_file <- "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/reference/Devil7.1.fa.fai"

# tumour_allelecounts_file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/Fibroblast_91H_DFT2_202T2_50-50MIX_DNA_alleleCount.txt"
# host_allelecounts_file = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/Fibroblast_91H_cell_line_DNA_alleleCount.txt"

# alternative input KG format
datafile <- "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/data/DFT/devil_dft1_data.tsv.gz"
sampleid <- "818"
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
# baflogr <- get_baf_logr(tumour_allelecounts_file = tumour_allelecounts_file,
#                         host_allelecounts_file = host_allelecounts_file,
#                         reference_alleles_file = snploci_file,
#                         genomeinfo = devilinfo_gr,
#                         min_contig_size = 1e6)

# debug(get_baf_logr_kg_format)
baflogr <- get_baf_logr_kg_format(datafile = datafile, sampleid = sampleid,
                        genomeinfo = devilinfo_gr,
                        min_contig_size = 1e6)


## adding the "position along the chromosome"
# snploci$REF <- factor(snploci$REF, levels = c("A", "C", "G", "T"))
# snploci$ALT <- factor(snploci$ALT, levels = c("A", "C", "G", "T"))

### quick check of total CN reconstruction
# debug(check_rho_psi_estimates)
check_rho_psi_estimates(sampleid = sampleid, baflogr = baflogr, rho = .75, psit = 2)

# debug(genotype_loci)
baflogr <- genotype_loci(baflogr = baflogr,
                         rho = .75, psit = 2,
                         sampleid = sampleid,
                         plotting = T)


### at this point, please carefully look at the histogram and smoothed nmin/ntot to select appropriate thresholds
### probably a manual check is safer here than automating everything at the moment

# findpeaks <- function (x, thresh = 0.00001, from = -.75, to = 1.75, mindens = .25) {
#   xdens <- stats::density(x, kernel = "gaussian", from = from, to = to)
#   pks <- which(diff(sign(diff(xdens$y, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
#   pks <- pks[xdens$y[pks - 1] - xdens$y[pks] > thresh & xdens$y[pks] > mindens]
#   return(xdens$x[pks])
# }
# 
# nminpeaks <- findpeaks(x = baflogr$nmin, thresh = .00001)
# nminonepeak <- nminpeaks[which.min(abs(nminpeaks - 1))]
# nminzeropeak <- nminpeaks[which.min(abs(nminpeaks - 0))]
## ntotpeaks <- findpeaks(x = baflogr$ntotsm, thresh = .00001, from = .25, to =1.75)

# debug(genotype_loci_2ndpass)
baflogr <- genotype_loci_2ndpass(baflogr = baflogr, sampleid = sampleid,
                                  ntotonepeak = 1, ntotoneerr = 1/3, cnoneseg_minsnps = 100,
                                  nminzeropeak = 0, nminonepeak = 1, nminzeroerr = 1/100, nminzeroseg_minsnps = 1000, 
                                  plotting = T)

outdir1 <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/", sampleid, "/ascatrun1")
dir.create(outdir1, recursive = T, showWarnings = F)
setwd(outdir1)
# debug(ascat_run1)
ascat.output <- ascat_run1(baflogr = baflogr, sampleid = sampleid)

setwd("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/")
# redo na, nb, ntot calculations using latest ASCAT rho/psi and get final genotypes/input
# debug(genotype_loci_finalpass)
baflogr <- genotype_loci_finalpass(baflogr = baflogr, sampleid = sampleid,
                                   lohseg_minlength = 1e5,
                                   ascatoutput = ascat.output,
                                   plotting = T)

####### run the above on every single tumour-host pair
####### collect the set of all loci which are */AB (i.e. het in tumour) in some low percent of all cases
####### add per tumour all loci in regions which are LOH in "all" cases (modify code below)
####### this creates a common set to run on which should, given # samples -> Inf go to a generally applicable set

outdir2 <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/", sampleid, "/ascatrun2")
dir.create(outdir2, recursive = T, showWarnings = F)
setwd(outdir2)
ascat.output <- ascat_run2(baflogr = baflogr, ascatoutput1 = ascat.output, sampleid = sampleid)

