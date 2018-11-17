### CN calling pipeline Transmissible cancers


get_devil_genome <- function(devilinfofile) {
  devilinfo <- readr::read_tsv(file = devilinfo_file, col_names = c("contig", "width"), col_types = "ci")
  devilinfo$chr <- toupper(substr(x = devilinfo$contig, start = 4, stop = 4))
  devilinfo$start <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = function(x) cumsum(c(1,x))[-length(x)]))
  devilinfo$end <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = cumsum))
  devilinfo_gr <- GRanges(seqnames = devilinfo$chr, ranges = IRanges(start = devilinfo$start, end = devilinfo$end), contig = devilinfo$contig)
  return(devilinfo_gr)  
}


get_baf_logr <- function(tumour_allelecounts_file, host_allelecounts_file, reference_alleles_file, genomeinfo, min_contig_size = 1e6) {
  
  reference_alleles <- readr::read_tsv(file = reference_alleles_file, col_names = c("scaffold", "pos", "ref", "alt"), col_types = "cicc-", skip = 1)
  
  tumour_allelecounts <- readr::read_tsv(file = tumour_allelecounts_file, col_names = c("scaffold", "pos", paste0("count_", c("A", "C", "G", "T"))), col_types = "ciiiii-", skip = 1)
  host_allelecounts <- readr::read_tsv(file = host_allelecounts_file, col_names = c("scaffold", "pos", paste0("count_", c("A", "C", "G", "T"))), col_types = "ciiiii-", skip = 1)
  
  baflogr <- reference_alleles
  baflogr$chr <- toupper(substr(baflogr$scaffold, 4, 4))
  baflogr$cumpos <- start(genomeinfo[match(x = baflogr$scaffold, table = mcols(genomeinfo)$contig)]) + baflogr$pos
  baflogr$h_ref <- ifelse(reference_alleles$ref == "A", host_allelecounts$count_A,
                          ifelse(reference_alleles$ref == "C", host_allelecounts$count_C, 
                                 ifelse(reference_alleles$ref == "G", host_allelecounts$count_G, host_allelecounts$count_T)))
  baflogr$h_alt <- ifelse(reference_alleles$alt == "A", host_allelecounts$count_A,
                          ifelse(reference_alleles$alt == "C", host_allelecounts$count_C, 
                                 ifelse(reference_alleles$alt == "G", host_allelecounts$count_G, host_allelecounts$count_T)))
  
  baflogr$t_ref <- ifelse(reference_alleles$ref == "A", tumour_allelecounts$count_A,
                          ifelse(reference_alleles$ref == "C", tumour_allelecounts$count_C, 
                                 ifelse(reference_alleles$ref == "G", tumour_allelecounts$count_G, tumour_allelecounts$count_T)))
  baflogr$t_alt <- ifelse(reference_alleles$alt == "A", tumour_allelecounts$count_A,
                          ifelse(reference_alleles$alt == "C", tumour_allelecounts$count_C, 
                                 ifelse(reference_alleles$alt == "G", tumour_allelecounts$count_G, tumour_allelecounts$count_T)))
  
  # some filtering on unreliable contigs
  allowed_contigs <- mcols(genomeinfo)[width(genomeinfo) >= min_contig_size, "contig"]
  # # require at least one read in the tumour and mindepth in normal
  # has_good_depth <- (baflogr$h_ref + baflogr$h_alt >= min_depth) & (baflogr$t_ref + baflogr$t_alt >= 1)
  baflogr <- baflogr[which(baflogr$scaffold %in% allowed_contigs), ]
  
  # dft2_dna <- dft2_dna[has_good_depth & dft2_dna$X.CHR %in% allowed_contigs, ]
  # fibro_dna <- fibro_dna[has_good_depth & fibro_dna$X.CHR %in% allowed_contigs, ]
  
  totalhost <- baflogr$h_ref + baflogr$h_alt
  totaltum <- baflogr$t_ref + baflogr$t_alt
  
  baflogr$logr <- log2( (totaltum/totalhost) / median((totaltum/totalhost), na.rm = T) )
  baflogr$baf <- baflogr$t_alt/(baflogr$t_alt + baflogr$t_ref)
  
  return(baflogr)
}


get_baf_logr_kg_format <- function(datafile, sampleid, genomeinfo, min_contig_size = 1e6) {
  
  inputdata <- readr::read_tsv(file = datafile)
  
  baflogr <- inputdata[, c("CONTIG", "OFFSET", "CHROM", "POS", "REF", "ALT", 
                           grep(pattern = paste0(sampleid, "H.*nr$"), x = colnames(inputdata), value = T),
                           grep(pattern = paste0(sampleid, "H.*nv$"), x = colnames(inputdata), value = T),
                           grep(pattern = paste0(sampleid, "T.*nr$"), x = colnames(inputdata), value = T),
                           grep(pattern = paste0(sampleid, "T.*nv$"), x = colnames(inputdata), value = T))]
  colnames(baflogr) <- c("scaffold", "pos", "chr", "cumpos", "ref", "alt", "h_ref", "h_alt", "t_ref", "t_alt")
  baflogr$chr <- sub(pattern = "Chr", replacement = "", x = baflogr$chr)
  baflogr$h_ref <- baflogr$h_ref - baflogr$h_alt
  baflogr$t_ref <- baflogr$t_ref - baflogr$t_alt
  baflogr <- baflogr[!is.na(baflogr$cumpos),]
  # baflogr$chr <- toupper(substr(baflogr$scaffold, 4, 4))
  # baflogr$cumpos <- start(genomeinfo[match(x = baflogr$scaffold, table = mcols(genomeinfo)$contig)]) + baflogr$pos
  # baflogr$h_ref <- ifelse(reference_alleles$ref == "A", host_allelecounts$count_A,
  #                         ifelse(reference_alleles$ref == "C", host_allelecounts$count_C, 
  #                                ifelse(reference_alleles$ref == "G", host_allelecounts$count_G, host_allelecounts$count_T)))
  # baflogr$h_alt <- ifelse(reference_alleles$alt == "A", host_allelecounts$count_A,
  #                         ifelse(reference_alleles$alt == "C", host_allelecounts$count_C, 
  #                                ifelse(reference_alleles$alt == "G", host_allelecounts$count_G, host_allelecounts$count_T)))
  # 
  # baflogr$t_ref <- ifelse(reference_alleles$ref == "A", tumour_allelecounts$count_A,
  #                         ifelse(reference_alleles$ref == "C", tumour_allelecounts$count_C, 
  #                                ifelse(reference_alleles$ref == "G", tumour_allelecounts$count_G, tumour_allelecounts$count_T)))
  # baflogr$t_alt <- ifelse(reference_alleles$alt == "A", tumour_allelecounts$count_A,
  #                         ifelse(reference_alleles$alt == "C", tumour_allelecounts$count_C, 
  #                                ifelse(reference_alleles$alt == "G", tumour_allelecounts$count_G, tumour_allelecounts$count_T)))
  
  # some filtering on unreliable contigs
  allowed_contigs <- mcols(genomeinfo)[width(genomeinfo) >= min_contig_size, "contig"]
  # # require at least one read in the tumour and mindepth in normal
  # has_good_depth <- (baflogr$h_ref + baflogr$h_alt >= min_depth) & (baflogr$t_ref + baflogr$t_alt >= 1)
  baflogr <- baflogr[which(baflogr$scaffold %in% allowed_contigs), ]
  
  # dft2_dna <- dft2_dna[has_good_depth & dft2_dna$X.CHR %in% allowed_contigs, ]
  # fibro_dna <- fibro_dna[has_good_depth & fibro_dna$X.CHR %in% allowed_contigs, ]
  
  totalhost <- baflogr$h_ref + baflogr$h_alt
  totaltum <- baflogr$t_ref + baflogr$t_alt
  
  baflogr$logr <- log2( (totaltum/totalhost) / median((totaltum/totalhost), na.rm = T) )
  baflogr$baf <- baflogr$t_alt/(baflogr$t_alt + baflogr$t_ref)
  
  return(baflogr)
}


check_rho_psi_estimates <- function(sampleid, baflogr, rho = .5, psit = 1.85) {
  baflogr$ntot <- ((2*(1-rho) + rho*psit)*2^baflogr$logr - 2*(1-rho))/rho
  # p1 <- ggplot(data = baflogr[sample(x = 1:nrow(baflogr), size = 1e6, replace = F), ], mapping = aes(x = cumpos, y = ntot)) + geom_point(shape = ".", alpha = .2) + facet_wrap(~chr) + ylim(c(-.2, 6))
  p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos, y = ntot)) + geom_point(shape = ".", alpha = .2) + facet_wrap(~chr) + ylim(c(-.2, 6)) + theme_minimal()
  # p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos, y = ntot)) + geom_point(shape = ".", alpha = .2) + facet_wrap(~chr) + ylim(c(-.2, 6))
  ggsave(filename = paste0(sampleid, "_initialRhoPsiT_totalCNplot.png"), plot = p1)
  return(NULL)
}



genotype_loci <- function(baflogr, rho = .5, psit = 2.1, sampleid, species = "devil", min_depth = 10, plotting = T) {
  
  ## Reasoning host loci:
  # host = hom if ref == 0 | alt== 0 (coverage normal ≥ 10, imposed in previous step mindepth)
  # host = het if ≥ 3 reads for both ALT and REF & .25 ≤ BAF ≤ .75
  
  ## Reasoning bulk tumour:
  # if host is hom alt|ref and the other allele is observed in the tumour, we can infer AA/B* or BB/A* genotype Host/Cancer
  # if host is het, no info at this point on genotype cancer: AB/*
  hostbaf <- baflogr$h_alt/(baflogr$h_alt+baflogr$h_ref)
  # baflogr$htgeno <- ifelse(baflogr$h_ref == 0 & baflogr$t_ref >= 3, "BB/A*", # host hom alt, ref reads tumour must come from cancer 
  #                          ifelse(baflogr$h_alt  == 0 & baflogr$t_alt >= 3, "AA/B*", # host hom ref, alt reads tumour must come from cancer
  #                                 ifelse(baflogr$h_ref >= 3 & baflogr$h_alt >= 3 & ( hostbaf >= .25 | hostbaf <= .75 ), "AB/*", NA)))
  
  baflogr$htgeno <- ifelse(baflogr$h_ref + baflogr$h_alt < min_depth, "*/*", # can we make a proper call, i.e. do we have enough coverage in the normal?
                           ifelse(baflogr$h_ref == 0 | hostbaf >= .95, ifelse(baflogr$t_ref >= 3, "BB/A*", "BB/*"), # host hom alt, ref reads tumour must come from cancer 
                           ifelse(baflogr$h_alt == 0 | hostbaf <= .05, ifelse(baflogr$t_alt >= 3, "AA/B*", "AA/*"), # host hom ref, alt reads tumour must come from cancer
                                  ifelse(baflogr$h_ref >= 3 & baflogr$h_alt >= 3 & ( hostbaf >= .25 | hostbaf <= .75 ), "AB/*", "*/*")))) #BB
  
  ### omit loci where host genotype can not clearly be assigned?
  # baflogr <- baflogr[!is.na(baflogr$htgeno), ]
  
  ## Reconstruct tumour-specific copy number of the A and B alleles using rho/psit estimates, Host/Cancer genotypes, BAF and LogR
  # see derivation of equations in notes to correct for number of copies present in the host
  # basically total CN - # copies contributed by host
  baflogr$na <- ifelse(baflogr$htgeno == "AB/*", ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr - (1-rho))/rho, # single copy of A in host
                       ifelse(baflogr$htgeno %in% c("AA/B*", "AA/*"), ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr - 2*(1-rho))/rho, # two copies of A in host
                              ifelse(baflogr$htgeno %in% c("BB/A*", "BB/*"), ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr)/rho, NA))) # 0 copies of A in host
  
  # same for the B allele
  baflogr$nb <- ifelse(baflogr$htgeno == "AB/*", ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr - (1-rho))/rho, # single copy of B in host
                       ifelse(baflogr$htgeno %in% c("AA/B*", "AA/*"), ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr)/rho, # 0 copies of B in host
                              ifelse(baflogr$htgeno %in% c("BB/A*", "BB/*"), ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr - 2*(1-rho))/rho, NA))) # two copies of B in host
  
  # get CN minor allele
  baflogr$nmin <- ifelse(baflogr$na < baflogr$nb, baflogr$na, baflogr$nb)
  # and total
  baflogr$ntot <- (psit*2^baflogr$logr-2*(1-rho))/rho
  
  ## smooth minor and total CN
  # total CN should be stable, so running median
  # the data needs to be sorted chr-pos-wise for this and NA's need to be skipped
  notna_idxs <- which(!is.na(baflogr$na))
  baflogr$ntotsm <- NA
  baflogr$nminsm <- NA
  baflogr$ntotsm[notna_idxs] <- unlist(by(data = baflogr$ntot[notna_idxs], INDICES = baflogr$chr[notna_idxs], FUN = function(x) runmed(x = x, k = 251, endrule = "constant")))
  # minor will oscillate between two states for any segment, reflecting hom/het SNPs, so use running mean to deect where truly zero, i.e. tumour LOH
  baflogr$nminsm[notna_idxs] <- unlist(by(data = baflogr$nmin[notna_idxs], INDICES = baflogr$chr[notna_idxs], FUN = function(x) as.numeric(runmean(x = Rle(x), k = 1001, endrule = "constant"))))
  
  if (plotting) {
    # checking
    # testdf <- baflogr[sample(x = 1:nrow(baflogr), size = 100000, replace = F), ]
    # p1 <- ggplot(data = baflogr[baflogr$chr == "19",], mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = 1, mapping = aes(y = ntot), colour = "goldenrod2")
    p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = ntot), colour = "goldenrod2")
    p1 <- p1 + geom_point(shape = ".", alpha = .2, mapping = aes(y = nmin), colour = "darkslategrey")
    p1 <- p1 + geom_path(mapping = aes(y = nminsm), colour = "red", alpha = .5)
    p1 <- p1 + geom_path(mapping = aes(y = ntotsm), colour = "black", alpha = .5)
    p1 <- p1 + facet_wrap(~chr)
    p1 <- p1 + ylim(c(-.2, 6)) + theme_minimal()
    ggsave(filename = paste0(sampleid, "_nmin-ntot-LOH_plot.png"), plot = p1)
    # p1 <- p1 + xlim(c(7.25e7,7.5e7))
    # p1 <- p1+ylim(c(0, 1))
    
    # ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180426_CN-transformed_BAFLogR.png", plot = p1, width = 16, height = 6)
    
    # plot(1:length(cnsegs$Chr2$yhat), cnsegs$Chr2$yhat)
    
    # p1 <- ggplot(data = baflogr, mapping = aes(x = ntotsm, fill = htgeno)) + geom_histogram(binwidth = .01)
    
    p1 <- ggplot(data = baflogr, mapping = aes(x = nmin, fill = htgeno)) + geom_histogram(binwidth = .01)
    # p1 <- p1 + geom_vline(xintercept = .1)
    p1 <- p1 + facet_wrap(~chr)
    p1 <- p1 + lims(x = c(-.2, 5), y = c(0, 5000)) + theme_minimal()
    ggsave(filename = paste0(sampleid, "_nmin-histogram-genotypes.png"), plot = p1)
  }
  
  return(baflogr)
}


genotype_loci_2ndpass <- function(baflogr, sampleid,
                                  ntotonepeak = 1, ntotoneerr = 1/3, cnoneseg_minsnps = 100,
                                  nminzeropeak = 0, nminonepeak = 1, nminzeroerr = 1/100, nminzeroseg_minsnps = 1000, 
                                  plotting = T) {
  
  # find LOH based on total CN = 1 regions in tumour
  # NAs should be scattered about, only where low coverage in normal, no long stretches, so safe to set them as T
  cnonesegs <- by(data = baflogr$ntotsm, INDICES = baflogr$chr, FUN = function(x) as(object = (x <= ntotonepeak + ntotoneerr | is.na(x)), "Rle"))
  cnonesegs <- as(object = unlist(lapply(cnonesegs, FUN = function(x) { runValue(x)[runLength(x) < cnoneseg_minsnps] <- F; return(as.vector(x))})), "Rle")
  
  segstarts <- start(cnonesegs)[runValue(cnonesegs)]
  segends <- end(cnonesegs)[runValue(cnonesegs)]
  
  # find LOH based on minor CN = 0 regions in tumour
  nminzerosegs <- by(data = baflogr$nminsm, INDICES = baflogr$chr, FUN = function(x) as(object = (abs(x) <= nminzeropeak + nminzeroerr | is.na(x)), "Rle"))
  nminzerosegs <- as(object = unlist(lapply(nminzerosegs, FUN = function(x) { runValue(x)[runLength(x) < nminzeroseg_minsnps] <- F; return(as.vector(x))})), "Rle")
  
  segstarts2 <- start(nminzerosegs)[runValue(nminzerosegs)]
  segends2 <- end(nminzerosegs)[runValue(nminzerosegs)]
  
  if (plotting) {
    p1 <- ggplot(data = baflogr, mapping = aes(x = 1:nrow(baflogr))) + geom_point(shape = ".", alpha = .2, mapping = aes(y = ntot), colour = "goldenrod2")
    p1 <- p1 + geom_point(shape = ".", alpha = .2, mapping = aes(y = nmin), colour = "darkslategrey")
    if (length(segstarts2 > 0))
      p1 <- p1 + geom_segment(data = data.frame(segstarts2, segends2), mapping = aes(x = segstarts2, xend = segends2, y = 0, yend = 0), colour = "red", size = 5, alpha = .75)
    if (length(segstarts > 0))
      p1 <- p1 + geom_segment(data = data.frame(segstarts, segends), mapping = aes(x = segstarts, xend = segends, y = 1, yend = 1), colour = "black", size = 5, alpha = .75)
    p1 <- p1 + ylim(c(-.2, 6)) + theme_minimal()
    ggsave(filename = paste0(sampleid, "_identifiedLOHsegments.png"), plot = p1)
  }
  # ggsave(filename = paste0(sampleid, "_nmin-ntot-LOH_plot.png"), plot = p1)
  
  
  ### selecting SNPs to run on:
  # get all hetloci (nmin ≥ 1) which are het in host as well (so htgeno = AB/*)
  hetidxs <- which(baflogr$htgeno == "AB/*")
  # nminoneidxs <- which(baflogr$nmin >= (.75*nminonepeak + .25*nminzeropeak)/2)
  nminoneidxs <- which(baflogr$nmin >= (.66*nminonepeak + .33*nminzeropeak))
  cnoneidxs <- unlist(mapply(segstarts, segends, FUN = ':'))
  nminzeroidxs <- unlist(mapply(segstarts2, segends2, FUN = ':'))
  # hostAAidxs <- which(baflogr$htgeno %in% c('AA/B*', 'AA/*'))
  # hostBBidxs <- which(baflogr$htgeno %in% c('BB/A*', 'BB/*'))
  
  # finalidxs <- intersect(hetidxs, c(nminoneidxs, cnoneidxs, nminzeroidxs))
  
  baflogr$type <- NA
  baflogr[intersect(hetidxs, c(cnoneidxs, nminzeroidxs)), "type"] <- 'AB/A'
  baflogr[intersect(hetidxs, nminoneidxs), "type"] <- "AB/AB"
  # baflogr[intersect(hostAAidxs, nminoneidxs), "type"] <- 'AA/AB'
  # baflogr[intersect(hostBBidxs, nminoneidxs), "type"] <- 'BB/AB'
  
  # finalbaflogr <- baflogr[finalidxs, ]
  if (plotting) {
    # p1 <- ggplot(data = baflogr[baflogr$type %in% c("AB/AB", "AB/A") & baflogr$chr != "Chrx", ], mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = baf), colour = "goldenrod2")
    p1 <- ggplot(data = baflogr[baflogr$type %in% c("AB/AB", "AB/A") & baflogr$chr != "X", ], mapping = aes(x = cumpos))
    p1 <- p1 + geom_point(shape = ".", alpha = .2, mapping = aes(y = baf), colour = "goldenrod2")
    p1 <- p1 + facet_wrap(~chr) + theme_minimal()
    # p1 <- p1+ylim(c(-.2, 6))
    ggsave(filename = paste0(sampleid, "_BAFinputToASCATrun1.png"), plot = p1)
    # ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180426_input_BAF_first_ASCAT_run.png", plot = p1, width = 16, height = 6)
  }
  
  # write_tsv(x = baflogr, path = paste0(sampleid, "_genotyped_BAF-Logr_all.txt"), col_names = T)
  
  return(baflogr)
}



genotype_loci_finalpass <- function(baflogr, sampleid,
                                    lohseg_minlength = 1e5,
                                    ascatoutput,
                                    plotting = T) {
  
  # str(ascat.output)
  psit <- ascatoutput$psi
  rho <- ascatoutput$aberrantcellfraction
  
  # recompute na, nb, ntot with ascat rho/psi
  baflogr$na <- ifelse(baflogr$htgeno == "AB/*", ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr - (1-rho))/rho, # single copy of A in host
                       ifelse(baflogr$htgeno %in% c("AA/B*", "AA/*"), ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr - 2*(1-rho))/rho, # two copies of A in host
                              ifelse(baflogr$htgeno %in% c("BB/A*", "BB/*"), ((2*(1-rho) + rho*psit)*(1-baflogr$baf)*2^baflogr$logr)/rho, NA))) # 0 copies of A in host
  
  # same for the B allele
  baflogr$nb <- ifelse(baflogr$htgeno == "AB/*", ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr - (1-rho))/rho, # single copy of B in host
                       ifelse(baflogr$htgeno %in% c("AA/B*", "AA/*"), ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr)/rho, # 0 copies of B in host
                              ifelse(baflogr$htgeno %in% c("BB/A*", "BB/*"), ((2*(1-rho) + rho*psit)*baflogr$baf*2^baflogr$logr - 2*(1-rho))/rho, NA))) # two copies of B in host
  
  # get CN minor allele
  baflogr$nmin <- ifelse(baflogr$na < baflogr$nb, baflogr$na, baflogr$nb)
  # and total
  baflogr$ntot <- (psit*2^baflogr$logr-2*(1-rho))/rho
  
  # get LOH segments from ASCAT
  ascat_lohsegs <- ascatoutput$segments[ascatoutput$segments$nMinor == 0, ]
  # ascat_lohsegs$chr <- toupper(sub(pattern = "chr", replacement = "", x = ascat_lohsegs$chr))
  
  # plotting for quick check
  if (plotting) {
    p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = ntot), colour = "goldenrod2")
    p1 <- p1 + geom_point(shape = ".", alpha = .2, mapping = aes(y = nmin), colour = "darkslategrey")
    p1 <- p1 + geom_segment(data = ascat_lohsegs, mapping = aes(x = startpos, xend = endpos, y = nMajor, yend = nMajor), colour = "red", size = 5, alpha = .75)
    p1 <- p1 + ylim(c(-.2, 6)) + theme_minimal() + facet_wrap(~chr)
    ggsave(filename = paste0(sampleid, "_identifiedLOHsegments_ascatrun1.png"), plot = p1)
  }
  
  # create lohseg granges
  ascat_lohsegs_gr <- GRanges(seqnames = ascat_lohsegs$chr, ranges = IRanges(start = ascat_lohsegs$startpos, end = ascat_lohsegs$endpos))
  ascat_lohsegs_gr <- ascat_lohsegs_gr[width(ascat_lohsegs_gr) >= lohseg_minlength]
  
  # create baflogr granges and overlap with LOH segments
  baflogr_gr <- GRanges(seqnames = baflogr$chr, ranges = IRanges(start = baflogr$cumpos, end = baflogr$cumpos))
  lohidxs <- which(baflogr_gr %over% ascat_lohsegs_gr)
  
  
  # reset typing
  baflogr$type <- NA

  # get all hetloci (nmin ≥ 1) which are het in host as well (so htgeno = AB/*)
  hetidxs <- which(baflogr$htgeno == "AB/*")
  nminoneidxs <- which(baflogr$nmin >= .66)
  hostAAidxs <- which(baflogr$htgeno %in% c('AA/B*', 'AA/*'))
  hostBBidxs <- which(baflogr$htgeno %in% c('BB/A*', 'BB/*'))
  
  baflogr[intersect(hetidxs, lohidxs), "type"] <- 'AB/A'
  baflogr[intersect(hetidxs, nminoneidxs), "type"] <- "AB/AB"
  baflogr[intersect(hostAAidxs, nminoneidxs), "type"] <- 'AA/AB'
  baflogr[intersect(hostBBidxs, nminoneidxs), "type"] <- 'BB/AB'
  
  write_tsv(x = baflogr, path = paste0(sampleid, "_genotyped_BAF-Logr_all.txt"), col_names = T)
  
  return(baflogr)
}



ascat_run1 <- function(baflogr, sampleid) { 
  
  ## create BAF/logR data ASCAT run 1
  baflogr <- baflogr[baflogr$type %in% c("AB/AB", "AB/A"), ]
  # baflogr$chr <- toupper(sub(pattern = "chr", replacement = "", x = baflogr$chr))
  outdf <- data.frame(baflogr$chr, baflogr$cumpos, baflogr$logr)
  colnames(outdf) <- c("chr", "pos", paste0(sampleid, "T"))
  write.table(x = outdf, file = paste0(sampleid, "_tumour_LogR.txt"), sep = "\t", row.names = T, quote = F)
  
  outdf[,3] <- 0
  colnames(outdf)[3] <- paste0(sampleid, "N")
  write.table(x = outdf, file = paste0(sampleid, "_host_LogR.txt"), sep = "\t", row.names = T, quote = F)
  
  randvect <- sample(x = c(0,1), size = nrow(baflogr), replace = T)
  outdf <- data.frame(baflogr$chr, baflogr$cumpos, abs(baflogr$baf - randvect))
  colnames(outdf) <- c("chr", "pos", paste0(sampleid, "T"))
  write.table(x = outdf, file = paste0(sampleid, "_tumour_BAF.txt"), sep = "\t", row.names = T, quote = F)
  
  outdf[,3] <- abs(baflogr$h_alt/(baflogr$h_ref+baflogr$h_alt) - randvect)
  colnames(outdf)[3] <- paste0(sampleid, "N")
  write.table(x = outdf, file = paste0(sampleid, "_host_BAF.txt"), sep = "\t", row.names = T, quote = F)
  
  
  ### running ascat first time
  require(ASCAT)
  ascat.bc <- ascat.loadData(paste0(sampleid, "_tumour_LogR.txt"), paste0(sampleid, "_tumour_BAF.txt"),paste0(sampleid, "_host_LogR.txt"),paste0(sampleid, "_host_BAF.txt"))
  ascat.plotRawData(ascat.bc)
  
  ascat.bc <- ascat.aspcf(ascat.bc, penalty = 250)
  ascat.plotSegmentedData(ascat.bc)
  ascat.output = ascat.runAscat(ascat.bc, gamma = 1)
  # View(ascat.output$segments)
  write.table(x = ascat.output$segments, file = paste0(sampleid, ".ASCATprofile.segments.txt"), sep = "\t", quote = F, row.names = F)

  return(ascat.output)
}


ascat_run2 <- function(baflogr, ascatoutput1, sampleid) {
  
  psi <- ascatoutput1$psi
  rho <- ascatoutput1$aberrantcellfraction
  
  baflogr <- baflogr[!is.na(baflogr$type), ]
  # baflogr$chr <- toupper(sub(pattern = "chr", replacement = "", x = baflogr$chr))
  # baflogr_run2$baf <- ifelse(baflogr_run2$type %in% c("AA/AB", "BB/AB"), (baflogr_run2$nb*rho + (1-rho))/(2*(1-rho) + rho*baflogr_run2$ntot), baflogr_run2$baf)
  # recalibrate BAF for tumour het, host hom SNPs:
  baflogr$baf <- ifelse(baflogr$type %in% c("AA/AB"), (psi*2^baflogr$logr*baflogr$baf+(1-rho))/(psi*2^baflogr$logr),
                        ifelse(baflogr$type %in% c("BB/AB"), (psi*2^baflogr$logr*baflogr$baf-1*(1-rho))/(psi*2^baflogr$logr),
                               baflogr$baf))
  baflogr$normalbaf <- ifelse(baflogr$type %in% c("AA/AB", "BB/AB"), .5, baflogr$h_alt/(baflogr$h_ref+baflogr$h_alt))
  
  ### write BAF/logR data ASCAT run 2
  outdf <- data.frame(baflogr$chr, baflogr$cumpos, baflogr$logr)
  colnames(outdf) <- c("chr", "pos", paste0(sampleid, "T"))
  write.table(x = outdf, file = paste0(sampleid, "_tumour_LogR.txt"), sep = "\t", row.names = T, quote = F)
  
  outdf[, 3] <- 0
  colnames(outdf)[3] <- paste0(sampleid, "N")
  write.table(x = outdf, file = paste0(sampleid, "_host_LogR.txt"), sep = "\t", row.names = T, quote = F)
  
  randvect <- sample(x = c(0,1), size = nrow(baflogr), replace = T)
  outdf <- data.frame(baflogr$chr, baflogr$cumpos, abs(baflogr$baf - randvect))
  colnames(outdf) <- c("chr", "pos", paste0(sampleid, "T"))
  write.table(x = outdf, file = paste0(sampleid, "_tumour_BAF.txt"), sep = "\t", row.names = T, quote = F)
  
  outdf[,3] <- abs(baflogr$normalbaf - randvect)
  colnames(outdf)[3] <- paste0(sampleid, "N")
  write.table(x = outdf, file = paste0(sampleid, "_host_BAF.txt"), sep = "\t", row.names = T, quote = F)
  
  
  ### running ascat
  ascat.bc <- ascat.loadData(paste0(sampleid, "_tumour_LogR.txt"), paste0(sampleid, "_tumour_BAF.txt"),paste0(sampleid, "_host_LogR.txt"),paste0(sampleid, "_host_BAF.txt"))
  ascat.plotRawData(ascat.bc)
  
  ascat.bc <- ascat.aspcf(ascat.bc, penalty = 250)
  ascat.plotSegmentedData(ascat.bc)
  ascat.output = ascat.runAscat(ascat.bc, gamma = 1)
  # View(ascat.output$segments)
  write.table(x = ascat.output$segments, file = paste0(sampleid, ".ASCATprofile.segments.txt"), sep = "\t", quote = F, row.names = F)
  
  return(ascat.output)
}
