### CN calling pipeline Transmissible cancers


get_devil_genome <- function(devilinfofile) {
  devilinfo <- readr::read_tsv(file = devilinfo_file, col_names = c("contig", "width"), col_types = "ci")
  devilinfo$chr <- tolower(substr(x = devilinfo$contig, start = 1, stop = 4))
  devilinfo$start <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = function(x) cumsum(c(1,x))[-length(x)]))
  devilinfo$end <- unlist(by(data = as.numeric(devilinfo$width), INDICES = devilinfo$chr, FUN = cumsum))
  devilinfo_gr <- GRanges(seqnames = devilinfo$chr, ranges = IRanges(start = devilinfo$start, end = devilinfo$end), contig = devilinfo$contig)
  return(devilinfo_gr)  
}


get_baf_logr <- function(tumour_allelecounts_file, host_allelecounts_file, reference_alleles_file, genomeinfo, min_contig_size = 1e6, min_depth = 10) {
  
  reference_alleles <- readr::read_tsv(file = reference_alleles_file, col_names = c("scaffold", "pos", "ref", "alt"), col_types = "cicc-", skip = 1)
  
  tumour_allelecounts <- readr::read_tsv(file = tumour_allelecounts_file, col_names = c("scaffold", "pos", paste0("count_", c("A", "C", "G", "T"))), col_types = "ciiiii-", skip = 1)
  host_allelecounts <- readr::read_tsv(file = host_allelecounts_file, col_names = c("scaffold", "pos", paste0("count_", c("A", "C", "G", "T"))), col_types = "ciiiii-", skip = 1)
  
  baflogr <- reference_alleles
  baflogr$chr <- tolower(substr(baflogr$scaffold, 1, 4))
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
  
  # some filtering on minimal depth and unreliable contigs
  allowed_contigs <- mcols(genomeinfo)[width(genomeinfo) >= min_contig_size, "contig"]
  # require at least one read in the tumour and mindepth in normal
  has_good_depth <- (baflogr$h_ref + baflogr$h_alt >= min_depth) & (baflogr$t_ref + baflogr$t_alt >= 1)
  
  # on_sel_chrom <- dft2_dna$chr == "Chr4"
  
  baflogr <- baflogr[has_good_depth & baflogr$scaffold %in% allowed_contigs, ]
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
  p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos, y = ntot)) + geom_point(shape = ".", alpha = .2) + facet_wrap(~chr) + ylim(c(-.2, 6))
  # p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos, y = ntot)) + geom_point(shape = ".", alpha = .2) + facet_wrap(~chr) + ylim(c(-.2, 6))
  ggsave(filename = paste0(sampleid, "_initialRhoPsiT_totalCNplot.png"), plot = p1)
  return(p1)
}



genotype_loci <- function(baflogr, rho = .5, psit = 2.1) {
  
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
  
  baflogr$htgeno <- ifelse(baflogr$h_ref == 0 | hostbaf >= .95, ifelse(baflogr$t_ref >= 3, "BB/A*", "BB/*"), # host hom alt, ref reads tumour must come from cancer 
                           ifelse(baflogr$h_alt == 0 | hostbaf <= .05, ifelse(baflogr$t_alt >= 3, "AA/B*", "AA/*"), # host hom ref, alt reads tumour must come from cancer
                                  ifelse(baflogr$h_ref >= 3 & baflogr$h_alt >= 3 & ( hostbaf >= .25 | hostbaf <= .75 ), "AB/*", "*/*"))) #BB
  
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
  baflogr$nmin <- ifelse(baflogr$na <= baflogr$nb, baflogr$na, baflogr$nb)
  
  ## smooth minor and total CN
  # total CN should be stable, so running median
  baflogr$ntotsm <- unlist(by(data = baflogr$ntot, INDICES = baflogr$chr, FUN = function(x) runmed(x = x, k = 251, endrule = "constant")))
  # minor will oscillate between two states for any segment, reflecting hom/het SNPs, so use running mean to deect where truly zero, i.e. tumour LOH
  baflogr$nminsm <- unlist(by(data = baflogr$nmin, INDICES = baflogr$chr, FUN = function(x) as.numeric(runmean(x = Rle(x), k = 1001, endrule = "constant"))))
  
  
  # checking
  # testdf <- baflogr[sample(x = 1:nrow(baflogr), size = 100000, replace = F), ]
  # p1 <- ggplot(data = baflogr[baflogr$chr == "19",], mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = 1, mapping = aes(y = ntot), colour = "goldenrod2")
  p1 <- ggplot(data = baflogr, mapping = aes(x = cumpos)) + geom_point(shape = ".", alpha = .2, mapping = aes(y = ntot), colour = "goldenrod2")
  p1 <- p1 + geom_point(shape = ".", alpha = .2, mapping = aes(y = nmin), colour = "darkslategrey") + geom_point(mapping = aes(y = nminsm), colour = "red", alpha = .5, shape = ".")
  p1 <- p1 + geom_point(mapping = aes(y = ntotsm), colour = "black", alpha = .5, shape = ".")
  p1 <- p1 + facet_wrap(~chr)
  p1 <- p1 + ylim(c(-.2, 6))
  # p1 <- p1 + xlim(c(7.25e7,7.5e7))
  # p1 <- p1+ylim(c(0, 1))
  
  
  # ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2018_Murchison/results/20180426_CN-transformed_BAFLogR.png", plot = p1, width = 16, height = 6)
  
  # plot(1:length(cnsegs$Chr2$yhat), cnsegs$Chr2$yhat)
  
  # p1 <- ggplot(data = baflogr, mapping = aes(x = ntotsm, fill = htgeno)) + geom_histogram(binwidth = .01)
  p1 <- ggplot(data = baflogr[baflogr$htgeno != "AA/*",], mapping = aes(x = nmin, fill = htgeno)) + geom_histogram(binwidth = .01)
  # p1 <- p1 + geom_vline(xintercept = .1)
  p1 <- p1 + facet_wrap(~chr)
  p1 <- p1 + xlim(c(-.2, 5))
  p1
  
  return(baflogr)
}

