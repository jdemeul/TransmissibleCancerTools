## mod getgermlinegenotypes

#' @title ascat.predictGermlineGenotypes
#' @description predicts the germline genotypes of samples for which no matched germline sample
#' is available
#' @param ASCATobj an ASCAT object
#' @param platform used array platform
#' @param img.dir directory in which figures will be written
#' @param img.prefix prefix for figure names
#' @details Currently possible values for platform:\cr
#' AffySNP6 (default)\cr
#' Custom10k\cr
#' Illumina109k\cr
#' IlluminaCytoSNP\cr
#' Illumina610k\cr
#' Illumina660k\cr
#' Illumina700k\cr
#' Illumina1M\cr
#' Illumina2.5M\cr
#' IlluminaOmni5\cr
#' Affy10k\cr
#' Affy100k\cr
#' Affy250k_sty\cr
#' Affy250k_nsp\cr
#' AffyOncoScan\cr
#' AffyCytoScanHD\cr
#' HumanCNV370quad\cr
#' HumanCore12\cr
#' HumanCoreExome24\cr
#' HumanOmniExpress12\cr
#' IlluminaOmniExpressExome\cr
#'
#' @return predicted germline genotypes
#'
#' @export
ascat.predictGermlineGenotypes.mod = function(ASCATobj, platform = "AffySNP6", img.dir=".", img.prefix="") {
  Homozygous = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  colnames(Homozygous) = colnames(ASCATobj$Tumor_LogR)
  rownames(Homozygous) = rownames(ASCATobj$Tumor_LogR)
  
  if (platform=="DevilNGS") {
    maxHomozygous = 0.05
    proportionHetero = 0.33
    proportionHomo = 0.59
    proportionOpen = 0.04
    segmentLength = 100
  } else {
    print("Error: platform unknown")
  }
  
  failedarrays = NULL
  
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    
    Tumor_BAF_noNA = ASCATobj$Tumor_BAF[!is.na(ASCATobj$Tumor_BAF[,i]),i]
    names(Tumor_BAF_noNA) = rownames(ASCATobj$Tumor_BAF)[!is.na(ASCATobj$Tumor_BAF[,i])]
    Tumor_LogR_noNA = ASCATobj$Tumor_LogR[names(Tumor_BAF_noNA),i]
    names(Tumor_LogR_noNA) = names(Tumor_BAF_noNA)
    
    chr_noNA = list()
    prev = 0
    for(j in 1:length(ASCATobj$chr)) {
      chrke = ASCATobj$chr[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,i]))
      chr_noNA[[j]] = (prev+1):next2
      prev = next2
    }
    
    ch_noNA = list()
    prev = 0
    for(j in 1:length(ASCATobj$ch)) {
      chrke = ASCATobj$ch[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,i]))
      ch_noNA[[j]] = (prev+1):next2
      prev = next2
    }
    
    tbsam = Tumor_BAF_noNA
    # sample, mirrored
    bsm = ifelse(tbsam<0.5, tbsam, 1-tbsam)
    
    homoLimit = max(sort(bsm)[round(length(bsm)*proportionHomo)],maxHomozygous)
    
    if(homoLimit>0.25) {
      failedarrays = c(failedarrays,ASCATobj$samples[i])
    }
    
    Hom = ifelse(bsm<homoLimit,T,NA)
    
    Homo = sum(Hom==T, na.rm=T)
    Undecided = sum(is.na(Hom))
    
    extraHetero = round(min(proportionHetero * length(Tumor_BAF_noNA),Undecided-proportionOpen*length(Tumor_BAF_noNA)))
    
    Hetero = 0
    
    if(extraHetero>0) {
      
      allProbes=1:length(Tumor_BAF_noNA)
      nonHomoProbes = allProbes[is.na(Hom)|Hom==F]
      
      lowestDist = NULL
      
      # bsm, with homozygous replaced by NA
      bsmHNA=bsm
      bsmHNA[!is.na(Hom)&Hom]=NA
      
      for (chrke in chr_noNA) {
        
        chrNonHomoProbes = intersect(nonHomoProbes,chrke)
        
        # there must be a minimum number of probes on the chromosome, otherwise these are called homozygous anyway
        if (length(chrNonHomoProbes)>5) {
          
          #make sure we're not going over any borders..
          segmentLength2 = min(length(chrNonHomoProbes)-1,segmentLength)
          
          chrNonHomoProbesStartWindowLeft = c(rep(NA,segmentLength2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2)])
          chrNonHomoProbesEndWindowLeft = c(NA,chrNonHomoProbes[1:(length(chrNonHomoProbes)-1)])
          chrNonHomoProbesStartWindowRight = c(chrNonHomoProbes[2:length(chrNonHomoProbes)],NA)
          chrNonHomoProbesEndWindowRight = c(chrNonHomoProbes[(segmentLength2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2))
          chrNonHomoProbesStartWindowMiddle = c(rep(NA,segmentLength2/2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2/2)])
          chrNonHomoProbesEndWindowMiddle = c(chrNonHomoProbes[(segmentLength2/2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2/2))
          
          chrLowestDist = NULL
          
          for (probeNr in 1:length(chrNonHomoProbes)) {
            probe = chrNonHomoProbes[probeNr]
            if(!is.na(chrNonHomoProbesStartWindowLeft[probeNr])&!is.na(chrNonHomoProbesEndWindowLeft[probeNr])) {
              medianLeft = median(bsmHNA[chrNonHomoProbesStartWindowLeft[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], na.rm=T)
            }
            else {
              medianLeft = NA
            }
            if(!is.na(chrNonHomoProbesStartWindowRight[probeNr])&!is.na(chrNonHomoProbesEndWindowRight[probeNr])) {
              medianRight = median(bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowRight[probeNr]], na.rm=T)
            }
            else {
              medianRight = NA
            }
            
            if(!is.na(chrNonHomoProbesStartWindowMiddle[probeNr])&!is.na(chrNonHomoProbesEndWindowMiddle[probeNr])) {
              medianMiddle = median(c(bsmHNA[chrNonHomoProbesStartWindowMiddle[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]],
                                      bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowMiddle[probeNr]]), na.rm=T)
            }
            else {
              medianMiddle = NA
            }
            
            chrLowestDist[probeNr] = min(abs(medianLeft-bsm[probe]),abs(medianRight-bsm[probe]),abs(medianMiddle-bsm[probe]),Inf,na.rm=T)
          }
        }
        
        # if too few probes on the chromosome
        else {
          chrLowestDist = NULL
          if (length(chrNonHomoProbes)>0) {
            # 1 is higher than any practical distance
            chrLowestDist[1:length(chrNonHomoProbes)] = 1
          }
        }
        
        lowestDist = c(lowestDist,chrLowestDist)
      }
      
      lowestDistUndecided = lowestDist[is.na(Hom[nonHomoProbes])]
      names(lowestDistUndecided)=names(Tumor_LogR_noNA)[nonHomoProbes[is.na(Hom[nonHomoProbes])]]
      
      sorted = sort(lowestDistUndecided)
      Hom[names(sorted[1:min(length(sorted),extraHetero)])] = F
      
      Hetero = sum(Hom==F, na.rm=T)
      Homo = sum(Hom==T, na.rm=T)
      Undecided = sum(is.na(Hom))
      
    }
    
    png(filename = file.path(img.dir,paste(img.prefix, "tumorSep",colnames(ASCATobj$Tumor_LogR)[i],".png",sep="")), width = 2000, height = 500, res = 200)
    title = paste(paste(colnames(ASCATobj$Tumor_BAF)[i], Hetero), Homo)
    ascat.plotGenotypes(ASCATobj,title,Tumor_BAF_noNA, Hom, ch_noNA)
    dev.off()
    
    # set all Undecided to homozygous
    Hom[is.na(Hom)] = T
    Homozygous[names(Hom),i] = Hom
  }
  
  return(list(germlinegenotypes = Homozygous, failedarrays = failedarrays))
  
}
