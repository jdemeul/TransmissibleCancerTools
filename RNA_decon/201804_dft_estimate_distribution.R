### check how error on ratio behaves

library(VGAM)

# get_deconv_distrib <- function(nAT, nBT, nAH, nBH, rho, eA, eB, size = 1E5) {
#   # xseq <- seq(stepsize, 1-stepsize, stepsize)
#   betaout <- rbeta(n = size, shape1 = eA+1, shape2 = eB+1)
#   outsample <- ( 1 - ((nAH+nBH)*rho*(nAT-(nAT+nBT)*betaout)) / ((nAT+nBT)*(1-rho)*(nAH-(nAH+nBH)*betaout) ) )^-1
#   return(data.frame(dens = outsample))
# }

get_deconv_distrib <- function(nAT, nBT, nAH, nBH, rho, eA, eB, size = 1E6, quant = c(.025, .975)) {
  Nt <- eA+eB
  eAm <- rbetabinom.ab(n = size, size = Nt, shape1 = eA+1, shape2 = eB+1)
  eT <- ((nAH+nBH)*eAm - nAH*Nt) / (rho*(nBH*nAT-nAH*nBT))
  eT[eT < 0 | eT*(nAT+nBT)*rho > Nt] <- NA
  if (all(is.na(eT))) {
    return(data.frame(eT = NA, eH = NA, ratio_man = NA))
  }
  if (!is.null(quant)) {
    mleT <- as.numeric(names(sort(table(eT), decreasing = T))[1])
    qeT <- quantile(x = eT, probs = c(.025,.975), na.rm = T)
    eT <- c(ml = mleT, qeT)
  }
  nTeT <- eT*(nAT+nBT)*rho
  nHeH <- Nt - nTeT
  eH <- nHeH/((nAH+nBH)*(1-rho))
  outratio_man <- nTeT/Nt
  # outratio <- ( 1 - ((nAH+nBH)*rho*(nAT*Nt-(nAT+nBT)*eAm)) / ((nAT+nBT)*(1-rho)*(nAH*Nt-(nAH+nBH)*eAm )) )^-1
  # outratio <- ((nAT+nBT)*((nAH+nBH)*eAm - nAH*Nt)) / (Nt*(nBH*nAT-nAH*nBT))
  # outsample[outsample > 1 | outsample < 0] <- NA
  return(data.frame(eT = eT, eH = eH, ratio_man = outratio_man))
}


df1 <- get_deconv_distrib(nAT = 1, nBT = 1, nAH = 0, nBH = 2, rho = .5, eA = 50, eB = 69, quant = NULL)

quantile(x = df1$ratio_man, probs = c(.025,.975), na.rm = T)

library(ggplot2)
p1 <- ggplot(data = df1, mapping = aes(x = eT)) + geom_histogram(binwidth = 1, mapping = aes(y = ..density..), fill = "red") + geom_histogram(binwidth = 1, mapping = aes(x = eH, y = ..density..), fill = "blue")
p1

(3*79.33333)/(3*79.3333+2*60)

solve_sys <- function(nAT, nBT, nAH, nBH, rho, eA, eB) {
  eT <- (nBH*eA-nAH*eB)/((nBH*nAT-nAH*nBT)*rho)
  eH <- (nAT*eB-nBT*eA)/((nBH*nAT-nAH*nBT)*(1-rho))
  return(c(eT = eT, eH = eH))
}

solve_sys(nAT = 1, nBT = 1, nAH = 0, nBH = 2, rho = .5, eA = 79, eB = 69)
