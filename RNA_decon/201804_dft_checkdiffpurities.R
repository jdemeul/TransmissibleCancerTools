### check solution to 25%, 50% and 75% purity sample

diffpuritydf <- comparisondf
diffpuritydf$eDFT2n <- diffpuritydf$eDFT2/rowSums(diffpuritydf[, c("nAT", "nBT")])
diffpuritydf$eFibron <- diffpuritydf$eFibro/rowSums(diffpuritydf[, c("nAH", "nBH")])

p1 <- ggplot(data = diffpuritydf, mapping = aes(x = log10(eN25+1), y = log(eN+1))) + geom_point(alpha = .50)
p1

p1 <- ggplot(data = diffpuritydf, mapping = aes(x = log10(Good_depth_50+1))) + geom_point(mapping = aes(y = eT/(eT+eN)), alpha = .50)
p1

p1 <- ggplot(data = diffpuritydf, mapping = aes(x = eT75/(eT75+eN75), y = eDFT2n/(eDFT2n+eFibron))) + geom_point(alpha = .5) + geom_density2d()
p1

p1 <- ggplot(data = diffpuritydf, mapping = aes(x = eDFT2n/(eDFT2n+eFibron))) + geom_histogram(bins = 100)
p1
