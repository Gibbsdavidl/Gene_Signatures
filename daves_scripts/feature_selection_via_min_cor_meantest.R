

# Find representative features by minimizing the within-set #
# Pearson correlation and using sampling distributions      #
# that are proportional to distance (1-corr).               #

# dgibbs@systemsbiology.org

require(ggplot2)
require(compiler)
require(parallel)

load("Signatures_for_Subtyping.rda")
set1 <- t(set1_1[,-c(1:2)])
set1spear <- cor(set1, method="spearman")

getSampDist <- function(idx, datCor) {
  # returns the sampling distribution based on the
  # median absolute distance between the selected set,
  # and all others
#  meanCorDist <- 1-abs(apply(set1spear[idx,], 2, median))
  meanCorDist <- 1-apply(abs(datCor)[idx,], 2, mean)
  meanCorDist[idx] <- 0
  meanCorDist/sum(meanCorDist)
}
getSampDistCmp <- cmpfun(getSampDist)


checkForImprovement <- function(i, idx, j, datCor, k) {
  # compare the median of pairwise correlations,
  # if the newly sampled feature reduced the median,
  # then take it.
  interPtCor <- abs(datCor[idx,idx])
  diag(interPtCor) <- 0
  avgCor <- mean(interPtCor[upper.tri(interPtCor)])
  avgCorVec <- apply(interPtCor, 2, mean)

  # remove the feature with the highest correlation
  l <- sample(1:k, size=1, prob=avgCorVec/sum(avgCorVec))

  jdx <- c(idx[-l], j) # try new index
  newInterPtCor <- abs(datCor[jdx,jdx])
  diag(newInterPtCor) <- 0
  newAvgCor <- mean(newInterPtCor[upper.tri(newInterPtCor)])

  # if we reduced the median correlation, then take it.
  if (newAvgCor < avgCor) { idx <- jdx; }
  idx
}
checkForImprovementCmp <- cmpfun(checkForImprovement)


featSel <- function(dat, datCor, iter, k) {
  # minimize pairwise correlation
  print(timestamp())

  # the list of selected features
  # start with two random features
  idx <- sample(1:ncol(dat), size=k, replace=F)

  for (i in 1:iter) {
    # first sample a distant point
    sampDist <- getSampDistCmp(idx, datCor)
    j <- sample(1:ncol(dat), size=1, replace=F, prob=sampDist)

    # if j is better, keep it.
    idx <- checkForImprovementCmp(i, idx, j, datCor, k)
  }

  return(idx)
}
featSelCmp <- cmpfun(featSel)


multiSel <- function(dat, datCor, reps, iter, k, cores) {

  # do feature selection `reps` number of times
  res0 <- mclapply(1:reps, function(a) featSelCmp(dat, datCor, iter, k), mc.cores=cores)

  # count number of times each feature was selected
  scores <- table(unlist(res0))

  # return the k best
  bestScores <- sort(scores, decreasing=T)[1:k]

  list(bestScores, scores)
}

################################
# run it a number of times... ##
soln <- multiSel(set1, set1spear, 100, 1000, 4, 4)

# what's the solution?
idx <- as.numeric(names(soln[[1]]))
print(set1_1[idx, 1:2])
write.table(set1_1[idx, 1:2], file="feature_selection.txt", quote=F)

################# Plotting ##
library(pheatmap)
pdf("selected_correlation_heatmap.pdf")
pheatmap(set1spear[idx,idx])
dev.off()

# plotting the number of times a feature was selected.
jdx <- as.numeric(names(soln[[2]]))
x <- data.frame(set1_1[jdx,c(1,2)], Scores=as.numeric(soln[[2]]))
pdf("selected_barplot.pdf")
ggplot(data=x, aes(x=factor(SetName), y=Scores)) +#
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
dev.off()
################### End
