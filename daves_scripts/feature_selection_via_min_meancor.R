

# Find representative features by minimizing the within-set #
# Pearson correlation and using sampling distributions      #
# that are proportional to distance (1-corr).               #

# dgibbs@systemsbiology.org

require(ggplot2)
require(compiler)
require(parallel)

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
# on a subsample


load("Signatures_for_Subtyping.rda")

dat <- set1_1
datm <- t(dat[,-c(1:2)])
reps <- 100
iters <- 1000
cores <- 4
k <- 4

datrcor <- cor(datm, method="spearman")
idx1 <- sample(1:nrow(datm), size=floor(nrow(datm)/2), replace=F)
idx2 <- (1:nrow(datm))[-idx1]
soln1 <- multiSel(datm[idx1,], datrcor, reps, iters, k, cores)
soln2 <- multiSel(datm[idx2,], datrcor, reps, iters, k, cores)

# what's the solution?
sidx1 <- as.numeric(names(soln1[[1]]))
sidx2 <- as.numeric(names(soln2[[1]]))
print(dat[sidx1, 1:2])
print(dat[sidx2, 1:2])
write.table(dat[sidx1, 1:2], file="feature_selection_soln1.txt", quote=F)
write.table(dat[sidx2, 1:2], file="feature_selection_soln2.txt", quote=F)

################# Plotting ##
#library(pheatmap)
#pdf("selected_correlation_heatmap.pdf")
#pheatmap(set1spear[idx,idx])
#dev.off()

# plotting the number of times a feature was selected.
jdx <- as.numeric(names(soln1[[2]]))
x <- data.frame(dat[jdx,c(1,2)], Scores=as.numeric(soln1[[2]]))
pdf("selected_barplot_soln1.pdf")
ggplot(data=x, aes(x=factor(SetName), y=Scores)) +#
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
dev.off()
################### End
