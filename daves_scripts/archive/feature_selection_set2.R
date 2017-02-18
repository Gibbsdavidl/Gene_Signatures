

# Find representative features by minimizing the within-set #
# Pearson correlation and using sampling distributions      #
# that are proportional to distance (1-corr).               #

# dgibbs@systemsbiology.org

require(ggplot2)
require(compiler)
require(parallel)
library(pheatmap)

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


getAvgCor <- function(idx, datCor) {
  interPtCor <- abs(datCor[idx,idx])
  diag(interPtCor) <- 0
  avgCor <- mean(interPtCor[upper.tri(interPtCor)])
}

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

resultsName <- "set2"
dat <- set2_1
datm <- t(dat[,-c(1:2)])
datmcor <- cor(datm, method="spearman")
reps <- 500
iters <- 1000
cores <- 5
kmax <- 22

meanCor1 <- c()
meanCor2 <- c()

for (k in 4:kmax) {
  print(paste("Working on k=",k,sep=""))

  # split the data into 2 parts
  idx1 <- sample(1:nrow(datm), size=floor(nrow(datm)/2), replace=F)
  idx2 <- (1:nrow(datm))[-idx1]
  datmcor1 <- cor(datm[idx1,], method="spearman")
  datmcor2 <- cor(datm[idx2,], method="spearman")
  soln1 <- multiSel(datm[idx1,], datmcor1, reps, iters, k, cores)
  print("    Solution 1 done.")
  soln2 <- multiSel(datm[idx2,], datmcor2, reps, iters, k, cores)
  print("    Solution 2 done.")

  # where we're going to write our data
  resDir <- paste("results/",resultsName,"/k",k,"/",sep="")
  print(paste("    Writing to: ", resDir))
  if (!dir.exists(resDir)) {dir.create(resDir,recursive=TRUE)}

  # what's the solution?
  sidx1 <- as.numeric(names(soln1[[1]]))
  sidx2 <- as.numeric(names(soln2[[1]]))
  write.table(dat[sidx1, 1:2], file=paste(resDir,"feature_selection_soln1.txt",sep=""), quote=F)
  write.table(dat[sidx2, 1:2], file=paste(resDir,"feature_selection_soln2.txt",sep=""), quote=F)

  # track the average inter-set correlation
  meanCor1 <- c(meanCor1, getAvgCor(sidx1, datmcor))
  meanCor2 <- c(meanCor2, getAvgCor(sidx2, datmcor))

  ## Plotting ##
  print("    Writing output.")
  pdf(paste(resDir,"selected_correlation_heatmap_1.pdf",sep=""), onefile=F)
  pheatmap(datmcor[sidx1,sidx1])
  dev.off()
  pdf(paste(resDir,"selected_correlation_heatmap_2.pdf",sep=""), onefile=F)
  pheatmap(datmcor[sidx2,sidx2])
  dev.off()

  # plotting the number of times a feature was selected.
  jdx <- as.numeric(names(soln1[[2]]))
  x <- data.frame(dat[jdx,c(1,2)], Scores=as.numeric(soln1[[2]]))
  g <- ggplot(data=x, aes(x=factor(SetName), y=Scores)) +#
       geom_bar(stat="identity")+
       theme(axis.text.x = element_text(angle = 50, hjust = 1))
  ggsave(paste(resDir,"selected_barplot_soln_1.pdf",sep=""), g, width = 8, height = 8)

  jdx <- as.numeric(names(soln2[[2]]))
  x <- data.frame(dat[jdx,c(1,2)], Scores=as.numeric(soln2[[2]]))
  g <- ggplot(data=x, aes(x=factor(SetName), y=Scores)) +#
       geom_bar(stat="identity")+
       theme(axis.text.x = element_text(angle = 50, hjust = 1))
  ggsave(paste(resDir,"selected_barplot_soln_2.pdf",sep=""), g, width = 8, height = 8)

} # end k loop

df <- data.frame(Ks=4:kmax, Corrs1=meanCor1, Corrs2=meanCor2)
write.table(df, file=paste("results/",resultsName,"/K_corr_table.txt",sep=""), quote=F)

print("DONE")
