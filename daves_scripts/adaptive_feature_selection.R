

# Find representative features by minimizing the within-set #
# Pearson correlation and using sampling distributions      #
# that are proportional to distance (1-corr).               #

# dgibbs@systemsbiology.org

load("data/Signatures_for_Subtyping_2.rda")

resultsName <- "set3"
dat <- set3_2
datm <- t(dat[,-c(1:2)])
datmcor <- cor(datm, method="spearman")
reps <- 1000
iters <- 20000
cores <- 12
k <- 15
learnRate <- 0.001
params <- list(reps, iters, cores, k, learnRate, resultsName)


require(ggplot2)
require(compiler)
require(parallel)
library(pheatmap)


getCors <- function(idx, datCor) {
  # return the mean corrs for the inter-set
  interPtCor <- abs(datCor[idx,idx])
  diag(interPtCor) <- 0
  1-apply(interPtCor, 2, mean)
}
getCorsCmp <- cmpfun(getCors)


updateDist <- function(idx, kcors, sampDist, learnRate) {
  # add in the new corrs
  sampDist[idx] <- sampDist[idx]+kcors*learnRate
  sampDist / sum(sampDist)
}
updateDistCmp <- cmpfun(updateDist)


featSel <- function(dat, datCor, iter, k, lrate) {
  # minimize pairwise correlation
  # the list of selected features
  # start with two random features
  print(timestamp())

  # start out with uniform distribution
  sampDist <- rep(1, ncol(datCor)) / ncol(datCor)

  for (i in 1:iter) {

    # sample an idx
    idx <- sample(1:ncol(dat), size=k, replace=F, prob=sampDist)

    # return the corr scores
    kcors <- getCorsCmp(idx, datCor)

    # update prob dist
    sampDist <- updateDistCmp(idx, kcors, sampDist, 0.01)
  }

  return(sampDist)
}
featSelCmp <- cmpfun(featSel)


multiSel <- function(dat, datCor, reps, iter, lrate, k, cores) {

  # do feature selection `reps` number of times
  # return reps number of probability distributions
  res0 <- mclapply(1:reps, function(a) featSelCmp(dat, datCor, iter, k, lrate), mc.cores=cores)

  # sample from each, and table the results
  res1 <- lapply(1:reps, function(a) sample(1:nrow(datCor), size=k, prob=res0[[a]]))

  #count up the selections
  list(Ps=res0, Ss=res1, Qs=table(unlist(res1)))
}


solnIdxToNames <- function(setdf, soln) {
  res0  <- soln$Qs[soln$Qs %in% sort(soln$Qs, decreasing=T)[1:k]]
  idx   <- as.numeric(names(res0))
  names <- as.character(setdf$SetName[idx])
  df <- data.frame(SetNames=names, Counts=as.numeric(res0))
  df[order(df$Counts, decreasing=T),]
}


solnNamesMerge <- function(setdf, soln1, soln2) {
  res0  <- soln1$Qs
  idx0   <- as.numeric(names(res0))
  names0 <- as.character(setdf$SetName[idx0])
  df0 <- data.frame(SetNames=names0, Counts1=as.numeric(res0))

  res1  <- soln2$Qs
  idx1   <- as.numeric(names(res1))
  names1 <- as.character(setdf$SetName[idx1])
  df1 <- data.frame(SetNames=names1, Counts2=as.numeric(res1))

  dfm <- merge(df0, df1)
  dfm$Dist = sqrt(dfm$Counts1^2 + dfm$Counts2^2)
  dfm[order(dfm$Dist, decreasing=T),]
}

################################
# run it a number of times... ##
# on a subsample

# split the data into 2 parts
idx1 <- sample(1:nrow(datm), size=floor(nrow(datm)/2), replace=F)
idx2 <- (1:nrow(datm))[-idx1]

# compute the solutions
print("solns")
soln1 <- multiSel(datm[idx1,], datmcor, reps, iters, learnRate, k, cores)
soln2 <- multiSel(datm[idx2,], datmcor, reps, iters, learnRate, k, cores)

# what's the solution?
names1 <- solnIdxToNames(dat, soln1)
names2 <- solnIdxToNames(dat, soln2)

# what's the solution?
nameCounts <- solnNamesMerge(dat, soln1, soln2)
write.table(nameCounts, file="Solution_Counts.txt", sep="\t", row.names=F, quote=F)

# compare solutions
g0 <- qplot(data=nameCounts, x=Counts1, y=Counts2, geom=c("point","smooth"))
ggsave(g0, file="solution_comparison.pdf")


#print("plotting")
# plotting the number of times a feature was selected.
jdx <- as.numeric(names(soln1$Qs))
x <- data.frame(dat[jdx,c(1,2)], Scores=as.numeric(soln1$Qs))
g1 <- ggplot(data=x, aes(x=factor(SetName), y=Scores)) +##
     geom_bar(stat="identity")+
     theme(axis.text.x = element_text(angle = 50, hjust = 1))
ggsave(g1, file="selected_features_barplot_soln1.pdf")

jdx <- as.numeric(names(soln2$Qs))
x <- data.frame(dat[jdx,c(1,2)], Scores=as.numeric(soln2$Qs))
g2 <- ggplot(data=x, aes(x=factor(SetName), y=Scores)) +#
     geom_bar(stat="identity")+
     theme(axis.text.x = element_text(angle = 50, hjust = 1))
ggsave(g2, file="selected_features_barplot_soln2.pdf")


save(names1, names2, soln1, soln2, idx1, idx2, params, nameCounts, file="run_solutions_and_params.rda")
