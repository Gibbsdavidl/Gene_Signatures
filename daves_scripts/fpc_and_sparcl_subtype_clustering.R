
library(ggplot2)
library(stringr)
library(pheatmap)
library(colorRamps)

load("Signatures_for_Subtyping.rda")

set1 <- t(set1_1[,-c(1:2)])
colnames(set1) <- set1_1$SetName

set2 <- t(set2_1[,-c(1,2)])
colnames(set2) <- set2_1$SetName

set1cor <- cor(set1, method="spearman")
pheatmap(mat=set1cor, fontsize_col=3, fontsize_row=3, border_color=NA, color=blue2green2red(400))

set2cor <- cor(set2, method="spearman")
pheatmap(mat=set2cor, fontsize_col=3, fontsize_row=3, border_color=NA, color=blue2green2red(400))


# sparse clustering
library(sparcl)
# choose tuning parameter
resp <- KMeansSparseCluster.permute(x=set1, K=3, wbounds=seq(3,7,len=15), nperms=10)
res0 <- KMeansSparseCluster(x, K=3, wbounds=km.perm$bestw)


# fpc
# predictive strength
library(fpc)
res0 <- prediction.strength(set1, Gmin=2, Gmax=10, M=50)
