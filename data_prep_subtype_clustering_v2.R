
library(ggplot2)
library(stringr)

# First get all the tables read in #

barcodes <- read.table("data/barcodes/EBpp_whitelist_no_blood_cancers.csv", sep=",", header=T, stringsAsFactors=F)
aliquotBarcodes <- unique(barcodes$AliquotBarcode)
aliquotBarcodes <- str_replace_all(string=aliquotBarcodes, pattern='-', replacement='\\.')
length(aliquotBarcodes)
#[1] 9163
participantBarcodes <- str_sub(aliquotBarcodes, start=1, end=12)

# COMBINED ssGSEA scores
#res0 <- read.table("data/Signatures/ssGSEA_Scores_Wolf68_Bindea_Yasin_C7.tsv.gz", sep="\t", header=T, stringsAsFactors=T)
load("data/Signatures/ssGSEA_Scores_Wolf68_Bindea_Yasin_C7.rda") # faster

# ssGSEA Wolf SCORES
ssWolf <- res0[res0$Source == "Wolf",]
ssWolfNames <- ssWolf$SetName
rownames(ssWolf) <- ssWolfNames

# CELL TYPE SCORES
bindea <- res0[res0$Source == "Bindea",]
bindeaNames <- bindea$SetName
bindeaNames <- as.character(bindea$SetName)

## MSKCC
yasin <- res0[res0$SetName %in% c('Angiogenesis', 'APM1'),]
yasinNames <- as.character(yasin$SetName)
# corrections
yasinNames[2] <- "Antigen presenting machinery"


## Wolf Scores
wolf <- read.table("data/Signatures/TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.csv", sep=",", header=T, stringsAsFactors=F)
wolfNames <- wolf$X
rownames(wolf) <- wolfNames
wolf <- wolf[, colnames(wolf) %in% aliquotBarcodes]
ssWolf <- ssWolf[, colnames(wolf)]
all(colnames(ssWolf) == colnames(wolf))
#TRUE
# 9129

## Attractors
anastRNA <- read.table("data/Signatures/mRNA_immune_attractors_EBpp_whitelist.tsv", sep="\t", header=T, stringsAsFactors=F)
anastNames <- rownames(anastRNA)

## Cibersort scores cell types
ciber <- read.table("data/Signatures/TCGA.Kallisto.cibersort.relative.tsv", sep="\t", header=T, stringsAsFactors=F)
ciberNames <- colnames(ciber)[-c(1,2,23,24,25,26)]
ciberSub <- ciber[,-c(1,2,23,24,25,26)]
ciberSub <- t(ciberSub)
colnames(ciberSub) <- ciber$SampleID


# The dictionary atoms separate better
c7atoms <- read.table("data/Signatures/c7_dict_transpose/modl_atoms.txt", sep="\t", header=F)
rownames(c7atoms) <- sapply(1:32, function(a) paste("C7_Atom_",a,sep=""))
c7atomsNames <- rownames(c7atoms)
# BUG:: colnames(c7atoms) <- aliquotBarcodes
colnames(c7atoms) <- colnames(res0)[-c(1,2)]


## ICR
icr <- read.table("data/Signatures/ICR_Scores.csv", sep=",", header=T, stringsAsFactors=F)
icr_ids <- icr$Sample_ID
icr_ids_dots <- str_replace_all(string=icr_ids, pattern='-', replacement='\\.')
icr_sub <- icr[icr_ids_dots %in% aliquotBarcodes,]
icr_sub_t <- t(icr_sub[,4:6])
colnames(icr_sub_t) <- icr_ids_dots[icr_ids_dots %in% aliquotBarcodes]
icr_sub_t <- icr_sub_t[, aliquotBarcodes]
icrNames <- rownames(icr_sub_t)

# Prep the barcodes
shortBarcodes <- colnames(wolf)
partiBarcodes <- str_sub(shortBarcodes, start=1, end=12)


#set 1
#Shared Starting Point 1
#4 Set origins
#	• 2 Yasin et al signatures, Angiogenesis and Antigen Presenting Machinery
#	• 68 Wolf et al signatures
#	• 9 Anastassiou
#	• 3 ICR

a <- yasin[,shortBarcodes]
b <- wolf[,shortBarcodes]
c <- anastRNA[,shortBarcodes]
d <- icr_sub_t[,shortBarcodes]
all(colnames(a) == colnames(b))
#[1] TRUE
all(colnames(a) == colnames(c))
#[1] TRUE
all(colnames(a) == colnames(d))
#[1] TRUE

set1 <- rbind(a,b,c,d)
rownames(set1) <- set1$SetName
source1 <- c(rep("Yasin", nrow(a)), rep("Wolf", nrow(b)), rep("Attractors", nrow(c)), rep("ICR", nrow(d)))
names1 <- c(yasinNames, wolfNames, anastNames, icrNames)

set1_2 <- cbind(data.frame(Source=source1, SetName=names1), set1)
write.table(set1_2, file="Signatures_Set1_2.tsv", sep="\t", row.names=F, quote=F)

#set 2
#Shared Starting Point 2
#5 Set origins
#	• 2 Yasin et al signatures, Angiogenesis and Antigen Presenting Machinery
#	• 68 Wolf et al signatures
#	• 9 Anastassiou
#	• 3 ICR
#	• ImmuneSigDB C7

a <- yasin[,shortBarcodes]
b <- wolf[,shortBarcodes]
c <- anastRNA[,shortBarcodes]
d <- icr_sub_t[,shortBarcodes]
e <- c7atoms[,shortBarcodes]
all(colnames(a) == colnames(b))
#[1] TRUE
all(colnames(a) == colnames(c))
#[1] TRUE
all(colnames(a) == colnames(d))
#[1] TRUE
all(colnames(a) == colnames(e))

set2 <- rbind(a,b,c,d,e)
rownames(set2) <- set2$SetName
source2 <- c(rep("Yasin", nrow(a)), rep("Wolf", nrow(b)), rep("Attractors", nrow(c)),
        rep("ICR", nrow(d)), rep("c7atoms", nrow(e)))
names2 <- c(yasinNames, wolfNames, anastNames, icrNames, c7atomsNames)

set2_2 <- cbind(data.frame(Source=source2, SetName=names2), set2)
write.table(set2_2, file="Signatures_Set2_2.tsv", sep="\t", row.names=F, quote=F)

# set 3
#Shared Starting Point 3
#7 Set origins
#	• 24 Bindea cell signature
#	• CIBERSORT
#	• 2 Yasin et al signatures, Angiogenesis and Antigen Presenting Machinery
#	• 68 Wolf et al signatures
#	• 9 Anastassiou
#	• 3 ICR
#	• ImmuneSigDB C7

a <- yasin[,shortBarcodes]
b <- wolf[,shortBarcodes]
c <- anastRNA[,shortBarcodes]
d <- icr_sub_t[,shortBarcodes]
e <- c7atoms[,shortBarcodes]
f <- bindea[,shortBarcodes]
g <- ciberSub[,partiBarcodes]

# corrections
colnames(g) <- colnames(f)
g <- g[1:20,]
ciberNames <- ciberNames[1:20]

all(colnames(a) == colnames(b))
#[1] TRUE
all(colnames(a) == colnames(c))
#[1] TRUE
all(colnames(a) == colnames(d))
#[1] TRUE
all(colnames(a) == colnames(e))
#[1] TRUE
all(colnames(a) == colnames(f))
#[1] TRUE

set3 <- rbind(a,b,c,d,e,f,g)
rownames(set3) <- set3$SetName
source3 <- c(rep("Yasin", nrow(a)), rep("Wolf", nrow(b)), rep("Attractors", nrow(c)),
        rep("ICR", nrow(d)), rep("c7atoms", nrow(e)), rep("Bindea",nrow(f)),
        rep("Cibersort", nrow(g)))
names3 <- c(yasinNames, wolfNames, anastNames, icrNames, c7atomsNames, bindeaNames, ciberNames)

set3_2 <- cbind(data.frame(Source=source3, SetName=names3), set3)
write.table(set3_2, file="Signatures_Set3_2.tsv", sep="\t", row.names=F, quote=F)


save(set1_2, set2_2, set3_2, file="Signatures_for_Subtyping_2.rda")
save(yasin, wolf, anastRNA, icr_sub_t, c7atoms, bindea, ciberSub, aliquotBarcodes, shortBarcodes, file="Signature_Components.rda")
