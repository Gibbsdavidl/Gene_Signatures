
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
res0 <- read.table("data/Signatures/ssGSEA_Scores_Wolf68_Bindea_Yasin_C7.tsv.gz", sep="\t", header=T, stringsAsFactors=T)

# ssGSEA Wolf SCORES
ssWolf <- res0[res0$Source == "Wolf",]
ssWolfNames <- ssWolf$SetName
rownames(ssWolf) <- ssWolfNames

# CELL TYPE SCORES
bindea <- res0[res0$Source == "Bindea",]
bindeaNames <- bindea$SetName
bindeaNames <- as.character(bindea$SetName)

## MSKCC
yasin <- res0[res0$SetName %in% c('Angiogenesis', 'APM1', 'APM2'),]
yasinNames <- as.character(yasin$SetName)

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

## Dictionary Learning of the c7 scores
c7codes <- read.table("data/Signatures/c7_dict_learning/modl_sparse_codes.txt", sep="\t", header=F)
colnames(c7codes) <- sapply(1:32, function(a) paste("C7_Code_",a,sep=""))
c7codesNames <- colnames(c7codes)
c7codes <- t(c7codes)
colnames(c7codes) <- aliquotBarcodes

# The dictionary atoms separate better
c7atoms <- read.table("data/Signatures/c7_dict_transpose/modl_atoms.txt", sep="\t", header=F)
rownames(c7atoms) <- sapply(1:32, function(a) paste("C7_Atom_",a,sep=""))
c7atomsNames <- rownames(c7atoms)
colnames(c7atoms) <- aliquotBarcodes


## ICR
icr <- read.table("data/Signatures/ICR_Scores.csv", sep=",", header=T, stringsAsFactors=F)
icr_ids <- icr$Sample_ID
icr_ids_dots <- str_replace_all(string=icr_ids, pattern='-', replacement='\\.')
icr_sub <- icr[icr_ids_dots %in% aliquotBarcodes,]
icr_sub_t <- t(icr_sub[,4:6])
colnames(icr_sub_t) <- icr_ids_dots[icr_ids_dots %in% aliquotBarcodes]
icr_sub_t <- icr_sub_t[, aliquotBarcodes]
icrNames <- rownames(icr_sub_t)

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
source1 <- c(rep("Yasin", nrow(a)), rep("Wolf", nrow(b)), rep("Attractors", nrow(c)), rep("ICR", nrow(d)))
names1 <- c(yasinNames, wolfNames, anastNames, icrNames)

set1_1 <- cbind(data.frame(Source=source1, SetName=names1), set1)
write.table(set1_1, file="Signatures_Set1.tsv", sep="\t", row.names=F, quote=F)

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
source2 <- c(rep("Yasin", nrow(a)), rep("Wolf", nrow(b)), rep("Attractors", nrow(c)),
        rep("ICR", nrow(d)), rep("c7atoms", nrow(c7codes)))
names2 <- c(yasinNames, wolfNames, anastNames, icrNames, c7atomsNames)

set2_1 <- cbind(data.frame(Source=source2, SetName=names2), set2)
write.table(set2_1, file="Signatures_Set2.tsv", sep="\t", row.names=F, quote=F)

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
colnames(g) <- colnames(f)

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
source3 <- c(rep("Yasin", nrow(a)), rep("Wolf", nrow(b)), rep("Attractors", nrow(c)),
        rep("ICR", nrow(d)), rep("c7atoms", nrow(c7codes)), rep("Bindea",nrow(bindea)),
        rep("Cibersort", nrow(ciberSub)))
names3 <- c(yasinNames, wolfNames, anastNames, icrNames, c7atomsNames, bindeaNames, ciberNames)

set3_1 <- cbind(data.frame(Source=source3, SetName=names3), set3)
write.table(set3_1, file="Signatures_Set3.tsv", sep="\t", row.names=F, quote=F)
