#### RNA-SEQ TPM DATA PRE-PROCESSING AND NORMALIZATION ####

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(reshape2)

#### DATA PREPARATION ####
## IMPORT TPM DATA
setwd("~/Dropbox/RNA seq data/DysiNKT iNKT_V")
tpm.dat <- read.csv("iNKT_V TPM.csv")
head(tpm.dat)
colnames(tpm.dat)[1] <- "Gene"
boxplot(log2(tpm.dat[, -1] + 0.25))

## CONVERT ENSEMBLE GENES TO SYMBOLS AND ENTREZIDS
gene.id <- clusterProfiler::bitr(tpm.dat[, "Gene"], fromType = "ENSEMBL", 
                                 toType = c("SYMBOL", "ENTREZID"), 
                                 OrgDb = 'org.Hs.eg.db', drop = TRUE)

#### DATA PRE-PROCESSING ####
## DETECT THE DUPLICATED ENSEMBLE GENES 
gene.id_dup.ens <- data.frame("ENSEMBL" = unique(gene.id$ENSEMBL[which(duplicated(gene.id$ENSEMBL))]))

## DETECT THE DUPLICATED SYMBOL AND ENTREZID GENES 
# ACCORDINGLY, THE POSITIONS OF DUPLICATED SYMBOLS AND ENTREZIDS ARE THE SAME, THIS CODE LINE WILL BE FOCUES ON THE EITHER SYMBOLS OR ENTREZIDS
gene.id_dup.se <- data.frame("SYMBOL"= unique(gene.id$SYMBOL[which(duplicated(gene.id$SYMBOL))]))

## REMOVE THE DUPLICATED GENES
gene.id_nodup <- gene.id %>% 
  dplyr::anti_join(gene.id_dup.ens) %>%
  dplyr::anti_join(gene.id_dup.se)

## INTERSECTION OF GENES IN THE EXPRESSION MATRIX
tpm.dat_ins <- tpm.dat %>% 
  subset(Gene %in% gene.id_nodup$ENSEMBL) %>% 
  dplyr::mutate(SYMBOL = gene.id_nodup$SYMBOL, .after = "Gene") %>%
  dplyr::mutate(ENTREZID = gene.id_nodup$ENTREZID, .after = "SYMBOL")
boxplot(log2(tpm.dat_ins[, -c(1:3)] + 0.25))

#iNKT_V
iNKT_V.dat <- tpm.dat_ins[,-c(6:7)]

## FILTER GENES WITH ZERO VALUE ACROSS ALL SAMPLES OUT
iNKT_V.dat_cut.zero <- iNKT_V.dat[-which(apply(iNKT_V.dat[, -c(1:3)], 
                                                 1, sum) == 0), ]
boxplot(log2(iNKT_V.dat_cut.zero[, -c(1:3)] + 0.25))

## REMOVE GENES WITH LESS EXPRESSION VALUES
# REF: https://doi.org/10.3389/fgene.2021.632620
qt.cutoff <- 0.25 # DEFINE THE PERCENTILE CUT-OFF VALUE
num.collect <- vector() # CREATE A VECTOR FOR COLLECTING THE NUMBER OF SAMPLES WHICH HIT THE CONDITION --> EXPRESSION VALUE > QUANTILE VALUE IN EACH GENE ACROSS ALL SAMPLES
# THIS FOR LOOP WILL BE PROCESSING FOR A WHILE
for (gene.name in iNKT_V.dat_cut.zero$Gene) {
  tmp <- as.numeric(iNKT_V.dat_cut.zero[
    iNKT_V.dat_cut.zero$Gene == gene.name, -c(1:3)])
  num.collect <- append(num.collect, length(which(tmp > quantile(tmp, qt.cutoff))))
}
table(num.collect) # CHECK THE DISTRIBUTION OF SAMPLES

## A SCORE WAS CONSTRUCTED FOR EACH GENE BY COUNTING THE NUMBER OF SAMPLES WITH EXPRESSION VALUES BELOW THE 25TH GENE PERCENTILE 
samp.cutoff <- (ncol(iNKT_V.dat_cut.zero)-3)*.5 # 50% OF NUMBER OF SAMPLES
iNKT_V.dat_cut.low <- iNKT_V.dat_cut.zero[which(num.collect > samp.cutoff), ]
boxplot(log2(iNKT_V.dat_cut.low[, -c(1:3)] + 0.25))
write.csv(iNKT_V.dat_cut.low, "iNKT_V_TPM_cutoff.csv")

library(edgeR)
library(limma)
library(factoextra)
library(ggplot2)
library(FactoMineR)

#Import data
counts <- read.csv("iNKT_V_TPM_cutoff.csv", row.names = 3)
counts <- counts[,-c(1:3)]

metadata <- read.csv("iNKT_V metadata.csv",header = T)
identical(metadata$X, colnames(counts))

#Make into table format
counts <- as.data.frame(counts)

#Create DGEList object
d0 <- DGEList(counts)

#Calculate normalization factors
d0 <- calcNormFactors(d0)
d0$samples

# Initial design matrix and linear model fit
#Specifies a model where each coefficient corresponds to a group mean
mm <- model.matrix(~0 + Condition, data = metadata)
mm

#filtering
keep <- filterByExpr(d0, mm)
sum(keep) # number of genes retained

d <- d0[keep,]

#obtain log-transformed normalized expression data
logcpm <- cpm(d, log=TRUE)

#Voom
y <- voom(d, mm, plot = T)

#filter more to get a smooth voom curve
tmp <- voom(d0, mm, plot = T)

#lmFit fits a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)
head(coef(fit))

# First contrast definition and application
#Comparisons between groups
contr <- makeContrasts(ConditioniNKTV - ConditioniNKT, levels = colnames(coef(fit)))

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smoothing of standard errors of log fold changes
tmp <- eBayes(tmp)

# most differentially expressed genes
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
head(top.table)

#number of DE genes
length(which(top.table$adj.P.Val < 0.05))

#Write top.table to a file
top.table$SYMBOL <- sapply(strsplit(rownames(top.table), split = ".", fixed = TRUE), `[`, 1)

# Further refinement adjust or refine based on initial results
contrast.matrix <- makeContrasts(ConditioniNKTV - ConditioniNKT, levels=colnames(coef(fit)))
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Top500 DE genes
top.table <- topTable(fit2, n = 500)
top.table$SYMBOL <- sapply(strsplit(rownames(top.table), split = ".", fixed = TRUE), `[`, 1)


#Filter rows where logFC is greater than 2 or less than -2
filtered_data <- top.table %>%
  filter(logFC >= 2 | logFC <= -2)
write.table(filtered_data, file = "iNKTVvsiNKT.txt", row.names = F, sep = "\t", quote = F)


#Filter genes across FuniNKT and DysiNKT
DysiNKTV.dat <- tpm.dat_ins[,-c(8:9)]

## FILTER GENES WITH ZERO VALUE ACROSS ALL SAMPLES OUT
DysiNKTV.dat_cut.zero <- DysiNKTV.dat[-which(apply(DysiNKTV.dat[, -c(1:3)], 
                                                       1, sum) == 0), ]
boxplot(log2(DysiNKTV.dat_cut.zero[, -c(1:3)] + 0.25))

## REMOVE GENES WITH LESS EXPRESSION VALUES
# REF: https://doi.org/10.3389/fgene.2021.632620
qt.cutoff <- 0.25 # DEFINE THE PERCENTILE CUT-OFF VALUE
num.collect <- vector() # CREATE A VECTOR FOR COLLECTING THE NUMBER OF SAMPLES WHICH HIT THE CONDITION --> EXPRESSION VALUE > QUANTILE VALUE IN EACH GENE ACROSS ALL SAMPLES
# THIS FOR LOOP WILL BE PROCESSING FOR A WHILE
for (gene.name in DysiNKTV.dat_cut.zero$Gene) {
  tmp <- as.numeric(DysiNKTV.dat_cut.zero[
    DysiNKTV.dat_cut.zero$Gene == gene.name, -c(1:3)])
  num.collect <- append(num.collect, length(which(tmp > quantile(tmp, qt.cutoff))))
}
table(num.collect) # CHECK THE DISTRIBUTION OF SAMPLES

## A SCORE WAS CONSTRUCTED FOR EACH GENE BY COUNTING THE NUMBER OF SAMPLES WITH EXPRESSION VALUES BELOW THE 25TH GENE PERCENTILE 
samp.cutoff <- (ncol(DysiNKTV.dat_cut.zero)-3)*.5 # 50% OF NUMBER OF SAMPLES
DysiNKTV.dat_cut.low <- DysiNKTV.dat_cut.zero[which(num.collect > samp.cutoff), ]
boxplot(log2(DysiNKTV.dat_cut.low[, -c(1:3)] + 0.25))
write.csv(DysiNKTV.dat_cut.low, "DysiNKTV_TPM_cutoff.csv")


library(dplyr)
#load tpm data
tpm.data <- read.csv("DysiNKTV_TPM_cutoff.csv")
# Merge the data frames based on SYMBOL
merged_data <- merge(filtered_data, tpm.data, by = "SYMBOL")

# View the merged data
print(merged_data)
merged_data <- merged_data[,-c(2:8)]
write.csv(merged_data, file = "iNKTVvsiNKT TPM.csv")

################

#DysiNKT
DysiNKT.dat <- tpm.dat_ins[,-c(4:5)]

## FILTER GENES WITH ZERO VALUE ACROSS ALL SAMPLES OUT
DysiNKT.dat_cut.zero <- DysiNKT.dat[-which(apply(DysiNKT.dat[, -c(1:3)], 
                                                 1, sum) == 0), ]
boxplot(log2(DysiNKT.dat_cut.zero[, -c(1:3)] + 0.25))

## REMOVE GENES WITH LESS EXPRESSION VALUES
# REF: https://doi.org/10.3389/fgene.2021.632620
qt.cutoff <- 0.25 # DEFINE THE PERCENTILE CUT-OFF VALUE
num.collect <- vector() # CREATE A VECTOR FOR COLLECTING THE NUMBER OF SAMPLES WHICH HIT THE CONDITION --> EXPRESSION VALUE > QUANTILE VALUE IN EACH GENE ACROSS ALL SAMPLES
# THIS FOR LOOP WILL BE PROCESSING FOR A WHILE
for (gene.name in DysiNKT.dat_cut.zero$Gene) {
  tmp <- as.numeric(DysiNKT.dat_cut.zero[
    DysiNKT.dat_cut.zero$Gene == gene.name, -c(1:3)])
  num.collect <- append(num.collect, length(which(tmp > quantile(tmp, qt.cutoff))))
}
table(num.collect) # CHECK THE DISTRIBUTION OF SAMPLES

## A SCORE WAS CONSTRUCTED FOR EACH GENE BY COUNTING THE NUMBER OF SAMPLES WITH EXPRESSION VALUES BELOW THE 25TH GENE PERCENTILE 
samp.cutoff <- (ncol(DysiNKT.dat_cut.zero)-3)*.5 # 50% OF NUMBER OF SAMPLES
DysiNKT.dat_cut.low <- DysiNKT.dat_cut.zero[which(num.collect > samp.cutoff), ]
boxplot(log2(DysiNKT.dat_cut.low[, -c(1:3)] + 0.25))
write.csv(DysiNKT.dat_cut.low, "DysiNKT_TPM_cutoff.csv")

#Import data
counts <- read.csv("DysiNKT_TPM_cutoff.csv", row.names = 3)
counts <- counts[,-c(1:3)]

metadata <- read.csv("DysiNKT meta.csv",header = T)
identical(metadata$X, colnames(counts))

#Make into table format
counts <- as.data.frame(counts)

#Create DGEList object
d0 <- DGEList(counts)

#Calculate normalization factors
d0 <- calcNormFactors(d0)
d0$samples

# Initial design matrix and linear model fit
#Specifies a model where each coefficient corresponds to a group mean
mm <- model.matrix(~0 + Condition, data = metadata)
mm

#filtering
keep <- filterByExpr(d0, mm)
sum(keep) # number of genes retained

d <- d0[keep,]

#obtain log-transformed normalized expression data
logcpm <- cpm(d, log=TRUE)

#Voom
y <- voom(d, mm, plot = T)

#filter more to get a smooth voom curve
tmp <- voom(d0, mm, plot = T)

#lmFit fits a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)
head(coef(fit))

# First contrast definition and application
#Comparisons between groups
contr <- makeContrasts(ConditionDysiNKT - ConditioniNKT, levels = colnames(coef(fit)))

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

#smoothing of standard errors of log fold changes
tmp <- eBayes(tmp)

# most differentially expressed genes
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
head(top.table)

#number of DE genes
length(which(top.table$adj.P.Val < 0.05))

#Write top.table to a file
top.table$SYMBOL <- sapply(strsplit(rownames(top.table), split = ".", fixed = TRUE), `[`, 1)

# Further refinement adjust or refine based on initial results
contrast.matrix <- makeContrasts(ConditionDysiNKT - ConditioniNKT, levels=colnames(coef(fit)))
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Top500 DE genes
top.table <- topTable(fit2, n = 500)
top.table$SYMBOL <- sapply(strsplit(rownames(top.table), split = ".", fixed = TRUE), `[`, 1)


#Filter rows where logFC is greater than 2 or less than -2
filtered_data <- top.table %>%
  filter(logFC >= 2 | logFC <= -2)
write.table(filtered_data, file = "DysiNKTvsiNKT.txt", row.names = F, sep = "\t", quote = F)

library(dplyr)
#load tpm data
tpm.data <- read.csv("DysiNKTV_TPM_cutoff.csv")
filtered_data <- read.table("DysiNKTvsiNKT.txt", sep = "\t", header = T)
# Merge the data frames based on SYMBOL
merged_data <- merge(filtered_data, tpm.data, by = "SYMBOL")

# View the merged data
print(merged_data)
merged_data <- merged_data[,-c(2:8)]
write.csv(merged_data, file = "DysiNKTvsiNKT TPM.csv")

##merge FuniNKTvsiNKT TPM and DysiNKTvsiNKT TPM
##Run PLSR in JMP software
##calculate VIP score
##VIP>0.8 genes >>> Enrichr (MSigDB hallmark pathway)

#Fig. 6c
#load hallmarks pathway data
Hallmarks <- read.csv("DysiNKTV Hallmarks.csv", header = TRUE, 
                      stringsAsFactors = FALSE)
#arrange descending order of Pvalue
df <- Hallmarks %>%
  arrange(desc("P.value"))

# Reorder the Pathways factor according to the sorted FDR values
df$Term <- factor(df$Term, levels = df$Term[order(-df$P.value)])

ggplot(df, aes(x= Combined.Score, y=Term, size= Overlap, color= P.value)) + geom_point() +
  theme_classic() +
  scale_color_gradient(low = "red", high = "black", breaks = c(0.01,0.05, 0.1))+ 
  labs(x = "Combined Score",y = "MSigDB_Hallmark 2020",size = "Overlap score",color = "P-value")


#Merge PLSR data with TPM
DGE.tpm <- read.csv("Top500 DysiNKTV PLSR VIP.csv")
#load tpm data
tpm.data <- read.csv("DysiNKTV_TPM_cutoff.csv")
# Merge the data frames based on SYMBOL
merged_data <- merge(DGE.tpm, tpm.data, by = "SYMBOL")
# View the merged data
print(merged_data)
merged_data <- merged_data[,-c(4:6)]
write.csv(merged_data, file = "DysiNKTV PLSR VIP TPM.csv")


#### DATA NORMALIZATION ####
# REF: https://github.com/sisyspharm/ccaexplorer/tree/master 
# @ LINE 102
## CONVERT TPM DATA TO MEDIAN-CENTERED LOG2-TPM DATA
# CALCULATE MEDIAN VALUE OF EACH GENE ACROSS ALL SAMPLES

VIP.tpm <- read.csv("DysiNKTV PLSR VIP TPM.csv", row.names = 2)
tpm.median  <- apply(log2(VIP.tpm[, -c(1:3)] + 0.25), 1, median)
# USE MEDIAN FOR SUBTRACTING THE LOG2-TPM VALUE
tpm.dat_median <- log2(VIP.tpm[, -c(1:3)] + 0.25) - tpm.median
# CONSTRUCT A NEW DATA FRAME WITH GENES AND EXPRESSIONS
tpm.dat_median <- cbind(VIP.tpm[,c(2:3)], tpm.dat_median)
boxplot(tpm.dat_median[, -c(1:3)])
write.csv(tpm.dat_median, "iNKTVvsDysiNKT_normalized.csv")

#Plot Fig. 2c heatmap in graphpad prism using FuniNKTvsDysiNKT_normalized, VIP score, and gene involved in each pathway
#Vendiagram plot using venny 2.1.0
