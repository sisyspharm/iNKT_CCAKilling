#### RNA-SEQ TPM DATA PRE-PROCESSING AND NORMALIZATION ####

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(reshape2)

#### DATA PREPARATION ####
## IMPORT TPM DATA
setwd("~/Dropbox/RNA seq data/9_CCA")
tpm.dat <- read.csv("9_CCA_TPM.csv")
head(tpm.dat)

#### DATA PRE-PROCESSING ####
## DETECT THE DUPLICATED ENSEMBLE GENES 
tpm.dat_nodup <- tpm.dat[!duplicated(tpm.dat$SYMBOL), ]

## FILTER GENES WITH ZERO VALUE ACROSS ALL SAMPLES OUT
tpm.dat_cut.zero <- tpm.dat_nodup[-which(apply(tpm.dat_nodup[, -c(1)], 
                                               1, sum) == 0), ]

## REMOVE GENES WITH LESS EXPRESSION VALUES
# REF: https://doi.org/10.3389/fgene.2021.632620
qt.cutoff <- 0.25 # DEFINE THE PERCENTILE CUT-OFF VALUE
num.collect <- vector() # CREATE A VECTOR FOR COLLECTING THE NUMBER OF SAMPLES WHICH HIT THE CONDITION --> EXPRESSION VALUE > QUANTILE VALUE IN EACH GENE ACROSS ALL SAMPLES
# THIS FOR LOOP WILL BE PROCESSING FOR A WHILE
for (gene.name in tpm.dat_cut.zero$SYMBOL) {
  tmp <- as.numeric(tpm.dat_cut.zero[
    tpm.dat_cut.zero$SYMBOL == gene.name, -c(1)])
  num.collect <- append(num.collect, length(which(tmp > quantile(tmp, qt.cutoff))))
}
table(num.collect) # CHECK THE DISTRIBUTION OF SAMPLES

## A SCORE WAS CONSTRUCTED FOR EACH GENE BY COUNTING THE NUMBER OF SAMPLES WITH EXPRESSION VALUES BELOW THE 25TH GENE PERCENTILE 
samp.cutoff <- (ncol(tpm.dat_cut.zero)-3)*.5 # 50% OF NUMBER OF SAMPLES
tpm.dat_cut.low <- tpm.dat_cut.zero[which(num.collect > samp.cutoff), ]
boxplot(log2(tpm.dat_cut.low[, -c(1)] + 0.25))
write.csv(tpm.dat_cut.low, "9_CCA_TPM_cutoff.csv")

library(edgeR)
library(limma)
library(factoextra)
library(ggplot2)
library(FactoMineR)
library(dplyr)

#Import data
counts <- read.csv("9_CCA_TPM_cutoff.csv", row.names = 2)
counts <- counts[,-1]

metadata <- read.csv("9_CCA_meta.csv",header = T)
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
contr <- makeContrasts(ConditioniNKTNR - ConditioniNKTR, levels = colnames(coef(fit)))

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
contrast.matrix <- makeContrasts(ConditioniNKTNR - ConditioniNKTR, levels=colnames(coef(fit)))
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Top500 DE genes
top.table <- topTable(fit2, n = 500)
top.table$SYMBOL <- sapply(strsplit(rownames(top.table), split = ".", fixed = TRUE), `[`, 1)
write.table(top.table, file = "9_CCA_iNKTNRvsiNKTR.txt", row.names = F, sep = "\t", quote = F)

#Filter rows where logFC is greater than 2 or less than -2
filtered_data <- top.table %>%
  filter(logFC >= 2 | logFC <= -2)
write.table(filtered_data, file = "9_CCA_iNKTNRvsiNKTR.txt", row.names = F, sep = "\t", quote = F)

library(dplyr)
#load tpm data
tpm.data <- read.csv("9_CCA_TPM_cutoff.csv")
# Merge the data frames based on SYMBOL
merged_data <- merge(filtered_data, tpm.data, by = "SYMBOL")

# View the merged data
print(merged_data)
merged_data <- merged_data[,-c(2:8)]
write.csv(merged_data, file = "9_CCA_iNKTNRvsiNKTR TPM.csv")

## using "9_CCA_iNKTNRvsiNKTR TPM.csv" data, we run PLSR in JMP software
## calculate VIP score
## VIP >0.8 genes >>> Enrichr (KEGG pathway analysis)

#load KEGG pathway data
KEGG <- read.csv("9_CCA KEGG.csv", header = TRUE, 
                 stringsAsFactors = FALSE)
#arrange descending order of Pvalue
df <- KEGG %>%
  arrange(desc("P-value"))

# Reorder the Pathways factor according to the sorted FDR values
df$Pathway <- factor(df$Pathway, levels = df$Pathway[order(-df$P.value)])

ggplot(df, aes(x= Combined.Score, y=Pathway, size= Overlap, color= P.value)) + geom_point() +
  theme_classic() +
  scale_color_gradient(low = "red", high = "black", breaks = c(0.01,0.03, 0.05, 0.07))+ 
  labs(x = "Combined Score",y = "KEGG 2021 Human",size = "Overlap score",color = "P-value")

## gene names and cutoff line added manually