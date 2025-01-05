setwd("D:/TEP/Ovarian Cancer/Validation")
# Compare with ImPlatelet ----
load("dgeGenesEnsembl75.rdata") #READ genes.csv file instead

ImPlatelet<-read.delim("ImPlatelet_samples.tsv",sep="\t")

library(bladderbatch)
library(magrittr)
library(dplyr)

ImPlatelet_OC=ImPlatelet%>%
  filter(OriginalGroup=="OC")

ImPlatelet_OC <- ImPlatelet_OC[-((nrow(ImPlatelet_OC) - 4):nrow(ImPlatelet_OC)), ]

#remove missing stage

ImPlatelet_HC=ImPlatelet%>%
  filter(OriginalGroup=="HC")

ImPlatelet_HC=ImPlatelet%>%
  filter(Gender=="2")

ImPlatelet_samples=rbind(ImPlatelet_HC, ImPlatelet_OC)

# Exclude Y-chr genes ----

ImPlatelet_counts<-read.delim("ImPlatelet_counts_raw.tsv",sep="\t")

Y_chr=filter(genes, genes$chromosome_name=='Y')
Y_rows=rownames(Y_chr)

Y_genes <- ImPlatelet_counts[rownames(ImPlatelet_counts) %in% Y_rows, ]
Y_genes =rownames(Y_genes)


ImPlatelet_samples$EarlyStage[is.na(ImPlatelet_samples$EarlyStage)] <- "Control"

ImPlatelet_samples$Stage[is.na(ImPlatelet_samples$Stage)] <- "Control"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IC", "IA")] <- "Stage1"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IIA", "IIB", "IIC")] <- "Stage2"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IIIA", "IIIB", "IIIC")] <- "Stage3"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IVB")] <- "Stage4"

ImPlatelet_samples=ImPlatelet_samples[, -c(8,9,10)]

selected_samples=ImPlatelet_samples$Id

ImPlatelet_counts <- ImPlatelet_counts[, selected_samples]

table(ImPlatelet_samples$Stage)

###Create DGE object ####
library(edgeR)

dge <- DGEList(counts = ImPlatelet_counts,
               group = ImPlatelet_samples$Stage,
               genes = genes[which(rownames(genes) %in% rownames(ImPlatelet_counts)),]
)           


ImPlatelet_samples=ImPlatelet_samples[,-c(2,3,4,9,6,7)]

dge$samples <- cbind(dge$samples, ImPlatelet_samples)

dge$samples$Stage<- factor(dge$samples$Stage, levels = c(levels(dge$samples$Stage),"Control","Stage1","Stage2","Stage3","Stage4"))

dge$samples$Stage[which(dge$samples$Stage %in% c("Control"))] <- "Control"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage1"))] <- "Stage1"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage2"))] <- "Stage2"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage3"))] <- "Stage3"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage4"))] <- "Stage4"

dim(dge)

summary(dge$samples$lib.size)
summary(dge$samples$group)
summary(dge$samples$Stage)

# Filter lowly expressed genes ----

keep.exprs <- filterByExpr(dge, group=dge$samples$Stage)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)
dge=calcNormFactors(dge)


library('EDASeq')
library(Biobase)
library(scater)

#RLE plot before RUVg correction ----

sce1 <- SingleCellExperiment(assays = list(counts =dge$counts))

plotRLE(sce1, exprs_values = "counts", exprs_logged=FALSE, style = "minimal")+
  ggtitle("Before RUVg correction")

# Make RUVg correction ----

library(scater)
library(EDASeq)

ord_count=dge$counts

library(ggplot2)
library(tidyr)

library(RUVSeq)
library(foreach)

filtered_cts <- ord_count

factor=as.factor(dge$samples$Stage)

set1 <- newSeqExpressionSet(as.matrix(ord_count),
                            phenoData = data.frame(factor, 
                                                   row.names=colnames(filtered_cts)))

design = model.matrix(~0+factor)
design=as.matrix(design)
design <- design[, c(2,3,4,5,1)]
y <- DGEList(counts=counts(set1), group=factor)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=ncol(fit$design))
top <- topTags(lrt, n=nrow(set1))$table
empirical <- rownames(set1)[which(!(rownames(set1) %in% rownames(top)[1:5000]))]

library(RUVSeq)

set2 <- RUVg(set1, empirical, k=3, isLog = FALSE)

W=pData(set2)

ruv_counts=as.data.frame(set2@assayData[["normalizedCounts"]])

rownames(W)=NULL

W=W[,-1]

W=as.matrix(W)

#https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html
#https://support.bioconductor.org/p/67791/#9155690

# Load required libraries
library(SingleCellExperiment)
library(scater)
library(ggplot2)

sce2 <- SingleCellExperiment(assays = list(counts = ruv_counts))

colData(sce2)$Stage <- dge$samples$Stage

plotRLE(sce2, exprs_values = "counts", exprs_logged = FALSE, style = "minimal") +
  ggtitle("After RUVg correction")


design <- model.matrix(~0+dge$samples$Stage)

colnames(design)[1] ="Control"
colnames(design)[2] ="Stage1"
colnames(design)[3] ="Stage2"
colnames(design)[4] ="Stage3"
colnames(design)[5] ="Stage4"

design <- design[, c(2,3,4,5,1)]

IMP_lcpm=cpm.DGEList(dge, log = TRUE)

dge$counts=IMP_lcpm

IMP_lcpm=removeBatchEffect(dge, design=design, covariates = as.matrix(W))

IMP_lcpm=as.data.frame(IMP_lcpm)

val_counts=IMP_lcpm
write.csv(val_counts,'external_validation_counts.csv')

sessionInfo()