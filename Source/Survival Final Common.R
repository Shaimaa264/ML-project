
#https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html#Transcriptomic_analysis

library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library(magrittr)
library(dplyr)
library(dbplyr)
library(tidygraph)
library(tidyr)
library(tidyselect)
library(tidytree)
library(tidyverse)
library(Trendy)
library(segmented)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TCGAWorkflowData)
library(DT)

clinical_ov <- GDCquery_clinic("TCGA-OV")

samples=clinical_ov$submitter_id

query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = samples,
  sample.type = "Primary Tumor"
)

GDCdownload(query = query)

data <- GDCprepare(query = query)

output_ov <- getResults(query)

dim(data)

colnames(colData(data))

saveRDS(object = data,
        file = "tcga_data.RDS",
        compress = FALSE)

###### LOAD DATA AGAIN ######
setwd("D:/TEP/Ovarian Cancer/Survival")

data = readRDS(file = "tcga_data.RDS")

clinical_data = as.data.frame(colData(data))
batch=get_IDs(data)
gene_info=rowRanges(data)
gene_info=gene_info@elementMetadata@listData

###### RNA Seq Normalization ######

exp_ov_preprocessed <- TCGAanalyze_Preprocessing(
  object = data,
  cor.cut = 0.6,    
  datatype = "unstranded",
  filename = "OV_IlluminaHiSeq_RNASeqV2.png"
)

exp_ov_normalized <- TCGAanalyze_Normalization(
  tabDF = exp_ov_preprocessed,
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "gcContent"
)

exp_ov_filtered <- TCGAanalyze_Filtering(
  tabDF = exp_ov_normalized,
  method = "quantile",
  qnt.cut =  0.25
) 


counts=as.data.frame(exp_ov_filtered)

rm("exp_ov_normalized")
rm("exp_ov_preprocessed")
rm(data)

table(rownames(clinical_data)==colnames(counts))

# Ensure the rownames of clinical_data match the colnames of counts
if (all(rownames(clinical_data) == colnames(counts))) {
  message("The rownames of clinical_data match the colnames of counts.")
} else {
  # Reorder the columns of counts to match the rownames of clinical_data
  counts <- counts[, rownames(clinical_data)]
  
  # Check if the rearrangement was successful
  if (all(rownames(clinical_data) == colnames(counts))) {
    message("Columns of counts have been successfully rearranged.")
  } else {
    message("There was an issue rearranging the columns of counts.")
  }
}


clinical_data <- clinical_data[clinical_data$vital_status != "Not Reported", ]


library(dtplyr)
library(dbplyr)
library(dplyr)


clinical_data <- clinical_data %>%
  left_join(clinical_ov %>% select(submitter_id, treatments_pharmaceutical_treatment_or_therapy),
            by = c("patient" = "submitter_id"))

clinical_data <- clinical_data %>%
  left_join(clinical_ov %>% select(submitter_id, treatments_radiation_treatment_or_therapy),
            by = c("patient" = "submitter_id"))


clinical_data$deceased = clinical_data$vital_status == "Dead"

clinical_data$overall_survival = ifelse(clinical_data$deceased,
                                        clinical_data$days_to_death,
                                        clinical_data$days_to_last_follow_up)

table(clinical_data$figo_stage)

clinical_data$figo_stage = gsub("Stage IC", "1", clinical_data$figo_stage)
clinical_data$figo_stage <- gsub("Stage IIA|Stage IIB|Stage IIC", "2", clinical_data$figo_stage)
clinical_data$figo_stage <- gsub("Stage IIIA|Stage IIIB|Stage IIIC", "3", clinical_data$figo_stage)
clinical_data$figo_stage <- gsub("Stage IV", "4", clinical_data$figo_stage)

table(clinical_data$figo_stage)
sum(is.na(clinical_data$figo_stage))
#remove samples with NA stage 
clinical_data <- clinical_data[!is.na(clinical_data$figo_stage), ]
table(clinical_data$figo_stage)
sum(is.na(clinical_data$figo_stage))
sum(is.na(clinical_data$overall_survival))
na_indices <- which(is.na(clinical_data$overall_survival))
clinical_data <- clinical_data[!is.na(clinical_data$overall_survival), ]
table(clinical_data$figo_stage)

clinical_data$figo_stage = gsub("1", "early", clinical_data$figo_stage)
clinical_data$figo_stage <- gsub("2", "early", clinical_data$figo_stage)
clinical_data$figo_stage <- gsub("3", "late", clinical_data$figo_stage)
clinical_data$figo_stage <- gsub("4", "late", clinical_data$figo_stage)


# Calculate the mean of age_at_diagnosis, excluding NA values
mean_age <- mean(clinical_data$age_at_diagnosis, na.rm = TRUE)
clinical_data <- clinical_data %>%
  mutate(age_at_diagnosis = ifelse(is.na(age_at_diagnosis), mean_age, age_at_diagnosis))


#### Attach wet lab information #####

clinical_data <- clinical_data %>%
  left_join(batch %>% select(patient, tss),
            by = c("patient"))

clinical_data <- clinical_data %>%
  left_join(batch %>% select(patient, center),
            by = c("patient"))

clinical_data <- clinical_data %>%
  left_join(batch %>% select(patient, portion),
            by = c("patient"))

clinical_data <- clinical_data %>%
  left_join(batch %>% select(patient, plate),
            by = c("patient"))


rm(batch, clinical_ov, counts)


samples=clinical_data$barcode
exp_ov_filtered=exp_ov_filtered[,samples]
clinical_data$figo_stage <- as.factor(clinical_data$figo_stage)

#### Process all samples then subset early and late ####

dge = DGEList(
  counts=exp_ov_filtered,
  samples=clinical_data)

# filtering
keep = filterByExpr(dge,design = clinical_data$figo_stage)

dge = dge[keep,,keep.lib.sizes=FALSE]

dge = calcNormFactors(dge,method="TMM")

dim(dge)

design=model.matrix(~clinical_data$age_at_diagnosis+
                      clinical_data$synchronous_malignancy+
                      clinical_data$race+
                      clinical_data$ethnicity+
                      clinical_data$treatments_radiation_treatment_or_therapy+clinical_data$figo_stage)

#RLE plot  ----

library('EDASeq')
library(Biobase)
library(scater)

sce1 <- SingleCellExperiment(assays = list(counts =dge$counts))

plotRLE(sce1, exprs_values = "counts", exprs_logged=FALSE, style = "minimal")+
  ggtitle("RLE")

norm_counts=cpm.DGEList(dge)

norm_counts=as.data.frame(norm_counts)

##### EXTRACT EARLY COUNTS ####
early_clinical=clinical_data%>%
  filter(figo_stage=="early")

rownames(early_clinical)=early_clinical$barcode

early_exp <- norm_counts[, colnames(norm_counts) %in% rownames(early_clinical)]

early_exp <- early_exp[, rownames(early_clinical)]

table(early_clinical$vital_status)

##### Survival function ####

perform_survival_analysis <- function(gene_list, norm_counts, survival_data) {
  results <- list()  # List to store results
  
  for (gene in gene_list) {
    # Check if the gene exists in norm_counts
    if (!(gene %in% rownames(norm_counts))) {
      warning(paste("Gene", gene, "not found in norm_counts. Skipping."))
      next
    }
    
    # Extract and prepare the gene expression data
    gene_expression <- norm_counts[gene, , drop = FALSE]
    gene_expression <- t(gene_expression)
    gene_expression <- as.data.frame(gene_expression)
    colnames(gene_expression)[1] <- gene
    
    # Check for missing values in gene_expression
    if (any(is.na(gene_expression))) {
      warning(paste("Gene expression data for", gene, "contains missing values. Skipping."))
      next
    }
    
    # Combine with clinical data
    gene_data <- cbind(gene_expression, survival_data)
    
    # Check if the merged data is valid
    if (nrow(gene_data) == 0 || any(is.na(gene_data$overall_survival)) || any(is.na(gene_data$deceased)) || any(is.na(gene_data[[gene]]))) {
      warning(paste("No valid data or missing values detected for gene", gene, ". Skipping."))
      next
    }
    
    # Compute median value and classify gene expression
    median_value <- median(gene_data[[gene]], na.rm = TRUE)
    gene_data$gene_class <- ifelse(gene_data[[gene]] >= median_value, "UP", "DOWN")
    
    # Perform Kaplan-Meier survival analysis
    fit <- tryCatch({
      survfit(Surv(overall_survival, deceased) ~ gene_class, data = gene_data)
    }, error = function(e) {
      warning(paste("Error in survfit for gene", gene, ":", e$message))
      return(NULL)
    })
    
    if (is.null(fit)) {
      next
    }
    
    # Extract p-value
    pval <- tryCatch({
      surv_pvalue(fit, data = gene_data)$pval
    }, error = function(e) {
      warning(paste("Error in surv_pvalue for gene", gene, ":", e$message))
      return(NA)
    })
    
    # Store results
    results[[gene]] <- list(
      p_value = pval,
      km_plot = tryCatch({
        ggsurvplot(fit, data = gene_data, pval = TRUE, risk.table = TRUE, title = paste(gene))
      }, error = function(e) {
        warning(paste("Error in ggsurvplot for gene", gene, ":", e$message))
        return(NULL)
      }),
      cox_model = tryCatch({
        summary(coxph(Surv(overall_survival, deceased) ~ gene_class, data = gene_data))
      }, error = function(e) {
        warning(paste("Error in coxph for gene", gene, ":", e$message))
        return(NULL)
      })
    )
    
    # Print p-value and Kaplan-Meier plot
    print(pval)
    if (!is.null(results[[gene]]$km_plot)) {
      print(results[[gene]]$km_plot)
    }
  }
  
  return(results)
}


###### LOAD RFE-SFS GENES ######
S1_genes=read_csv("D:/TEP/Ovarian cancer/Survival/final_common_S1.csv") 
S2_genes=read_csv("D:/TEP/Ovarian cancer/Survival/final_common_S2.csv")  
S3_genes=read_csv("D:/TEP/Ovarian cancer/Survival/final_common_S3.csv") 
S4_genes=read_csv("D:/TEP/Ovarian cancer/Survival/final_common_S4.csv") 

results1 <- perform_survival_analysis(S1_genes$EnsemblID, early_exp, early_clinical)

results1[["ENSG00000142546"]]
results1[["ENSG00000163219"]]


results2 <- perform_survival_analysis(S2_genes$EnsemblID, early_exp, early_clinical)


results3 <- perform_survival_analysis(S3_genes$EnsemblID, late_exp, late_clinical)
results3[["ENSG00000169895"]]
results3[["ENSG00000008130"]]

results4 <- perform_survival_analysis(S4_genes$EnsemblID, late_exp, late_clinical)

results4[["ENSG00000133742"]]
results4[["ENSG00000083312"]]

############ NOSIP #############

NOSIP <- early_exp[c("ENSG00000142546"),]
NOSIP = t(NOSIP)
NOSIP = as.data.frame(NOSIP)
colnames(NOSIP)[1] = "NOSIP"
NOSIP$NOSIP = as.numeric(NOSIP$NOSIP)
NOSIP_median_value = median(NOSIP$NOSIP)
print(NOSIP_median_value)
NOSIP = cbind(NOSIP, early_clinical)
NOSIP$group = ifelse(NOSIP$NOSIP >= NOSIP_median_value, "High Expression", "Low Expression")
fit = survfit(Surv(overall_survival, deceased) ~ group, data=NOSIP)
pval = surv_pvalue(fit, data=NOSIP)$pval
print(pval)
ggsurvplot(fit, 
           data = NOSIP, 
           pval = TRUE, 
           risk.table = TRUE, 
           title = "NOSIP",
           palette = c("red", "blue"),
           xlab = "Time(days)",          
           xscale = 1)

############ ARHGAP25 #############

ARHGAP25 <- early_exp[c("ENSG00000163219"),]
ARHGAP25 = t(ARHGAP25)
ARHGAP25 = as.data.frame(ARHGAP25)
colnames(ARHGAP25)[1] = "ARHGAP25"
ARHGAP25$ARHGAP25 = as.numeric(ARHGAP25$ARHGAP25)
ARHGAP25_median_value = median(ARHGAP25$ARHGAP25)
print(ARHGAP25_median_value)
ARHGAP25 = cbind(ARHGAP25, early_clinical)
ARHGAP25$group = ifelse(ARHGAP25$ARHGAP25 >= ARHGAP25_median_value, "High Expression", "Low Expression")
fit = survfit(Surv(overall_survival, deceased) ~ group, data=ARHGAP25)
pval = surv_pvalue(fit, data=ARHGAP25)$pval
print(pval)
ggsurvplot(fit, 
           data = ARHGAP25, 
           pval = TRUE, 
           risk.table = TRUE, 
           title = "ARHGAP25",
           palette = c("red", "blue"),
           xlab = "Time(days)",          
           xscale = 1)

###########  SYAP1 ############

SYAP1 <- late_exp[c("ENSG00000169895"),]
SYAP1 = t(SYAP1)
SYAP1 = as.data.frame(SYAP1)
colnames(SYAP1)[1] = "SYAP1"
SYAP1$SYAP1 = as.numeric(SYAP1$SYAP1)
SYAP1_median_value = median(SYAP1$SYAP1)
print(SYAP1_median_value)
SYAP1 = cbind(SYAP1, late_clinical)
SYAP1$group = ifelse(SYAP1$SYAP1 >= SYAP1_median_value, "High Expression", "Low Expression")
fit = survfit(Surv(overall_survival, deceased) ~ group, data=SYAP1)
pval = surv_pvalue(fit, data=SYAP1)$pval
print(pval)
ggsurvplot(fit, 
           data = SYAP1, 
           pval = TRUE, 
           risk.table = TRUE, 
           title = "SYAP1",
           palette = c("red", "blue"),
           xlab = "Time(days)",          
           xscale = 1)

############ NADK ############

NADK <- late_exp[c("ENSG00000008130"),]
NADK = t(NADK)
NADK = as.data.frame(NADK)
colnames(NADK)[1] = "NADK"
NADK$NADK = as.numeric(NADK$NADK)
NADK_median_value = median(NADK$NADK)
print(NADK_median_value)
NADK = cbind(NADK, late_clinical)
NADK$group = ifelse(NADK$NADK >= NADK_median_value, "High Expression", "Low Expression")
fit = survfit(Surv(overall_survival, deceased) ~ group, data=NADK)
pval = surv_pvalue(fit, data=NADK)$pval
print(pval)
ggsurvplot(fit, 
           data = NADK, 
           pval = TRUE, 
           risk.table = TRUE, 
           title = "NADK",
           palette = c("red", "blue"),
           xlab = "Time(days)",          
           xscale = 1)

############## CA1  ############

CA1 <- late_exp[c("ENSG00000133742"),]
CA1 = t(CA1)
CA1 = as.data.frame(CA1)
colnames(CA1)[1] = "CA1"
CA1$CA1 = as.numeric(CA1$CA1)
CA1_median_value = median(CA1$CA1)
print(CA1_median_value)
CA1 = cbind(CA1, late_clinical)
CA1$group = ifelse(CA1$CA1 >= CA1_median_value, "High Expression", "Low Expression")
fit = survfit(Surv(overall_survival, deceased) ~ group, data=CA1)
pval = surv_pvalue(fit, data=CA1)$pval
print(pval)
ggsurvplot(fit, 
           data = CA1, 
           pval = TRUE, 
           risk.table = TRUE, 
           title = "CA1",
           palette = c("red", "blue"),
           xlab = "Time(days)",          
           xscale = 1)

##############

TNPO1 <- late_exp[c("ENSG00000083312"),]  
TNPO1 = t(TNPO1)
TNPO1 = as.data.frame(TNPO1)
colnames(TNPO1)[1] = "TNPO1"
TNPO1$TNPO1 = as.numeric(TNPO1$TNPO1)
TNPO1_median_value = median(TNPO1$TNPO1)
print(TNPO1_median_value)
TNPO1 = cbind(TNPO1, late_clinical)
TNPO1$group = ifelse(TNPO1$TNPO1 >= TNPO1_median_value, "High Expression", "Low Expression")
fit = survfit(Surv(overall_survival, deceased) ~ group, data=TNPO1)
pval = surv_pvalue(fit, data=TNPO1)$pval
print(pval)
ggsurvplot(fit, 
           data = TNPO1, 
           pval = TRUE, 
           risk.table = TRUE, 
           title = "TNPO1",
           palette = c("red", "blue"),
           xlab = "Time(days)",          
           xscale = 1)


