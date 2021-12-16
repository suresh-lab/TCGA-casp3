library(TCGAbiolinks)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(survival)
# load the pre-processed counts matrix

load("./TCGA-MK2_bootstrappedHR/LUAD_Norm_Filtered_IlluminaHiSeq.rda")

# load the clinical data

clinical <<- read.csv("rawdata/tcga-luad-clinical.processed_forCoxPH.csv",sep="\t")
dataFilt_LUAD <- as.data.frame(dataFilt_LUAD)
dataFilt_LUAD$gene <- rownames(dataFilt_LUAD)



RunModelWithGene <- function(gene){
  

Exp <- as.data.frame(t(dataFilt_LUAD[dataFilt_LUAD$gene==gene,]))
Exp$sample <- rownames(Exp)
colnames(Exp) <- c("gene_of_interest", "sample")
Exp <- Exp[!grepl("11A", Exp$sample),]
Exp <- Exp[grepl("01A", Exp$sample),]
Exp$sample.abbrev <- substring(Exp$sample, 1,12)
clinical_temp <- left_join(clinical, Exp, by=c("submitter_id"="sample.abbrev"))
clinical_temp$gene_of_interest <- as.numeric(clinical_temp$gene_of_interest)
bottom_quant <<- 0.33
top_quant <<- 0.66
clinical_temp <- clinical_temp[!is.na(clinical_temp$gene_of_interest),]
clinical_temp$gene_topQ <- sapply(clinical_temp$gene_of_interest, function (x) if(x > quantile(clinical_temp$gene_of_interest, top_quant))  return(1) else return(0))

return(summary(coxph(Surv(time, censor.2y) ~ gene_topQ+gender+smoking+age_at_diagnosis, id=submitter_id, data=clinical_temp[clinical_temp$ajcc_pathologic_stage_combined=="Early Stage",])))
}



genes <- dataFilt_LUAD$gene


# now we try all genes!
genes <- as.data.frame(genes)
genes$model_results <- lapply(genes$genes, RunModelWithGene)

save(genes, file=paste0("bootstrappedHR_2y.rda"))



