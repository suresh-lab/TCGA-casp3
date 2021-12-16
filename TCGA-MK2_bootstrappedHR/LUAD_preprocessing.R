# first load the RNAseq data from a rda file. See APIAccess.rmd for details for how this was obtained using GDC
load("TCGA-MK2_bootstrappedHR/LUADIllumina_HiSeq.rda")

library(TCGAbiolinks)
dataPrep_LUAD  <- TCGAanalyze_Preprocessing(object=data, cor.cut=0.6, datatype="raw_count")
dataNorm_LUAD <- TCGAanalyze_Normalization(tabDF=dataPrep_LUAD, geneInfo=TCGAbiolinks::geneInfo, method="gcContent")
dataFilt_LUAD <- TCGAanalyze_Filtering(tabDF=dataNorm_LUAD, method="quantile", qnt.cut=0.25)
# write the normalized, filtered data
save(dataFilt_LUAD, file=paste0("LUAD_Norm_Filtered_IlluminaHiSeq.rda"))

clinical_data <- read.csv("rawdata/tcga-luad-clinical.processed.csv")

