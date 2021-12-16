
# cutoff sensitivity analysis
library(survival)

clinical <- read.csv("./rawdata/tcga-luad-clinical.processed_forCoxPH.csv", sep="\t")

ChangeCutoff <- function(cutoff) {
  
 clinical$MK2_Expression_topQvar <- sapply(clinical$MK2.Expression, function (x) if(x > quantile(clinical$MK2.Expression, cutoff))  return(1) else return(0)) 
 clinical$MK2_Expression_topQvar <- as.factor(clinical$MK2_Expression_topQvar)
 return(summary(coxph(Surv(time, censor.2y) ~ MK2_Expression_topQvar+gender+smoking+age_at_diagnosis, id=submitter_id, data=clinical[clinical$ajcc_pathologic_stage_combined=="Early Stage",])))

  
}

cutoff <- seq(0.2,0.8,by=0.001)
cutoff <- as.data.frame(cutoff)
cutoff$model_results <- lapply(cutoff$cutoff, ChangeCutoff)
cutoff$hr <- sapply(cutoff$model_results, function(x) return(x[[7]][5]))
cutoff$lci <- sapply(cutoff$model_results, function(x) return(x[[8]][9]))
cutoff$uci <- sapply(cutoff$model_results, function(x) return(x[[8]][13]))

library(ggplot2)

ggplot(cutoff, aes(x=cutoff,y=hr, ymin=lci, ymax=uci))+geom_point()+theme_bw()+geom_hline(yintercept = 1)+geom_errorbar(width=0) 
ggplot(cutoff, aes(x=cutoff,y=hr, ymin=lci, ymax=uci))+geom_point()+theme_bw()+geom_errorbar(width=0, color="grey", alpha=0.3)+geom_smooth(method="loess")+geom_vline(xintercept=0.66, color="red", linetype=2)+xlim(0.3,0.8)
