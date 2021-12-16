# process bootstrapped HRs
library(ggplot2)
library(ggrepel)
library(ggridges)
load("TCGA-MK2_bootstrappedHR/bootstrappedHR.rda")
genes$model_results_hr <- sapply(genes$model_results, function(x) return(x[[8]][1]))
genes$model_results_hr_lci <- sapply(genes$model_results, function(x) return(x[[8]][9]))
genes$model_results_hr_uci <- sapply(genes$model_results, function(x) return(x[[8]][13]))

genes <- genes[order(genes$model_results_hr),]
# let's plot!

genes$genes <- factor(genes$genes, levels=genes$genes)
genes$label <- sapply(genes$genes, function(x) if(x=="MAPKAPK2") return("MK2") else return(""))
genes$MK2_est <- mapply(function(gene, hr) if(gene=="MAPKAPK2") return(hr) else return(NA), genes$genes, genes$model_results_hr)
genes$MK2_lci <- mapply(function(gene, lci) if(gene=="MAPKAPK2") return(lci) else return(NA), genes$genes, genes$model_results_hr_lci)
genes$MK2_uci <- mapply(function(gene, uci) if(gene=="MAPKAPK2") return(uci) else return(NA), genes$genes, genes$model_results_hr_uci)

genes$sddevbelow <-sapply(genes$model_results_hr, function(x) if( x < quantile(genes$model_results_hr, prob=0.025)) return(1) else return(0))



# overall plot
ggplot(genes, aes(x=genes,y=model_results_hr,  ymin=model_results_hr_lci,label=label, ymax=model_results_hr_uci))+
geom_errorbar(color="grey")+geom_point(size=1)+geom_hline(yintercept = 1)+
geom_hline(yintercept=quantile(genes$model_results_hr, prob=0.025), linetype=2)+
geom_hline(yintercept=quantile(genes$model_results_hr, prob=0.975), linetype=2)+ 
theme(axis.text = element_text(size = 0.4))+
geom_label_repel(max.overlaps=10000)+
geom_point(inherit.aes=FALSE, aes(x=genes, y=MK2_est, color="red"), size=3)+
geom_errorbar(inherit.aes=FALSE, aes(color="red", x=genes, ymin=MK2_lci, ymax=MK2_uci), width=0)+
coord_flip()+
xlab("Gene (n=14,899)")+ylab("HR (95% CI) for survival at 1 year in multivariate Cox PH model")


# plot genes that are in the < 2.5% 
ggplot(genes[genes$sddevbelow==1,], aes(x=genes,y=model_results_hr,  ymin=model_results_hr_lci,label=label, ymax=model_results_hr_uci))+
geom_errorbar(color="grey")+geom_point(size=1)+geom_hline(yintercept = 1)+
geom_hline(yintercept=quantile(genes$model_results_hr, prob=0.025), linetype=2)+
geom_hline(yintercept=quantile(genes$model_results_hr, prob=0.975), linetype=2)+ 
theme(axis.text = element_text(size = 0.4))+
geom_label_repel(max.overlaps=10000)+
geom_point(inherit.aes=FALSE, aes(x=genes, y=MK2_est, color="red"), size=3)+
geom_errorbar(inherit.aes=FALSE, aes(color="red", x=genes, ymin=MK2_lci, ymax=MK2_uci), width=0)+
coord_flip()+
xlab("Gene (n=14,899)")+ylab("HR (95% CI) for survival at 1 year in multivariate Cox PH model")


# genes that are significant in the < 2.5 percentile
siggens_below2.5 <- genes[genes$sddevbelow==1 & genes$model_results_hr_uci < 1,]
write.table(as.matrix(siggens_below2.5[,c(1,3,4,5)]), quote=FALSE, col.names = TRUE, row.names = TRUE, file="./TCGA-MK2_bootstrappedHR/siggenes_2.5.csv")

# distribution plot
ggplot(genes, aes(x=model_results_hr))+geom_histogram(alpha=0.3, binwidth=0.01)+geom_density(aes(y=0.01 * ..count..))+
geom_segment(x=quantile(genes$model_results_hr, prob=0.025), xend=quantile(genes$model_results_hr, prob=0.025), y=0,yend=140*approxfun(density(genes$model_results_hr))(quantile(genes$model_results_hr, prob=0.025)), linetype=2)+
geom_segment(x=quantile(genes$model_results_hr, prob=0.975), xend=quantile(genes$model_results_hr, prob=0.975), y=0,yend=20*approxfun(density(genes$model_results_hr))(quantile(genes$model_results_hr, prob=0.0975)), linetype=2)+
geom_segment(x=genes[genes$genes=="MAPKAPK2",]$model_results_hr, xend=genes[genes$genes=="MAPKAPK2",]$model_results_hr,y=0,yend=100*approxfun(density(genes$model_results_hr))(genes[genes$genes=="MAPKAPK2",]$model_results_hr), linetype=1, color="red")+
geom_rug(alpha=0.5, size=0.4)+
xlab("HR in Cox PH Model")+
ylab("count")+
theme_bw()

sig_genes <- genes[(genes$model_results_hr < 1 & genes$model_results_hr_uci < 1) | (genes$model_results_hr > 1 & genes$model_results_hr_lci > 1), ]
ggplot(sig_genes, aes(x=genes,y=model_results_hr,  ymin=model_results_hr_lci,label=label, ymax=model_results_hr_uci))+geom_errorbar(color="grey")+geom_point(size=1)+geom_hline(yintercept = 1)+theme(axis.text = element_text(size = 0.9))+geom_label_repel(max.overlaps=10000)+coord_flip()+xlab("Gene (n=1524)")+ylab("HR (95% CI) for survival at 1 year in multivariate Cox PH model")
 
write.table(as.matrix(sig_genes[,c(1,3,4,5)]), col.names = TRUE, quote=FALSE, sep=",", file="genes_with_sig_HR.csv")
