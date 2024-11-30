options(stringsAsFactors=F)
library(data.table)
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

covariate_filename <- sprintf("/data/gc_sceqtl/sceqtl/eQTL/%s/%s_hp/covar/%s.combined_covariates15.txt", cellLabel, typeLabel, typeLabel)
covar<-read.table(covariate_filename,h=T)
covar<-as.data.frame(t(covar))
write.table(covar,sprintf("/data/gc_sceqtl/twas/covar/%s_covar_tmp",typeLabel),row.names = T,col.names = F,quote = F)
covar <- fread(sprintf("/data/gc_sceqtl/twas/covar/%s_covar_tmp",typeLabel))
covar <- cbind(FID = covar$ID, IID = covar$ID, covar[,2:ncol(covar)])
write.table(covar, sprintf("/data/gc_sceqtl/twas/covar/%s_covar",typeLabel), row.names = F, col.names = T, quote = F, sep = "\t")
