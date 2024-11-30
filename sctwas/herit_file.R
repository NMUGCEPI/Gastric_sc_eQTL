library(tidyverse)
library(broom)
library(ggplot2)
library(qvalue)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

# Example
# cellLabel <- "epi"
# typeLabel <- "pit"

# Directory paths
exp.dir <- '/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig'
output.dir <- "/data/gc_sceqtl/herit"
cov.dir <- "/data/gc_sceqtl/sceqtl/eQTL"

# Input filenames
expression_filename <- sprintf("%s/%s/%s/exp_%s_mean_filter.xls", exp.dir, cellLabel, typeLabel, typeLabel)
pheno_filename <- sprintf("%s/%s/pheno/", output.dir, typeLabel)
region_filename <- sprintf("%s/%s/gene_region/", output.dir, typeLabel)
covariate_filename <- sprintf("%s/%s/%s_hp/covar/%s.combined_covariates.txt", cov.dir , cellLabel, typeLabel, typeLabel)

options(stringsAsFactors=F)
exp_mean_filter<-read.table(expression_filename,sep='\t',header=T,row.names=1)
exp_mean_filter$gene_name<-rownames(exp_mean_filter)
gtf_gene<-read.csv('/Public/wtp/blj/sceqtl/eQTL/GTF/gtf_gene_auto.csv',header=T)
exp_mean_filter1<-merge(exp_mean_filter,gtf_gene,by='gene_name') 
exp_bed<-exp_mean_filter1[,c((ncol(exp_mean_filter1)-3):ncol(exp_mean_filter1),2:(ncol(exp_mean_filter1)-4))] 
colnames(exp_bed)[1:4]<-c('chr','start','end','gene_id')
exp_bed$end<-exp_bed$start
exp_bed$start<-exp_bed$start-1
exp_bed$bps<-exp_bed$start-1000000
exp_bed$bpe<-exp_bed$end+1000000
exp_bed$bps1<-ifelse(exp_bed$bps<0,0,exp_bed$bps)
write.table(exp_bed[,c("chr","bps1","bpe","gene_id")],sprintf("%s/%s/Gene_region_for_GCTA_1MB", output.dir, typeLabel),row.names = F,col.names = F,quote = F)

#phenotype file
exp<-exp_bed[,c(4:(ncol(exp_bed)-3))]
texp<-as.data.frame(t(exp))
write.table(texp,sprintf("%s/%s/texp", output.dir, typeLabel),row.names = T,col.names = F,quote = F)
texp<-fread(sprintf("%s/%s/texp", output.dir, typeLabel))
fam<-read.table(sprintf("%s/%s/%s_info05_filter05.fam", output.dir, typeLabel, typeLabel))
colnames(fam)[1:2]<-c('array_id','sample_id')
colnames(texp)[1]<-"sample_id"
exp_fam<-merge(fam,texp,by='sample_id')
#write.table(exp_fam,sprintf("%s/%s/Gene_expression_used", output.dir, typeLabel),row.names=F,quote=F)

n<-ncol(exp_fam)
for(i in 7:n){
	write.table(exp_fam[,c(2,1,i)],file=paste(pheno_filename,colnames(exp_fam)[i],sep=''),row.names=F,col.names=F,quote=F)
}

gene<-read.table(sprintf("%s/%s/Gene_region_for_GCTA_1MB", output.dir, typeLabel))
gene$V1<-as.numeric(gene$V1)

n<-nrow(gene)
for(i in 1:n){
     write.table(gene[i,],file=paste(region_filename,gene[i,4],sep=''),col.names=F,row.names=F,quote=F)
}

gene1<-gene$V4
gene1<-t(gene1)
write.table(gene1,sprintf("%s/%s/genelist", output.dir, typeLabel),col.names=F,row.names=F,quote=F)

##covar
covar<-read.table(covariate_filename,h=T)
covar<-as.data.frame(t(covar))
write.table(covar,sprintf("%s/%s/covar", output.dir, typeLabel),row.names = T,col.names = F,quote = F)
covar<-fread(sprintf("%s/%s/covar", output.dir, typeLabel))
colnames(covar)[1]<-'sample_id'
covar<-merge(fam,covar,by='sample_id')
ccovar<-covar[,c(2,1,32,34)]
write.table(ccovar,sprintf("%s/%s/ccovar", output.dir, typeLabel),col.names = F,row.names = F,quote = F)
qcovar<-covar[,c(2,1,7:16,33)]
write.table(qcovar,sprintf("%s/%s/qcovar", output.dir, typeLabel),col.names = F,row.names = F,quote = F)

print('JOB IS DONE!')

quit()
