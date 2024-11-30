#exp input
# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

# Example
# cellLabel <- "Tcell"
# typeLabel <- "Tcell"

# Import libraries
library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
exp.dir <- '/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig'
rdata_filename <- sprintf("%s/%s/%s/%s.RDS", exp.dir, cellLabel, typeLabel, typeLabel)
##lognorm
type <- readRDS(rdata_filename)
lognorm_data<-GetAssayData(type,slot="data")
at_least_one <- apply(lognorm_data, 1, function(x) sum(x>0))
at_least_one <- as.data.frame(at_least_one)
write.table(at_least_one, sprintf("%s/%s/%s/gene_exp_all%s.xls", exp.dir, cellLabel, typeLabel, typeLabel), row.names=T, col.names=T, sep='\t')
#all<-read.table(sprintf("%s/%s/unlog/%s/gene_exp_all%s.xls", exp.dir, cellLabel, typeLabel, typeLabel),h=T,stringsAsFactors = F)
keep_gene = subset(at_least_one, at_least_one > ncol(lognorm_data)/100)
lognorm_data<-as.data.frame(lognorm_data)
lognorm_data<-as.data.frame(t(lognorm_data))
lognorm_data$id<-rownames(lognorm_data)
meta <- type@meta.data[,c("ID", "seurat_clusters")]
meta$id<-rownames(meta)
lognorm_data<-merge(lognorm_data,meta,by='id')
N=ncol(lognorm_data)-2
exp_mean<-aggregate(lognorm_data[,2:N], by=list(group=lognorm_data$ID),mean)
rownames(exp_mean)<-exp_mean$group
exp_mean<-as.data.frame(t(exp_mean[,2:ncol(exp_mean)]))
fwrite(exp_mean, sprintf("%s/%s/%s/exp_%s_mean.xls", exp.dir, cellLabel, typeLabel, typeLabel), row.names=T, col.names=T, sep='\t')


#keep_sample = read.table(sprintf("%s/%s/unlog/%s/%s_keep_sample.xls", exp.dir, cellLabel, typeLabel, typeLabel),h=T,stringsAsFactors = F)
a<-table(type$ID) >= 5
keep_sample<-as.data.frame(a)
write.table(keep_sample, sprintf("%s/%s/%s/%s_keep_sample.xls", exp.dir, cellLabel, typeLabel, typeLabel), row.names=T, col.names=T, sep='\t')
keep_sample<-subset(keep_sample,a=='TRUE')
exp_mean_filter<-exp_mean[rownames(keep_gene),rownames(keep_sample)]
fwrite(exp_mean_filter, sprintf("%s/%s/%s/exp_%s_mean_filter.xls", exp.dir, cellLabel, typeLabel, typeLabel), row.names=T, col.names=T, sep='\t')

##output
exp_filename <- sprintf("%s/%s/%s/exp_%s_mean_filter.xls", exp.dir, cellLabel, typeLabel, typeLabel)
exp_mean_filter <- read.table(exp_filename,h=T,row.names=1,stringsAsFactors=F)
means = apply(exp_mean_filter,1,mean)
all_mean <- as.data.frame(means)
all_mean$gene<-rownames(all_mean)
all_mean$cell_labels <- typeLabel
fwrite(all_mean,sprintf("%s/%s/%s/expall_mean.xls", exp.dir, cellLabel, typeLabel),sep="\t",quote=F)
exp_mean_filter$gene_name<-rownames(exp_mean_filter)
gtf_gene<-read.csv('/data/gc_sceqtl/sceqtl/eQTL/GTF/gtf_gene_auto.csv',header=T)
exp_mean_filter1<-merge(exp_mean_filter,gtf_gene,by='gene_name') 
exp_bed<-exp_mean_filter1[,c((ncol(exp_mean_filter1)-3):ncol(exp_mean_filter1),2:(ncol(exp_mean_filter1)-4))] 
colnames(exp_bed)[1:4]<-c('chr','start','end','ID')
exp_bed$end<-exp_bed$start
exp_bed$start<-exp_bed$start-1
exp_bed$chr<-as.numeric(exp_bed$chr)
exp_bed<-exp_bed[order(exp_bed$chr,exp_bed$start),]
write.table(exp_bed[,c(4,1:3)], sprintf("%s/%s/%s/%s_geneloc.txt", exp.dir, cellLabel, typeLabel, typeLabel),row.names=F, col.names=T,quote =F, sep='\t')
write.table(exp_bed[,c(4: ncol(exp_bed))], sprintf("%s/%s/%s/%s_exp.txt", exp.dir, cellLabel, typeLabel, typeLabel), row.names=F, col.names=T,quote =F, sep='\t')
colnames(exp_bed)[1]<-'#Chr'
write.table(exp_bed, sprintf("%s/%s/%s/exp_%s_bed.bed", exp.dir, cellLabel, typeLabel, typeLabel), row.names=F, col.names=T,quote =F, sep='\t')
dim(exp_bed)

print('JOB IS DONE!')

quit()
