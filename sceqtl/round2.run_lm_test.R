# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

# Example
# cellLabel <- "epi"
# typeLabel <- "parietal"

# Import libraries
library(tidyverse)
library(broom)
library(ggplot2)
library(qvalue)
library(matrixStats)
library(reshape2)
library(data.table)
library(dplyr)
library(magrittr)
library(future) 

# Directory paths
exp.dir <- '/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig'
genotype.dir <- '/data/gc_sceqtl/sceqtl/plinkfile208'
output.dir <- "/data/gc_sceqtl/sceqtl/eQTL"

# Input filenames
r1_residual_filename <- sprintf("%s/%s/%s/round1_residual_HP_df.xls", exp.dir, cellLabel, typeLabel)
r1_significant_SNPs <- sprintf("%s/%s/%s_hp/%s_qtl_sig_FDR05.xls",output.dir,cellLabel,typeLabel,typeLabel)
geneLoc_filename <- sprintf("%s/%s/%s/%s_geneloc.txt", exp.dir, cellLabel, typeLabel, typeLabel)
genotype_filename <- sprintf("%s/%s/%s_indiv.txt", genotype.dir, cellLabel, typeLabel)
snpLoc_filename <- sprintf("%s/%s/%s_snploc.txt", genotype.dir, cellLabel, typeLabel)

# Read in files
## Count matrix
r1_residual_df <- fread(r1_residual_filename)
print("r1_residaul_df info ...")
dim(r1_residual_df)
print(r1_residual_df[1:5,1:5])

## Significant SNPs file
r1_significant_SNPs_df <- fread(r1_significant_SNPs)
print("r1_significant_SNPs_df info...")
dim(r1_significant_SNPs_df)
print(r1_significant_SNPs_df[1:5,1:5])

## Genotype file
genotype_df <- fread(genotype_filename)
print("genotype_df info...")
dim(genotype_df)
print(genotype_df[1:5,1:5])

## Gene location file
geneLoc_df <- fread(geneLoc_filename)
print("geneLoc_df info...")
dim(geneLoc_df)
print(geneLoc_df[1:5,1:4])

## SNP location file
snpLoc_df <- fread(snpLoc_filename)
print("snpLoc_df info...")
dim(snpLoc_df)
print(snpLoc_df[1:5,1:3])

# Identify the top eSNP for each eGene 
eSNP1 <- r1_significant_SNPs_df %>%
  group_by(gene_id) %>%
  arrange(qvalue) %>%
  filter(row_number()==1)

print("eSNP1 info...")
dim(eSNP1)

fwrite(eSNP1,sprintf("%s/%s/%s_hp/%s_lead_lm_eSNP1.xls",output.dir,cellLabel,typeLabel,typeLabel),sep="\t",quote=F)

# Remaning eSNPs to test
eSNPs_to_test <- r1_significant_SNPs_df %>% 
    group_by(gene_id) %>%
    arrange(qvalue) %>%
    filter(row_number()!=1)  

## Subset residuals for the genes to be tested
dim(r1_residual_df)
sample_ids <- r1_residual_df$sampleid
gene_ids <- eSNP1$gene_id
r1_residual_df <- r1_residual_df %>% select (all_of(gene_ids))
r1_residual_df$sampleid <- sample_ids
dim(r1_residual_df)

# Subset genotype file for the significant SNPs
dim(genotype_df)
snpid <- r1_significant_SNPs_df$variant_id
genotype_df1 <- genotype_df[(genotype_df$SNP %in% snpid),]
genotype_df2 <- as.data.frame(t(genotype_df1))
colnames(genotype_df2)<-genotype_df1$SNP
genotype_df2<-genotype_df2[-1,]
genotype_df3=apply(genotype_df2,2,function(x) as.numeric(as.character(x)))
genotype_df3<-as.data.frame(genotype_df3)
genotype_df3$sampleid<-substr(rownames(genotype_df2),21,28)
dim(genotype_df3)
genotype_df3 <- genotype_df3 %>% ungroup()
snp_ids <- colnames(genotype_df3[-1])

# Find residuals after adjustment of lead SNP
calculate_adjusted_residuals <- function (x) {
  gene <- x

  # select gene to regress
  exprs_val <- r1_residual_df %>% select("sampleid", all_of(gene))

  # select SNP to add
  snp = as.character(eSNP1$variant_id[(eSNP1$gene_id==gene)])
  snp_genotype = genotype_df3 %>% select("sampleid", all_of(snp))

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,snp_genotype,by="sampleid")
  colnames(test_df)[2] <- "expression"
  colnames(test_df)[3] <- "genotype"

  # Generate model
  model <- lm(expression ~ genotype , data=test_df)
  residuals=resid(model)
  residuals  
}

options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 4)
adjusted_residual_mat <- future(sapply(gene_ids,calculate_adjusted_residuals))
adjusted_residual_mat <- value(adjusted_residual_mat)
mat<-unlist(adjusted_residual_mat)
mat<-as.data.frame(mat)
mat$gene<-substr(rownames(mat),1,15)
mat$sample<-substring(rownames(mat),17)
mat$order<-1:nrow(mat)
sid<-as.data.frame(sample_ids)
sid$sample<-rownames(sid)
mat1<-merge(mat,sid,by='sample',all.x=T,all.y=T)
mat1 <- mat1[order(mat1$order),]
colnames(mat1)[5]<-'sampleid'
fwrite(mat1,sprintf("%s/%s/%s/round2_adjusted_lm_residual_HP.xls", exp.dir, cellLabel, typeLabel),sep="\t",quote=F)

# Spearman's rank correlation 
# x is the data frame with chr and pos, y is snpid and geneid
lm_correlation <- function (x,y) {
  geneid <- y$gene_id
  snp <- y$variant_id
  
  # Select values to test
  res_val <- mat1 %>% select("sampleid","gene","mat") %>% subset(gene==geneid)
  genotype_val <- genotype_df3 %>% select("sampleid", all_of(snp))
  
  # Create a test matrix
  test_df <- left_join(res_val,genotype_val,by="sampleid")
  colnames(test_df) <- c("sampleid","gene","residual", "SNP")
  
  
  # Generate model
  model <- lm(residual ~ SNP, data=test_df)
  beta <- as.numeric(coefficients(model)[2])
  t_stat <- as.numeric(summary(model)$coefficients[,3][2])
  pvalue <- as.numeric(summary(model)$coefficients[,4][2])
  summary <- data.frame(geneid,snp,beta,t_stat,pvalue)
  summary
}

gene_snp_test_df <- eSNPs_to_test %>% select(variant_id,gene_id)
options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 8)
adjusted_lm_df <- future(gene_snp_test_df %>% group_by(variant_id,gene_id) %>% group_modify(lm_correlation))
adjusted_lm_df <- value(adjusted_lm_df)

# Calculate the qvalues for pvalues
pvalues <- adjusted_lm_df$pvalue
qobj <- qvalue(p = pvalues, pi0=1)

# Save regardless of significance
adjusted_lm_df <- adjusted_lm_df %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
gtf_gene<-read.csv('/data/gc_sceqtl/sceqtl/eQTL/GTF/gtf_gene_auto.csv',header=T)
adjusted_lm_df<-merge(adjusted_lm_df,gtf_gene,by='gene_id')
dim(adjusted_lm_df)
fwrite(adjusted_lm_df,sprintf("%s/%s/%s_hp/%s_round2_lm_results.xls",output.dir,cellLabel,typeLabel,typeLabel),sep="\t",quote=F)


# Save only significant
adjusted_lm_df_significant <- adjusted_lm_df[(adjusted_lm_df$qvalue < 0.05),]
nrow(adjusted_lm_df_significant)
fwrite(adjusted_lm_df_significant,sprintf("%s/%s/%s_hp/%s_round2_lm_sig_results.xls",output.dir,cellLabel,typeLabel,typeLabel),sep="\t",quote=F)

print('JOB IS DONE!')

quit()
