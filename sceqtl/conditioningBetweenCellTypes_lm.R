# Call variables
args = commandArgs(trailingOnly=TRUE)
cellLabel1 <- args[1]
print(cellLabel1)
typeLabel1 <- args[2]
print(typeLabel1)
cellLabel2 <- args[3]
print(cellLabel2)
typeLabel2 <- args[4]
print(typeLabel2)

# Example
# cellLabel1 <- "epi"
# typeLabel1 <- "chief"
# cellLabel2 <- "epi"
# typeLabel2 <- "pit"

# Import libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(broom)
library(future)

# Directory paths
exp.dir <- '/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig'
genotype.dir <- '/data/gc_sceqtl/sceqtl/plinkfile208'
output.dir <- "/data/gc_sceqtl/sceqtl/eQTL"

# Input filenames
celltype1_filename <- sprintf("%s/%s/%s_hp/%s_lead_eSNP1.xls",output.dir,cellLabel1,typeLabel1,typeLabel1)
celltype2_filename <- sprintf("%s/%s/%s_hp/%s_lead_eSNP1.xls",output.dir,cellLabel2,typeLabel2,typeLabel2)
residuals_filename <- sprintf("%s/%s/%s/round1_residual_HP_df.xls", exp.dir, cellLabel1, typeLabel1) ##校正协变量表达残差
genotype_filename <- sprintf("/data/gc_sceqtl/sceqtl/plinkfile208/ALLCHR_info05_indiv.txt")  ##读取所有基因型

# Read in files
celltype1_df <- fread(celltype1_filename)
celltype1_df %>% head()
celltype2_df <- fread(celltype2_filename)
celltype2_df %>% head()
residuals_df <- fread(residuals_filename)
sample_ids <- residuals_df$sampleid
 
# Add prefix to colnames
colnames(celltype1_df) <- paste("celltype1", colnames(celltype1_df), sep = "_")
colnames(celltype2_df) <- paste("celltype2", colnames(celltype2_df), sep = "_")
colnames(celltype1_df)[1] <- "geneid"
colnames(celltype2_df)[1] <- "geneid"

# Merge top snps to see the overlap
top_snps_df <- full_join(celltype1_df, celltype2_df, by="geneid")
top_snps_df %>% head()

# Check if celltype1_snpid and celltype2_snpid are the same
top_snps_df <- top_snps_df %>% mutate(agreement = celltype1_variant_id==celltype2_variant_id)

# Number of eGenes with same lead SNP in two cell types 
print(nrow(top_snps_df %>% filter(agreement==TRUE))) 
same_top_snps_df <- top_snps_df %>% filter(agreement==TRUE)
same_egene_df <- top_snps_df %>% filter(agreement %in% c('TRUE','FALSE'))
fwrite(same_top_snps_df, sprintf("%s/type_share_lm/%s_vs_%s_same_top_snps.xls", output.dir, typeLabel1, typeLabel2), sep="\t", row.names=F, col.names=T, quote=F)
fwrite(same_egene_df, sprintf("%s/type_share_lm/%s_vs_%s_same_egene.xls", output.dir, typeLabel1, typeLabel2), sep="\t", row.names=F, col.names=T, quote=F)

# Subset those agree=0 to calculate new betas
gene_ids <- same_egene_df$geneid
No_eGenes <- print(nrow(same_egene_df))
print(sprintf("There is %s overlapping SNPs!", No_eGenes))
stopifnot(No_eGenes!= 0)

# Read the genotype file
genotype_df <- fread(genotype_filename)
snpid1 <- same_egene_df$celltype1_variant_id
snpid1 <- as.data.frame(snpid1)
colnames(snpid1) <- 'snpid'
snpid2 <- same_egene_df$celltype2_variant_id
snpid2 <- as.data.frame(snpid2)
colnames(snpid2) <- 'snpid'
snpid <- rbind(snpid1,snpid2)
genotype_df1 <- genotype_df[(genotype_df$SNP %in% snpid$snpid),]
genotype_df2 <- as.data.frame(t(genotype_df1))
colnames(genotype_df2)<-genotype_df1$SNP
genotype_df2<-genotype_df2[-1,]
genotype_df3=apply(genotype_df2,2,function(x) as.numeric(as.character(x)))
genotype_df3<-as.data.frame(genotype_df3)
genotype_df3$sampleid<-substr(rownames(genotype_df2),21,28)
dim(genotype_df3)
rm(genotype_df)


# Rename the files for conditioning function
# c1_exp=expression1_df
# df_snps=same_egene_df
# c2_geno=genotype2_df
# c1_covs=covariates1_df

calculate_pval_snp2 <- function(x) {
    genename <- x
    c1 <- typeLabel1
    c2 <- typeLabel2
    c1_exp <- residuals_df
    df_snps <- same_egene_df 
    df_geno <- genotype_df3
	
    # get expression of the gene in cell type1
    expression_celltype1 <- c1_exp %>% select("sampleid",all_of(genename))
    colnames(expression_celltype1)[2]<- "expression"
	
    # identify the snps for testing
    celltype1_lead_snp <- df_snps$celltype1_variant_id[df_snps$geneid==genename]
    celltype2_lead_snp <- df_snps$celltype2_variant_id[df_snps$geneid==genename]
	
    # get the genotype of celltype2_snp in celltype1
    genotype_celltype2_snp <- df_geno %>% select("sampleid",all_of(celltype2_lead_snp))
    colnames(genotype_celltype2_snp)[2] <- "genotype_c2snp"
    df <- merge(expression_celltype1,genotype_celltype2_snp, by="sampleid")
    model <- lm(expression ~ genotype_c2snp, data=df)
    beta_c2snp <- as.numeric(coefficients(model)[2])
    t_stat_c2snp <- as.numeric(summary(model)$coefficients[,3][2])
    pvalue_c2snp <- as.numeric(summary(model)$coefficients[,4][2])
    summary <- data.frame(genename,typeLabel1,typeLabel2,celltype1_lead_snp,celltype2_lead_snp,beta_c2snp,t_stat_c2snp,pvalue_c2snp)
    summary
}

# Apply regression function to all eGenes
options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 4)
exp_c1_snp_c2 <- future(lapply(gene_ids,calculate_pval_snp2))
exp_c1_snp_c2 <- value(exp_c1_snp_c2)
results_exp_c1_snp_c2<- do.call(rbind.data.frame,exp_c1_snp_c2)
fwrite(results_exp_c1_snp_c2, sprintf("%s/type_share_lm/%s_vs_%s_linear_model_summary.xls", output.dir, typeLabel1, typeLabel2), sep="\t", row.names=F, col.names=T, quote=F)

# Calculate residuals for each individual
calculate_residuals <- function(x) {
    genename <- x
    c1 <- typeLabel1
    c2 <- typeLabel2
    c1_exp <- residuals_df
    df_snps <- same_egene_df 
    df_geno <- genotype_df3
	
    # get expression of the gene in cell type1
    expression_celltype1 <- c1_exp %>% select("sampleid",all_of(genename))
    colnames(expression_celltype1)[2]<- "expression"
	
    # identify the snps for testing
    celltype1_lead_snp <- df_snps$celltype1_variant_id[df_snps$geneid==genename]
    celltype2_lead_snp <- df_snps$celltype2_variant_id[df_snps$geneid==genename]
	
    # get the genotype of celltype2_snp in celltype1
    genotype_celltype2_snp <- df_geno %>% select("sampleid",all_of(celltype2_lead_snp))
    colnames(genotype_celltype2_snp)[2] <- "genotype_c2snp"
    df <- merge(expression_celltype1,genotype_celltype2_snp, by="sampleid")
    model <- lm(expression ~ genotype_c2snp, data=df)
    residuals_all <- resid(model)
    residuals_all 
}


res_expression <- sapply(gene_ids, calculate_residuals)
mat<-unlist(res_expression)
mat<-as.data.frame(mat)
mat$gene<-substr(rownames(mat),1,15)
mat$sample<-substring(rownames(mat),17)
mat$order<-1:nrow(mat)
sid<-as.data.frame(sample_ids)
sid$sample<-rownames(sid)
mat1<-merge(mat,sid,by='sample',all.x=T,all.y=T)
mat1 <- mat1[order(mat1$order),]
colnames(mat1)[5]<-'sampleid'
fwrite(mat1, sprintf("%s/type_share_lm/%s_vs_%s_linear_model_residuals.xls", output.dir, typeLabel1, typeLabel2), sep="\t", row.names=F, col.names=T, quote=F)

# lm test 
# x is the data frame with chr and pos, y is snpid and geneid
lm_correlation <- function (x,y) {
  gene_id <- y$geneid
  snp <- y$celltype1_variant_id
  
  # Select values to test
  res_val <- mat1 %>% select("sampleid","gene","mat") %>% subset(gene==gene_id)
  genotype_val <- genotype_df3 %>% select("sampleid",all_of(snp))
  
  # Create a test matrix
  test_df <- left_join(res_val,genotype_val,by="sampleid")
  colnames(test_df) <- c("sampleid","gene","residual", "SNP")
  
  # Generate model
  model <- lm(residual ~ SNP, data=test_df)
  beta <- as.numeric(coefficients(model)[2])
  t_stat <- as.numeric(summary(model)$coefficients[,3][2])
  pvalue <- as.numeric(summary(model)$coefficients[,4][2])
  summary <- data.frame(gene_id,snp,beta,t_stat,pvalue)
  summary
}

options(future.globals.maxSize = 1 * 1024^3)
plan("multicore", workers = 4)
adjusted_lm_df <- future(same_egene_df %>% group_by(celltype1_variant_id,geneid) %>% group_modify(lm_correlation))
adjusted_lm_df <- value(adjusted_lm_df)
colnames(adjusted_lm_df)[5] <- "c1_new_beta"
colnames(adjusted_lm_df)[6] <- "c1_new_S_statistics"
colnames(adjusted_lm_df)[7] <- "c1_new_p.value"
adjusted_lm_df <- adjusted_lm_df %>% ungroup() %>% select("geneid", "c1_new_beta",
    "c1_new_S_statistics", "c1_new_p.value")
same_egene_df <- full_join(same_egene_df,adjusted_lm_df,by="geneid")
write.table(same_egene_df, sprintf("%s/type_share_lm/%s_vs_%s_lm_summary.xls", output.dir, typeLabel1, typeLabel2), sep="\t", row.names=F, col.names=T, quote=F)

print('JOB IS DONE!')

quit()

