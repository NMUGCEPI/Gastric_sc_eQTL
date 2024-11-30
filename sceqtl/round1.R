# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

# Example
# cellLabel <- "epi"
# typeLabel <- "neck"

# Import libraries
library(tidyverse)
library(broom)
library(ggplot2)
library(qvalue)
library(matrixStats)
library(data.table)
library(dplyr)
library(magrittr)

# Directory paths
exp.dir <- '/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig'
genotype.dir <- '/data/gc_sceqtl/sceqtl/plinkfile208'
output.dir <- "/data/gc_sceqtl/sceqtl/eQTL"

# Input filenames
expression_filename <- sprintf("%s/%s/%s/%s_exp.txt", exp.dir, cellLabel, typeLabel, typeLabel)
geneLoc_filename <- sprintf("%s/%s/%s/%s_geneloc.txt", exp.dir, cellLabel, typeLabel, typeLabel)
genotype_filename <- sprintf("%s/%s/%s_indiv.txt", genotype.dir, cellLabel, typeLabel)
snpLoc_filename <- sprintf("%s/%s/%s_snploc.txt", genotype.dir, cellLabel, typeLabel)
covariate_filename <- sprintf("%s/%s/%s_hp/covar/%s.combined_covariates15.txt", output.dir, cellLabel, typeLabel, typeLabel)


# Read in files
##covar
covariate_df <- fread(covariate_filename,header=T)
covariate_df1<-as.data.frame(t(covariate_df))
colnames(covariate_df1)<-covariate_df$ID
covariate_df2<-covariate_df1[-1,]
covariate_df3=apply(covariate_df2,2,function(x) as.numeric(as.character(x)))
dim(covariate_df3)
covariate_df3<-as.data.frame(covariate_df3)
covariate_df3$ID<-rownames(covariate_df2)
print(covariate_df3[1:5,1:5])

## Count matrix
expression_df <- fread(expression_filename,header=T)
expression_df1<-as.data.frame(t(expression_df))
colnames(expression_df1)<-expression_df$ID
expression_df2<-expression_df1[-1,]
expression_df3=apply(expression_df2,2,function(x) as.numeric(as.character(x)))
dim(expression_df3)
expression_df3<-as.data.frame(expression_df3)
expression_df3$ID<-rownames(expression_df2)
print(expression_df3[1:5,1:5])

covariate_ids <- colnames(covariate_df3[-ncol(covariate_df3)])
gene_ids <- colnames(expression_df3[-ncol(expression_df3)])
sample_ids <- expression_df3$ID

# Find residuals for log transformed expressions
calculate_residuals <- function (x) {
  gene <- x

  # select gene to regress
  exprs_val <- expression_df3[,c("ID",gene)]

  # Create a test df by adding covariates
  test_df <- left_join(exprs_val,covariate_df3,by="ID")
  colnames(test_df)[2] <- "expression"

  # Generate model
  model <- lm(expression ~ pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + InferredCov1 + InferredCov2 + age + sex + HP, data=test_df)
  residuals=resid(model)
  residuals  
}


residual_mat <- sapply(gene_ids,calculate_residuals) 
rownames(residual_mat) <- sample_ids
residual_df <- data.frame(residual_mat, check.names=FALSE)
residual_df$sampleid <- sample_ids
print(residual_df[1:5,1:5])
dim(residual_df)
write.table(residual_df, sprintf("%s/%s/%s/round1_residual_HP_df.xls", exp.dir, cellLabel, typeLabel), row.names=T, col.names=T, sep='\t')

print('JOB IS DONE!')

quit()
