#!/bin/bash

type=$1
celltype=$2

/data/gc_sceqtl/sceqtl/eQTL/gtex-pipeline/qtl/src/combine_covariates.py /data/gc_sceqtl/sceqtl/eQTL/${type}/${celltype}_hp/covar/${celltype}_peer.PEER_covariates.txt /data/gc_sceqtl/sceqtl/eQTL/${type}/${celltype}_hp/covar/${celltype} \
--genotype_pcs /data/gc_sceqtl/sceqtl/eQTL/${type}/${celltype}_hp/covar/genotype_pca.txt \
--add_covariates /data/gc_sceqtl/sceqtl/eQTL/${type}/${celltype}_hp/covar/add_covariates1.txt

echo ${type}
echo ${celltype}
