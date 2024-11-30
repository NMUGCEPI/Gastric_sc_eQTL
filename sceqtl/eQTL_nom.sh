#!/bin/bash

type=$1
celltype=$2

/data/gc_sceqtl/sceqtl/fastqtl/python/run_FastQTL_threaded.py /data/gc_sceqtl/sceqtl/plinkfile208/${type}/${celltype}_filter05_info05.vcf.gz /data/gc_sceqtl/sceqtl/seuratall/allfilter_orig/${type}/${celltype}/exp_${celltype}_bed.bed.gz ${celltype}_eqtl \
--covariates /data/gc_sceqtl/sceqtl/eQTL/${type}/${celltype}_hp/covar/${celltype}.combined_covariates15.txt \
--window 1e6 --chunks 100 --threads 16 \
--output_dir /data/gc_sceqtl/sceqtl/eQTL/${type}/${celltype}_hp