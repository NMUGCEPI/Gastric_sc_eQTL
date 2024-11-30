#!/bin/bash

type=$1
celltype=$2
Rscript /data/gc_sceqtl/twas/code/covar_twas.R ${type} ${celltype}
