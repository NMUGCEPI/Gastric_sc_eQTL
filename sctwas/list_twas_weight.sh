#!/bin/bash

type=$1
celltype=$2
Rscript /data/gc_sceqtl/twas/code/list_twas_weight.R ${type} ${celltype}
