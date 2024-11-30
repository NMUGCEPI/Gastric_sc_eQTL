#!/bin/bash

type=$1
celltype=$2
Rscript /data/gc_sceqtl/sceqtl/eQTL/round2.run_lm_test.R  ${type} ${celltype} 
