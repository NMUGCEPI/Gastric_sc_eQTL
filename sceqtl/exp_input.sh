#!/bin/bash

type=$1
celltype=$2
Rscript /data/gc_sceqtl/sceqtl/eQTL/exp_input.R  ${type} ${celltype}
