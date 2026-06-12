#!/bin/bash


HEIGHT_MEAN=$1
HEIGHT_SD=$2

Rscript BePhyNE_Height_collect_priors.R "$HEIGHT_MEAN" "$HEIGHT_SD"

#Rscript posterior_correlation_summary.R "$HEIGHT_MEAN" "$HEIGHT_SD"
