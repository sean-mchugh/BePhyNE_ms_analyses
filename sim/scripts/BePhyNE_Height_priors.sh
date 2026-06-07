#!/bin/bash

RUN_ID=$1
HEIGHT_MEAN=$2
HEIGHT_SD=$3

Rscript Height_BePhyNE_sim.R "$RUN_ID" "$HEIGHT_MEAN" "$HEIGHT_SD"
