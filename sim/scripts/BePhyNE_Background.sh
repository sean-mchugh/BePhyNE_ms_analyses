#PATH=$PATH:/opt/R/4.2.3/lib/R
#
#
#
# 
# 
#Rscript  Background_BePhyNE_sim_new.R


#!/bin/bash

IFS="_" read -ra arr <<< "$LSB_JOBNAME"
unset IFS

BOOL_FLAG="${arr[0]}"
RUNNUM="${arr[1]}"

echo "Running job with RUN = $RUNNUM and FLAG = $BOOL_FLAG"

# Add R to PATH if needed
PATH=$PATH:/opt/R/4.2.3/lib/R

# Run the R script with extracted arguments
Rscript Background_BePhyNE_sim.R "$BOOL_FLAG"
