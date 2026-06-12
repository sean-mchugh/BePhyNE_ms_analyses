PATH=$PATH:/opt/R/4.2.3/lib/R



HEIGHT_MEAN=$1
HEIGHT_SD=$2


Rscript  posterior_correlation_summary.R  "$HEIGHT_MEAN" "$HEIGHT_SD"
