library(stringr)
library(coda)
library(BePhyNE, lib = "packages")

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)

tolerance_mean <- as.numeric(args[1])
tolerance_sd   <- as.numeric(args[2])

cat("tolerance_mean:", tolerance_mean, "\n")
cat("tolerance_sd:", tolerance_sd, "\n")

tolerance_mean_label <- gsub(pattern = ".", replacement = "p", tolerance_mean, fixed = TRUE)
tolerance_sd_label   <- gsub(pattern = ".", replacement = "p", tolerance_sd,   fixed = TRUE)

run_dir <- paste0(
  "Neww_New_New_Height_runs/Height_runs_mean_",
  tolerance_mean_label,
  "_sd_",
  tolerance_sd_label,
  "/"
)

save_dir <- "no_burnin_Neww_New_New_Height_runs"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

n_spp <- 200

comparisons <- c(
  "optimum_vs_tolerance",
  "breadth_vs_tolerance",
  "optimum_vs_breadth",
  "optimum_vs_tolerance_ft",
  "breadth_vs_tolerance_ft"
)

mcmc_files <- list.files(run_dir)
mcmc_files_list <- mcmc_files[grep("log.pars.log", mcmc_files)]

if (length(mcmc_files_list) == 0) {
  stop("No log.pars.log files found in: ", run_dir)
}

make_trait_vectors <- function(logdf, pred, sp) {

  opt_col       <- paste0("pred_", pred, "_dat.opt_t", sp)
  breadth_col   <- paste0("pred_", pred, "_dat.brdth_t", sp)
  tolerance_col <- paste0("pred_", pred, "_dat.tol_t", sp)

  opt          <- logdf[[opt_col]]
  breadth      <- exp(logdf[[breadth_col]])
  tolerance_ft <- logdf[[tolerance_col]]
  tolerance    <- 0.05 + 0.95 * (exp(-1 * tolerance_ft) / (1 + exp(-1 * tolerance_ft)))

  list(
    opt = opt,
    breadth = breadth,
    tolerance = tolerance,
    tolerance_ft = tolerance_ft
  )
}

get_xy_for_comparison <- function(traits, comp) {

  if (comp == "optimum_vs_tolerance") {
    return(list(x = traits$opt, y = traits$tolerance))
  }

  if (comp == "breadth_vs_tolerance") {
    return(list(x = traits$breadth, y = traits$tolerance))
  }

  if (comp == "optimum_vs_breadth") {
    return(list(x = traits$opt, y = traits$breadth))
  }

  if (comp == "optimum_vs_tolerance_ft") {
    return(list(x = traits$opt, y = traits$tolerance_ft))
  }

  if (comp == "breadth_vs_tolerance_ft") {
    return(list(x = traits$breadth, y = traits$tolerance_ft))
  }

  stop("Unknown comparison: ", comp)
}

posterior_cor_list <- list()

## ------------------------------------------------------------
## First pass: compute within-chain correlations
##
## Each row is one:
##   MCMC replicate x predictor x species x comparison
## ------------------------------------------------------------

for (k in seq_along(mcmc_files_list)) {

  print(k)

  log_filename <- mcmc_files_list[[k]]
  logdf <- read_BePhyNE_log(paste0(run_dir, log_filename))

  #exclude burnin
  logdf = logdf[-c(1:(nrow(logdf)/2)),]
  cor_rows <- list()
  row_i <- 1

  for (pred in 1:2) {

    for (sp in 1:n_spp) {

      traits <- make_trait_vectors(logdf, pred, sp)

      for (comp in comparisons) {

        xy <- get_xy_for_comparison(traits, comp)

        cor_rows[[row_i]] <- data.frame(
          rep = k,
          pred = pred,
          species = sp,
          comparison = comp,
          cor = suppressWarnings(cor(xy$x, xy$y, use = "complete.obs"))
        )

        row_i <- row_i + 1
      }
    }
  }

  posterior_cor_list[[k]] <- do.call(rbind, cor_rows)
}

posterior_cor_summary <- do.call(rbind, posterior_cor_list)

## ------------------------------------------------------------
## Pick one representative chain per comparison
##
## For each comparison, choose the rep/pred/species whose within-chain
## correlation is closest to the median correlation for that treatment.
## ------------------------------------------------------------

representative_rows <- list()

for (comp in comparisons) {

  d <- posterior_cor_summary[
    posterior_cor_summary$comparison == comp &
      is.finite(posterior_cor_summary$cor),
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    next
  }

  med_cor <- median(d$cor, na.rm = TRUE)
  pick_i <- which.min(abs(d$cor - med_cor))

  representative_rows[[comp]] <- data.frame(
    comparison = comp,
    median_cor = med_cor,
    chosen_cor = d$cor[pick_i],
    rep = d$rep[pick_i],
    pred = d$pred[pick_i],
    species = d$species[pick_i]
  )
}

representative_correlations <- do.call(rbind, representative_rows)
row.names(representative_correlations) <- NULL

print(representative_correlations)

## ------------------------------------------------------------
## Second pass: extract the full posterior chain for each representative
## ------------------------------------------------------------

posterior_sample_list <- list()

for (i in seq_len(nrow(representative_correlations))) {

  comp <- representative_correlations$comparison[i]
  k    <- representative_correlations$rep[i]
  pred <- representative_correlations$pred[i]
  sp   <- representative_correlations$species[i]

  log_filename <- mcmc_files_list[[k]]
  logdf <- read_BePhyNE_log(paste0(run_dir, log_filename))
  
  logdf = logdf[-c(1:(nrow(logdf)/2)),]
  
  traits <- make_trait_vectors(logdf, pred, sp)
  xy <- get_xy_for_comparison(traits, comp)

  posterior_sample_list[[i]] <- data.frame(
    rep = k,
    pred = pred,
    species = sp,
    comparison = comp,
    draw = seq_along(xy$x),
    x = xy$x,
    y = xy$y,
    median_cor = representative_correlations$median_cor[i],
    chosen_cor = representative_correlations$chosen_cor[i]
  )
}

posterior_representative_samples <- do.call(rbind, posterior_sample_list)

## ------------------------------------------------------------
## Save RDS files
## ------------------------------------------------------------

saveRDS(
  posterior_cor_summary,
  paste0(
    save_dir,
    "/",
    tolerance_mean_label,
    "_",
    tolerance_sd_label,
    "_trait_height_posterior_samples_corr_summary.rds"
  )
)

saveRDS(
  posterior_representative_samples,
  paste0(
    save_dir,
    "/",
    tolerance_mean_label,
    "_",
    tolerance_sd_label,
    "_trait_height_posterior_representative_samples.rds"
  )
)

saveRDS(
  representative_correlations,
  paste0(
    save_dir,
    "/",
    tolerance_mean_label,
    "_",
    tolerance_sd_label,
    "_trait_height_median_correlation_representatives.rds"
  )
)

cat("Saved correlation summary and median-correlation representative posterior samples to: ",
    save_dir, "\n", sep = "")
