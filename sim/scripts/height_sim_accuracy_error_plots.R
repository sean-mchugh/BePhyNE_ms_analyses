


## Pooled-across-predictor treatment panels
## Assumes MCMC_summary_list objects have already been processed and saved.
## Produces one combined line/error PNG per pooled parameter.
## For each parameter:
##   rows = sd_vec
##   columns = mean values, with line plot then error-density plot for each mean
##   pooled over predictor 1 and predictor 2

height_vec <- c("0.5", "0.7", "0.95")
sd_vec     <- c("0.15", "0.5", "1")

n_per <- 200
dir <- "~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/plots/"
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

## If your processed RDS files use a different pattern, change only this function.
summary_file_for_treatment <- function(height_mean, height_sd) {
  height_mean_label <- gsub(pattern = ".", replacement = "p", height_mean, fixed = TRUE)
  height_sd_label   <- gsub(pattern = ".", replacement = "p", height_sd,   fixed = TRUE)

  candidates <- c(
    paste0("Height_runs_mean_", height_mean_label, "_sd_", height_sd_label,
           "_MCMC_summary_list_height.rds"),
    paste0("Height_runs_mean_", height_mean_label, "_sd_", height_sd_label,
           "/Height_runs_mean_", height_mean_label, "_sd_", height_sd_label,
           "_MCMC_summary_list_height.rds"),
    paste0("Height_runs_mean_", height_mean_label, "_sd_", height_sd_label,
           "/MCMC_summary_list_height.rds"),
    paste0("MCMC_summary_list_height_mean_", height_mean_label,
           "_sd_", height_sd_label, ".rds")
  )
  
  candidates <- c(
    paste0("outfiles/Height_runs_mean_", height_mean_label, "_sd_", height_sd_label,
           "_MCMC_summary_list_height.rds")
    #paste0("MCMC_summary_list_height_mean_", height_mean_label,
    #       "_sd_", height_sd_label, ".rds")
  )
  

  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Could not find processed summary RDS for mean=", height_mean,
         ", sd=", height_sd,
         ". Checked:\n", paste(candidates, collapse = "\n"))
  }

  hit[[1]]
}

## Same Greek/plotmath labels as the line plots, but pooled across predictors.
pooled_labels <- list(
  h     = expression(italic(h)),
  w     = expression(italic(omega)),
  c     = expression(italic(theta)),
  A_C   = expression(italic(A)[italic(theta)]),
  Rsd_C = expression(italic(sigma)[italic(theta)]),
  A_W   = expression(italic(A)[italic(omega)]),
  Rsd_W = expression(italic(sigma)[italic(omega)]),
  Rcor  = expression(italic(R)[italic(COR)])
)

## File-safe names, in the plotting order you want.
pooled_order <- c("t", "w", "c", "A_C", "Rsd_C", "A_W", "Rsd_W", "Rcor")

pool_post_out <- function(MCMC_summary, par_name) {
  x <- MCMC_summary[[1]]

  if (par_name == "h") {
    return(rbind(x$full_traits_h_1, x$full_traits_h_2))
  }
  if (par_name == "w") {
    return(rbind(x$full_traits_w_1, x$full_traits_w_2))
  }
  if (par_name == "c") {
    return(rbind(x$full_traits_c_1, x$full_traits_c_2))
  }
  if (par_name == "A_C") {
    return(rbind(x$full_traits_A_C_1, x$full_traits_A_C_2))
  }
  if (par_name == "Rsd_C") {
    return(rbind(x$full_traits_Rsd_C_1, x$full_traits_Rsd_C_2))
  }
  if (par_name == "A_W") {
    return(rbind(x$full_traits_A_W_1, x$full_traits_A_W_2))
  }
  if (par_name == "Rsd_W") {
    return(rbind(x$full_traits_Rsd_W_1, x$full_traits_Rsd_W_2))
  }
  if (par_name == "Rcor") {
    return(rbind(x$full_traits_Rcor_1, x$full_traits_Rcor_2))
  }

  stop("Unknown par_name: ", par_name)
}

percent_error_from_post <- function(post_out) {
  (post_out[, 1] - post_out[, 2]) #/ post_out[, 1]
}

read_treatment_summary <- function(height_mean, height_sd) {
  readRDS(summary_file_for_treatment(height_mean, height_sd))
}

plot_post_line <- function(post_out, treatment_title, point_size = 0.35) {
  good <- is.finite(post_out[, 1]) & is.finite(post_out[, 2])
  post_out <- post_out[good, , drop = FALSE]

  if (nrow(post_out) == 0) {
    plot.new()
    title(main = treatment_title)
    text(0.5, 0.5, "no finite points")
    return(invisible(NULL))
  }

  min_now <- min(c(post_out[, 1], post_out[, 2]), na.rm = TRUE) - 0.1
  max_now <- max(c(post_out[, 1], post_out[, 2]), na.rm = TRUE) + 0.1
  lim_now <- c(min_now, max_now)

  plot(NULL, NULL,
       ylab = "true",
       xlab = "posterior median",
       main = paste0(treatment_title, "\nline"),
       xlim = lim_now,
       ylim = lim_now,
       cex.main = 0.95)

  points(post_out[, 1], post_out[, 2], pch = 16, cex = point_size, col = 3)

  if (nrow(post_out) > 1) {
    abline(lm(post_out[, 2] ~ post_out[, 1]), col = 3)
  }

  abline(0, 1, col = "grey")
}

plot_error_density <- function(post_out, treatment_title) {
  ## Percent error is still (posterior median - true) / posterior median,
  ## matching the old script. The old barplot grouped errors by n_per = 200;
  ## that collapses global parameters like Rcor into one giant bar. This plots
  ## the distribution directly, so trait and global parameters both work.
  x <- percent_error_from_post(post_out)
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    plot.new()
    title(main = paste0(treatment_title, "\nerror"))
    text(0.5, 0.5, "no finite errors")
    return(invisible(NULL))
  }

  ## Trim only for display so rare extreme ratios do not destroy the panel.
  ## The vertical dashed lines mark the untrimmed median and mean.
  x_med  <- median(x, na.rm = TRUE)
  x_mean <- mean(x, na.rm = TRUE)

  x_plot <- x
  if (length(x_plot) > 20) {
    q <- quantile(x_plot, probs = c(0.01, 0.99), na.rm = TRUE)
    x_plot <- x_plot[x_plot >= q[[1]] & x_plot <= q[[2]]]
  }

  if (length(unique(x_plot)) < 2) {
    plot(NULL, NULL,
         xlim = range(x_plot) + c(-0.1, 0.1),
         ylim = c(0, 1),
         xlab = "error",
         ylab = "density",
         main = paste0(treatment_title, "\nerror"),
         cex.main = 0.95)
    rug(x_plot)
  } else {
    h <- hist(x_plot, breaks = 30, plot = FALSE)
    plot(h,
         freq = FALSE,
         main = paste0(treatment_title, "\nerror"),
         xlab = "error",
         ylab = "density",
         col = "grey85",
         border = "grey40",
         cex.main = 0.95)
  }

  abline(v = 0, col = "grey", lwd = 1.5)
  abline(v = x_med, col = 3, lwd = 1.5)
  abline(v = x_mean, col = 3, lty = 2, lwd = 1.5)
}


plot_error_density <- function(post_out, treatment_title,
                                bad_cut = 0.25,
                                extreme_cut = 1,
                                breaks = 50) {
  x <- percent_error_from_post(post_out)
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    plot.new()
    title(main = paste0(treatment_title, "\nerror"))
    text(0.5, 0.5, "no finite errors")
    return(invisible(NULL))
  }
  
  n <- length(x)
  
  bad     <- abs(x) > bad_cut
  extreme <- abs(x) > extreme_cut
  
  med <- median(x)
  mn  <- mean(x)
  
  ## Keep all points, but cap plotting window.
  xlim <- quantile(x, c(0.01, 0.99), na.rm = TRUE)
  xlim <- range(c(xlim, -bad_cut, bad_cut, 0))
  
  hist(x,
       breaks = breaks,
       freq = FALSE,
       xlim = xlim,
       main = paste0(treatment_title, "\nerror distribution"),
       xlab = "percent error",
       ylab = "density",
       col = "grey85",
       border = "grey40",
       cex.main = 0.95)
  
  abline(v = 0, col = "grey40", lwd = 1.5)
  abline(v = med, col = 3, lwd = 1.5)
  abline(v = mn, col = 3, lty = 2, lwd = 1.5)
  
  abline(v = c(-bad_cut, bad_cut), col = 2, lty = 2, lwd = 1.25)
  abline(v = c(-extreme_cut, extreme_cut), col = 2, lty = 3, lwd = 1.25)
  
  rug(x[x >= xlim[1] & x <= xlim[2]])
  
  usr <- par("usr")
  text(x = usr[1], y = usr[4],
       labels = paste0(
         "n = ", n,
         "\n|err| > ", bad_cut, ": ", sum(bad), " (", round(mean(bad) * 100, 1), "%)",
         "\n|err| > ", extreme_cut, ": ", sum(extreme), " (", round(mean(extreme) * 100, 1), "%)"
       ),
       adj = c(0, 1),
       cex = 0.8)
  
  invisible(list(
    n = n,
    median = med,
    mean = mn,
    n_bad = sum(bad),
    prop_bad = mean(bad),
    n_extreme = sum(extreme),
    prop_extreme = mean(extreme),
    errors = x
  ))
}


plot_error_density <- function(post_out, treatment_title,
                               bad_cut = 0.25,
                               breaks = 25) {
  
  err_df <- percent_error_from_post(post_out)
  
  ## Expected columns:
  ## err_df$run
  ## err_df$error
  
  err_df <- err_df[is.finite(err_df$error), ]
  
  if (nrow(err_df) == 0) {
    plot.new()
    title(main = paste0(treatment_title, "\nmax error per run"))
    text(0.5, 0.5, "no finite errors")
    return(invisible(NULL))
  }
  
  run_max <- tapply(abs(err_df$error), err_df$run, max, na.rm = TRUE)
  run_max <- run_max[is.finite(run_max)]
  
  n_runs <- length(run_max)
  n_bad_runs <- sum(run_max > bad_cut)
  
  hist(run_max,
       breaks = breaks,
       freq = FALSE,
       col = "grey85",
       border = "grey40",
       main = paste0(treatment_title, "\nmax error per run"),
       xlab = "max |percent error| in run",
       ylab = "density",
       cex.main = 0.95)
  
  abline(v = bad_cut, col = 2, lty = 2, lwd = 1.5)
  abline(v = median(run_max), col = 3, lwd = 1.5)
  rug(run_max)
  
  usr <- par("usr")
  text(usr[1], usr[4],
       labels = paste0(
         "runs = ", n_runs,
         "\nmax |err| > ", bad_cut, ": ",
         n_bad_runs, " (",
         sprintf("%.1f", 100 * n_bad_runs / n_runs), "%)"
       ),
       adj = c(0, 1),
       cex = 0.75)
  
  invisible(run_max)
}


plot_error_density <- function(post_out, treatment_title,
                               n_species = 200,
                               n_predictors = 2,
                               bad_cut = 0.25,
                               breaks = 25) {
  
  x <- percent_error_from_post(post_out)
  x <- x[is.finite(x)]
  
  n_per_run <- n_species * n_predictors
  
  if (length(x) == 0) {
    plot.new()
    title(main = paste0(treatment_title, "\nmax error per run"))
    text(0.5, 0.5, "no finite errors")
    return(invisible(NULL))
  }
  
  n_runs <- floor(length(x) / n_per_run)
  
  if (n_runs < 1) {
    plot.new()
    title(main = paste0(treatment_title, "\nmax error per run"))
    text(0.5, 0.5, "not enough rows for one run")
    return(invisible(NULL))
  }
  
  ## Drop incomplete final run, if present.
  x <- x[seq_len(n_runs * n_per_run)]
  
  run_id <- rep(seq_len(n_runs), each = n_per_run)
  
  run_max <- tapply(abs(x), run_id, max, na.rm = TRUE)
  run_max <- as.numeric(run_max)
  
  bad_runs <- run_max > bad_cut
  
  hist(run_max,
       breaks = breaks,
       freq = FALSE,
       col = "grey85",
       border = "grey40",
       main = paste0(treatment_title, "\nmax error per run"),
       xlab = "max |error| in run",
       ylab = "density",
       cex.main = 0.95)
  
  abline(v = bad_cut, col = 2, lty = 2, lwd = 1.5)
  abline(v = median(run_max), col = 3, lwd = 1.5)
  rug(run_max)
  
  usr <- par("usr")
  text(usr[1], usr[4],
       labels = paste0(
         "runs = ", n_runs,
         "\nrows/run = ", n_per_run,
         "\nmax |err| > ", bad_cut, ": ",
         sum(bad_runs), " (",
         sprintf("%.1f", mean(bad_runs) * 100), "%)"
       ),
       adj = c(0, 1),
       cex = 0.72)
  
  invisible(list(
    n_runs = n_runs,
    n_per_run = n_per_run,
    run_max = run_max,
    bad_runs = bad_runs
  ))
}



plot_error_density <- function(post_out, treatment_title,
                               n_species = 200,
                               n_predictors = 2,
                               breaks = 25) {
  
  x <- percent_error_from_post(post_out)
  x <- x[is.finite(x)]
  
  n_per_run <- n_species * n_predictors
  
  if (length(x) == 0) {
    plot.new()
    title(main = paste0(treatment_title, "\nmax error per run"))
    text(0.5, 0.5, "no finite errors")
    return(invisible(NULL))
  }
  
  n_runs <- floor(length(x) / n_per_run)
  
  if (n_runs < 1) {
    plot.new()
    title(main = paste0(treatment_title, "\nmax error per run"))
    text(0.5, 0.5, "not enough rows for one run")
    return(invisible(NULL))
  }
  
  x <- x[seq_len(n_runs * n_per_run)]
  
  run_id <- rep(seq_len(n_runs), each = n_per_run)
  
  run_max <- tapply(abs(x), run_id, max, na.rm = TRUE)
  run_max <- as.numeric(run_max)
  run_max <- run_max[is.finite(run_max)]
  
  med_run_max <- median(run_max)
  mean_run_max <- mean(run_max)
  q95_run_max <- quantile(run_max, 0.95, na.rm = TRUE)
  worst_run_max <- max(run_max)
  
  hist(run_max,
       breaks = breaks,
       freq = FALSE,
       col = "grey85",
       border = "grey40",
       main = paste0(treatment_title, "\nmax error per run"),
       xlab = "max |error| within run",
       ylab = "density",
       cex.main = 0.95)
  
  abline(v = med_run_max, col = 3, lwd = 1.5)
  abline(v = mean_run_max, col = 3, lty = 2, lwd = 1.5)
  abline(v = q95_run_max, col = 2, lty = 2, lwd = 1.25)
  
  rug(run_max)
  
  usr <- par("usr")
  text(usr[1], usr[4],
       labels = paste0(
         "runs = ", length(run_max),
         "\nmedian = ", round(med_run_max, 3),
         "\nmean = ", round(mean_run_max, 3),
         "\n95% = ", round(q95_run_max, 3),
         "\nworst = ", round(worst_run_max, 3)
       ),
       adj = c(0, 1),
       cex = 0.72)
  
  invisible(list(
    n_runs = length(run_max),
    n_per_run = n_per_run,
    run_max = run_max,
    median = med_run_max,
    mean = mean_run_max,
    q95 = q95_run_max,
    worst = worst_run_max
  ))
}
## Combined line/error treatment panel per pooled parameter.
## Layout is 3 rows x 6 columns:
##   for each mean column: line panel, then matching error panel.
#for (par_name in pooled_order) {
#
#  print(par_name)
#
#  png(filename = file.path(dir, paste0("pooled_", par_name,
#                                      "_line_error_by_treatment.png")),
#      width = 18, height = 12, units = "in", res = 300)
#
#  oldpar <- par(no.readonly = TRUE)
#
#  par(mfrow = c(length(sd_vec), length(height_vec) * 2),
#      mar = c(4.2, 4.2, 3.1, 0.8),
#      oma = c(0, 0, 3.0, 0))
#
#  for (height_sd in sd_vec) {
#    for (height_mean in height_vec) {
#      MCMC_summary_list <- read_treatment_summary(height_mean, height_sd)
#      post_out <- pool_post_out(MCMC_summary_list, par_name)
#      treatment_title <- paste0("mean = ", height_mean, ", sd = ", height_sd)
#
#      plot_post_line(post_out, treatment_title)
#      plot_error_density(post_out, treatment_title)
#    }
#  }
#
#  mtext(pooled_labels[[par_name]], outer = TRUE, line = 0.8, cex = 1.8)
#
#  dev.off()
#  par(oldpar)
#}
#
#cat("Wrote combined line/error treatment-panel PNGs to: ", dir, "\n", sep = "")



## ------------------------------------------------------------
## Plot pooled parameter panels
##
## Trait parameters with enough rows:
##   h, w, c -> line + error panels
##
## mvBM/global parameters:
##   A_C, Rsd_C, A_W, Rsd_W, Rcor -> line-only panels
## ------------------------------------------------------------

trait_pars <- c("h", "w", "c")
mvbm_pars  <- c("A_C", "Rsd_C", "A_W", "Rsd_W", "Rcor")


## Trait parameters: line + error panels
for (par_name in trait_pars) {
  
  print(par_name)
  
  png(filename = file.path(dir, paste0("pooled_", par_name,
                                       "_line_error_by_treatment.png")),
      width = 18, height = 12, units = "in", res = 300)
  
  oldpar <- par(no.readonly = TRUE)
  
  par(mfrow = c(length(sd_vec), length(height_vec) * 2),
      mar = c(4.2, 4.2, 3.1, 0.8),
      oma = c(0, 0, 3.0, 0))
  
  for (height_sd in sd_vec) {
    for (height_mean in height_vec) {
      
      MCMC_summary_list <- read_treatment_summary(height_mean, height_sd)
      post_out <- pool_post_out(MCMC_summary_list, par_name)
      treatment_title <- paste0("mean = ", height_mean, ", sd = ", height_sd)
      
      plot_post_line(post_out, treatment_title)
      plot_error_density(post_out, treatment_title)
    }
  }
  
  mtext(pooled_labels[[par_name]], outer = TRUE, line = 0.8, cex = 1.8)
  
  dev.off()
  par(oldpar)
}


## mvBM/global parameters: line-only panels
for (par_name in mvbm_pars) {
  
  print(par_name)
  
  png(filename = file.path(dir, paste0("pooled_", par_name,
                                       "_line_only_by_treatment.png")),
      width = 12, height = 10, units = "in", res = 300)
  
  oldpar <- par(no.readonly = TRUE)
  
  par(mfrow = c(length(sd_vec), length(height_vec)),
      mar = c(4.2, 4.2, 3.1, 0.8),
      oma = c(0, 0, 3.0, 0))
  
  for (height_sd in sd_vec) {
    for (height_mean in height_vec) {
      
      MCMC_summary_list <- read_treatment_summary(height_mean, height_sd)
      post_out <- pool_post_out(MCMC_summary_list, par_name)
      treatment_title <- paste0("mean = ", height_mean, ", sd = ", height_sd)
      
      plot_post_line(post_out, treatment_title)
    }
  }
  
  mtext(pooled_labels[[par_name]], outer = TRUE, line = 0.8, cex = 1.8)
  
  dev.off()
  par(oldpar)
}

cat("Wrote pooled treatment-panel PNGs to: ", dir, "\n", sep = "")

#####################################

## ------------------------------------------------------------
## Point-level tolerance diagnostics
##
## Writes four PNGs:
##   1. tolerance error by optimum error
##   2. tolerance error by breadth error
##   3. tolerance posterior median by optimum error
##   4. tolerance posterior median by breadth error
##
## Each PNG is a 3 x 3 treatment panel.
## Rows    = height_sd
## Columns = height_mean
##
## Assumes:
##   post_out[, 1] = posterior median
##   post_out[, 2] = true value
## ------------------------------------------------------------

make_point_trait_df <- function(post_out) {
  
  good <- is.finite(post_out[, 1]) & is.finite(post_out[, 2])
  post_out <- post_out[good, , drop = FALSE]
  
  if (nrow(post_out) == 0) {
    return(data.frame())
  }
  
  post_val <- post_out[, 1]
  true_val <- post_out[, 2]
  
  data.frame(
    post_median = post_val,
    abs_error   = abs(post_val - true_val)
  )
}


collect_tolerance_point_diagnostics <- function(height_mean, height_sd) {
  
  MCMC_summary_list <- read_treatment_summary(height_mean, height_sd)
  
  h_df <- make_point_trait_df(pool_post_out(MCMC_summary_list, "h"))
  c_df <- make_point_trait_df(pool_post_out(MCMC_summary_list, "c"))
  w_df <- make_point_trait_df(pool_post_out(MCMC_summary_list, "w"))
  
  n <- min(nrow(h_df), nrow(c_df), nrow(w_df))
  
  if (n == 0) {
    return(data.frame())
  }
  
  data.frame(
    height_mean = height_mean,
    height_sd   = height_sd,
    
    tolerance_abs_error    = h_df$abs_error[seq_len(n)],
    tolerance_post_median  = h_df$post_median[seq_len(n)],
    
    optimum_abs_error      = c_df$abs_error[seq_len(n)],
    breadth_abs_error      = w_df$abs_error[seq_len(n)]
  )
}


all_tol_diag <- do.call(
  rbind,
  lapply(sd_vec, function(height_sd) {
    do.call(
      rbind,
      lapply(height_vec, function(height_mean) {
        collect_tolerance_point_diagnostics(height_mean, height_sd)
      })
    )
  })
)


plot_one_scatter <- function(x, y, xlab, ylab, main,
                             point_cex = 0.18,
                             point_col = make.transparent("black", 1.0)) {
  
  good <- is.finite(x) & is.finite(y)
  x <- x[good]
  y <- y[good]
  
  if (length(x) < 2) {
    plot.new()
    title(main = main)
    text(0.5, 0.5, "not enough data")
    return(invisible(NULL))
  }
  
  xlim <- range(x, finite = TRUE)
  ylim <- range(y, finite = TRUE)
  
  if (diff(xlim) == 0) xlim <- xlim + c(-0.05, 0.05)
  if (diff(ylim) == 0) ylim <- ylim + c(-0.05, 0.05)
  
  ylim[1] <- min(0, ylim[1])
  
  plot(x, y,
       pch = 16,
       cex = point_cex,
       col = point_col,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xlim = xlim,
       ylim = ylim,
       cex.main = 0.9)
  
  if (length(x) > 2) {
    abline(lm(y ~ x), col = 2, lwd = 1.5)
  }
  
  cc <- suppressWarnings(cor(x, y, use = "complete.obs"))
  
  usr <- par("usr")
  text(usr[1], usr[4],
       labels = paste0("cor = ", round(cc, 3)),
       adj = c(0, 1),
       cex = 0.8)
  
  invisible(cc)
}


plot_treatment_grid <- function(filename,
                                y_var,
                                x_var,
                                ylab,
                                xlab,
                                outer_title) {
  
  png(filename = file.path(dir, filename),
      width = 12, height = 10, units = "in", res = 300)
  
  oldpar <- par(no.readonly = TRUE)
  
  par(mfrow = c(length(sd_vec), length(height_vec)),
      mar = c(4.2, 4.2, 3.0, 0.8),
      oma = c(0, 0, 3, 0))
  
  for (height_sd in sd_vec) {
    for (height_mean in height_vec) {
      
      df_now <- all_tol_diag[
        all_tol_diag$height_mean == height_mean &
          all_tol_diag$height_sd == height_sd,
      ]
      
      treatment_title <- paste0("mean = ", height_mean, ", sd = ", height_sd)
      
      plot_one_scatter(
        x = df_now[[x_var]],
        y = df_now[[y_var]],
        xlab = xlab,
        ylab = ylab,
        main = treatment_title
      )
    }
  }
  
  mtext(outer_title, outer = TRUE, line = 1, cex = 1.5)
  
  dev.off()
  par(oldpar)
}


plot_treatment_grid(
  filename = "tolerance_error_by_optimum_error.png",
  y_var = "tolerance_abs_error",
  x_var = "optimum_abs_error",
  ylab = "tolerance |error|",
  xlab = "optimum |error|",
  outer_title = "Tolerance error by optimum error"
)

plot_treatment_grid(
  filename = "tolerance_error_by_breadth_error.png",
  y_var = "tolerance_abs_error",
  x_var = "breadth_abs_error",
  ylab = "tolerance |error|",
  xlab = "breadth |error|",
  outer_title = "Tolerance error by breadth error"
)

plot_treatment_grid(
  filename = "tolerance_postmedian_by_optimum_error.png",
  y_var = "tolerance_post_median",
  x_var = "optimum_abs_error",
  ylab = "tolerance posterior median",
  xlab = "optimum |error|",
  outer_title = "Tolerance posterior median by optimum error"
)

plot_treatment_grid(
  filename = "tolerance_postmedian_by_breadth_error.png",
  y_var = "tolerance_post_median",
  x_var = "breadth_abs_error",
  ylab = "tolerance posterior median",
  xlab = "breadth |error|",
  outer_title = "Tolerance posterior median by breadth error"
)

cat("Wrote four point-level tolerance diagnostic PNGs to: ",
    dir, "\n", sep = "")


plot_treatment_grid <- function(filename_base,
                                y_var,
                                x_var,
                                ylab,
                                xlab,
                                outer_title) {
  
  for (dev_type in c("png", "pdf")) {
    
    if (dev_type == "png") {
      png(filename = file.path(dir, paste0(filename_base, ".png")),
          width = 12, height = 10, units = "in", res = 300)
    }
    
    if (dev_type == "pdf") {
      pdf(file = file.path(dir, paste0(filename_base, ".pdf")),
          width = 12, height = 10)
    }
    
    oldpar <- par(no.readonly = TRUE)
    
    par(mfrow = c(length(sd_vec), length(height_vec)),
        mar = c(4.2, 4.2, 3.0, 0.8),
        oma = c(0, 0, 3, 0))
    
    for (height_sd in sd_vec) {
      for (height_mean in height_vec) {
        
        df_now <- all_tol_diag[
          all_tol_diag$height_mean == height_mean &
            all_tol_diag$height_sd == height_sd,
        ]
        
        treatment_title <- paste0("mean = ", height_mean, ", sd = ", height_sd)
        
        plot_one_scatter(
          x = df_now[[x_var]],
          y = df_now[[y_var]],
          xlab = xlab,
          ylab = ylab,
          main = treatment_title
        )
      }
    }
    
    mtext(outer_title, outer = TRUE, line = 1, cex = 1.5)
    
    dev.off()
    par(oldpar)
  }
}


plot_treatment_grid(
  filename_base = "tolerance_error_by_optimum_error",
  y_var = "tolerance_abs_error",
  x_var = "optimum_abs_error",
  ylab = "tolerance |error|",
  xlab = "optimum |error|",
  outer_title = "Tolerance error by optimum error"
)

plot_treatment_grid(
  filename_base = "tolerance_error_by_breadth_error",
  y_var = "tolerance_abs_error",
  x_var = "breadth_abs_error",
  ylab = "tolerance |error|",
  xlab = "breadth |error|",
  outer_title = "Tolerance error by breadth error"
)

plot_treatment_grid(
  filename_base = "tolerance_postmedian_by_optimum_error",
  y_var = "tolerance_post_median",
  x_var = "optimum_abs_error",
  ylab = "tolerance posterior median",
  xlab = "optimum |error|",
  outer_title = "Tolerance posterior median by optimum error"
)

plot_treatment_grid(
  filename_base = "tolerance_postmedian_by_breadth_error",
  y_var = "tolerance_post_median",
  x_var = "breadth_abs_error",
  ylab = "tolerance posterior median",
  xlab = "breadth |error|",
  outer_title = "Tolerance posterior median by breadth error"
)

cat("Wrote point-level tolerance diagnostic PNGs and PDFs to: ",
    dir, "\n", sep = "")



## ------------------------------------------------------------
## Breadth error vs optimum error
## Point-level pooled across predictors
## Writes PNG and PDF
## ------------------------------------------------------------

{
make_error_df <- function(post_out) {
  
  good <- is.finite(post_out[,1]) & is.finite(post_out[,2])
  post_out <- post_out[good, , drop = FALSE]
  
  data.frame(
    abs_error = abs(post_out[,1] - post_out[,2])
  )
}


## ---------- PNG ----------

png(
  filename = file.path(dir, "breadth_error_by_optimum_error.png"),
  width = 12,
  height = 10,
  units = "in",
  res = 300
)

oldpar <- par(no.readonly = TRUE)

par(mfrow = c(length(sd_vec), length(height_vec)),
    mar = c(4.2, 4.2, 3.0, 0.8),
    oma = c(0, 0, 3, 0))

for (height_sd in sd_vec) {
  for (height_mean in height_vec) {
    
    MCMC_summary_list <- read_treatment_summary(height_mean,
                                                height_sd)
    
    c_df <- make_error_df(
      pool_post_out(MCMC_summary_list, "c")
    )
    
    w_df <- make_error_df(
      pool_post_out(MCMC_summary_list, "w")
    )
    
    n <- min(nrow(c_df), nrow(w_df))
    
    treatment_title <- paste0(
      "mean = ", height_mean,
      ", sd = ", height_sd
    )
    
    plot_one_scatter(
      x = c_df$abs_error[seq_len(n)],
      y = w_df$abs_error[seq_len(n)],
      xlab = "optimum |error|",
      ylab = "breadth |error|",
      main = treatment_title
    )
  }
}

mtext("Breadth error by optimum error",
      outer = TRUE,
      line = 1,
      cex = 1.5)

dev.off()
par(oldpar)


## ---------- PDF ----------

pdf(
  file = file.path(dir, "breadth_error_by_optimum_error.pdf"),
  width = 12,
  height = 10
)

oldpar <- par(no.readonly = TRUE)

par(mfrow = c(length(sd_vec), length(height_vec)),
    mar = c(4.2, 4.2, 3.0, 0.8),
    oma = c(0, 0, 3, 0))

for (height_sd in sd_vec) {
  for (height_mean in height_vec) {
    
    MCMC_summary_list <- read_treatment_summary(height_mean,
                                                height_sd)
    
    c_df <- make_error_df(
      pool_post_out(MCMC_summary_list, "c")
    )
    
    w_df <- make_error_df(
      pool_post_out(MCMC_summary_list, "w")
    )
    
    n <- min(nrow(c_df), nrow(w_df))
    
    treatment_title <- paste0(
      "mean = ", height_mean,
      ", sd = ", height_sd
    )
    
    plot_one_scatter(
      x = c_df$abs_error[seq_len(n)],
      y = w_df$abs_error[seq_len(n)],
      xlab = "optimum |error|",
      ylab = "breadth |error|",
      main = treatment_title
    )
  }
}

mtext("Breadth error by optimum error",
      outer = TRUE,
      line = 1,
      cex = 1.5)

dev.off()
par(oldpar)

cat("Wrote breadth-vs-optimum error PNG and PDF to: ",
    dir, "\n", sep = "")

}
