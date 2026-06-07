## Pooled-across-predictor treatment panels
## Assumes MCMC_summary_list objects have already been processed and saved.
## Produces one combined line/error PDF per pooled parameter.
## For each parameter:
##   rows = sd_vec
##   columns = mean values, with line plot then error-density plot for each mean
##   pooled over predictor 1 and predictor 2

height_vec <- c("0.5", "0.7", "0.95")
sd_vec     <- c("0.15", "0.5", "1")

n_per <- 200
dir <- "New_Height_runs_plots"
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
    paste0("New_Height_runs/Height_runs_mean_", height_mean_label, "_sd_", height_sd_label,
           "_MCMC_summary_list_height.rds")
    #paste0("MCMC_summary_list_height_mean_", height_mean_label,
    #       "_sd_", height_sd_label, ".rds")
  )
  
  hit <- candidates[file.exists(candidates)]
  hit
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
pooled_order <- c("Rcor", "h", "w", "c", "A_C", "Rsd_C", "A_W", "Rsd_W")

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
          #xlab = "percent error",
         ylab = "density",
         main = paste0(treatment_title, "\nerror"),
         cex.main = 0.95)
    rug(x_plot)
  } else {
    h <- hist(x_plot, breaks = 30, plot = FALSE)
    plot(h,
         freq = FALSE,
         main = paste0(treatment_title, "\nerror"),
         #xlab = "percent error",
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

## Combined line/error treatment panel per pooled parameter.
## Layout is 3 rows x 6 columns:
##   for each mean column: line panel, then matching error panel.
for (par_name in pooled_order) {

  print(par_name)

  pdf(file = file.path(dir, paste0("pooled_", par_name,
                                   "_line_error_by_treatment.pdf")),
      width = 18, height = 12)

  oldpar <- par(no.readonly = TRUE)

  par(mfrow = c(length(sd_vec), length(height_vec) * 2),
      mar = c(4.2, 4.2, 3.1, 0.8),
      oma = c(0, 0, 3.0, 0))

  for (height_sd in sd_vec) {
    for (height_mean in height_vec) {
      MCMC_summary_list <- read_treatment_summary(height_mean, height_sd)
      post_out <- pool_post_out(MCMC_summary = MCMC_summary_list, par_name)
      treatment_title <- paste0("mean = ", height_mean, ", sd = ", height_sd)

      plot_post_line(post_out, treatment_title)
      plot_error_density(post_out, treatment_title)
    }
  }

  mtext(pooled_labels[[par_name]], outer = TRUE, line = 0.8, cex = 1.8)

  dev.off()
  par(oldpar)
}
