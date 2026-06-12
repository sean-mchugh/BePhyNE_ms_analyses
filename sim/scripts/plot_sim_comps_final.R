library(vioplot)

## ------------------------------------------------------------
## Labels and parameters
## ------------------------------------------------------------

label_name <- list(
  expression(theta),
  expression(omega),
  expression(A[theta]),
  expression(A[omega]),
  expression(sigma[theta]),
  expression(sigma[omega]),
  expression(R[COR])
)

par_names <- c(
  "full_traits_c",
  "full_traits_w",
  "full_traits_A_C",
  "full_traits_A_W",
  "full_traits_Rsd_C",
  "full_traits_Rsd_W",
  "full_traits_Rcor"
)

## ------------------------------------------------------------
## Helpers
## ------------------------------------------------------------

open_plot_device <- function(file_base, dev_type,
                             width = 10,
                             height = 10,
                             res = 300) {
  
  if (dev_type == "png") {
    png(
      filename = paste0(file_base, ".png"),
      width = width,
      height = height,
      units = "in",
      res = res
    )
  } else if (dev_type == "pdf") {
    pdf(
      file = paste0(file_base, ".pdf"),
      width = width,
      height = height
    )
  } else {
    stop("Unknown dev_type: ", dev_type)
  }
}

get_sim_true_vardifs <- function(simvstrue, treatments, par_names) {
  
  var_diff_df_list <- list()
  
  for (i in seq_along(par_names)) {
    
    var_diffs <- lapply(seq_along(treatments), function(treat) {
      simvstrue[[treat]][[par_names[[i]]]][[1]][, 1] -
        simvstrue[[treat]][[par_names[[i]]]][[1]][, 2]
    })
    
    min_length <- min(unlist(lapply(var_diffs, length)))
    
    var_diff_df <- do.call(
      cbind,
      lapply(var_diffs, function(x) x[sample(seq_along(x), min_length, replace = FALSE)])
    )
    
    colnames(var_diff_df) <- treatments
    
    var_diff_df_list[[i]] <- var_diff_df
  }
  
  names(var_diff_df_list) <- par_names
  
  var_diff_df_list
}

plot_vardiff_boxplots <- function(var_diff_df_list,
                                  file_base,
                                  label_name,
                                  cols) {
  
  for (dev_type in c("png", "pdf")) {
    
    open_plot_device(
      file_base = file_base,
      dev_type = dev_type,
      width = 10,
      height = 10
    )
    
    oldpar <- par(no.readonly = TRUE)
    
    par(
      mfrow = c(4, 2),
      mar = c(2.5, 2.5, 3, 1),
      oma = c(0, 0, 0, 0),
      mgp = c(1.5, 0.5, 0),
      cex.main = 1.5
    )
    
    for (par_i in seq_along(label_name)) {
      
      ymax <- max(abs(var_diff_df_list[[par_i]]), na.rm = TRUE)
      ylim <- c(-ymax, ymax)
      
      boxplot(
        var_diff_df_list[[par_i]],
        ylim = ylim,
        col = cols,
        main = label_name[[par_i]],
        cex.main = 2
      )
      
      abline(h = 0, col = "grey", lty = 2)
    }
    
    plot.new()
    
    dev.off()
    par(oldpar)
  }
}

plot_sim_true_lines <- function(simvstrue,
                                file_base,
                                label_name,
                                cols,
                                legend_labels,
                                legend_title) {
  
  for (dev_type in c("png", "pdf")) {
    
    open_plot_device(
      file_base = file_base,
      dev_type = dev_type,
      width = 10,
      height = 10
    )
    
    oldpar <- par(no.readonly = TRUE)
    
    par(
      mfrow = c(4, 2),
      mar = c(2.5, 2.5, 3, 1),
      oma = c(0, 0, 0, 0),
      mgp = c(1.5, 0.5, 0),
      cex.main = 1.5
    )
    
    for (par_i in seq_along(label_name)) {
      
      for (i in seq_along(simvstrue)) {
        
        par_simvstrue <- simvstrue[[i]][[par_i]][[1]]
        
        if (i == 1) {
          
          plot(
            par_simvstrue[, 2],
            par_simvstrue[, 1],
            ylab = "True",
            xlab = "Posterior Median",
            main = label_name[[par_i]],
            col = cols[[i]],
            cex = 0.3,
            pch = 19
          )
          
        } else {
          
          points(
            par_simvstrue[, 2],
            par_simvstrue[, 1],
            col = cols[[i]],
            cex = 0.2,
            pch = 19
          )
        }
        
        abline(
          lm(par_simvstrue[, 1] ~ par_simvstrue[, 2]),
          col = cols[[i]]
        )
        
        abline(a = 0, b = 1, lwd = 0.5, col = "grey")
      }
    }
    
    plot.new()
    par(xpd = NA)
    
    legend(
      "left",
      col = c(cols, "grey"),
      legend = c(legend_labels, "1:1 trend"),
      lwd = 5,
      cex = 1.4,
      horiz = FALSE,
      title = legend_title,
      ncol = 2
    )
    
    dev.off()
    par(oldpar)
  }
}

## ------------------------------------------------------------
## Kappa plots
## ------------------------------------------------------------

Kappa_simvstrue <- readRDS("~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/outfiles/Kappa_simmedian_vs_true_pars.RDS")

kappa_pars <- paste0(c("0.0", "0.25", "0.5", "0.75", "1.0"))

kappa_var_diff_df_list <- get_sim_true_vardifs(
  simvstrue = Kappa_simvstrue,
  treatments = kappa_pars,
  par_names = par_names
)

plot_vardiff_boxplots(
  var_diff_df_list = kappa_var_diff_df_list,
  file_base = "~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/plots/Kappa_boxplots",
  label_name = label_name,
  cols = 2:6
)

#plot_sim_true_lines(
#  simvstrue = Kappa_simvstrue,
#  file_base = "~/Projects/BePhyNE/BePhyNE_ms_analyses/Kappa_plots",
#  label_name = label_name,
#  cols = 2:6,
#  legend_labels = c(0, 0.25, 0.5, 0.75, 1.0),
#  legend_title = expression(kappa~"Transformation")
#)

## ------------------------------------------------------------
## Lambda plots
## ------------------------------------------------------------

lambda_simvstrue <- readRDS("~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/outfiles/Lambda_simmedian_vs_true_pars.RDS")

lambda_pars <- paste0(c("1.0", "0.0"))

lambda_var_diff_df_list <- get_sim_true_vardifs(
  simvstrue = lambda_simvstrue,
  treatments = lambda_pars,
  par_names = par_names
)

plot_vardiff_boxplots(
  var_diff_df_list = lambda_var_diff_df_list,
  file_base = "~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/plots/Lambda_boxplots",
  label_name = label_name,
  cols = 6:7
)

#plot_sim_true_lines(
#  simvstrue = lambda_simvstrue,
#  file_base = "~/Downloads/Lambda_lines",
#  label_name = label_name,
#  cols = 6:7,
#  legend_labels = c(1.0, 0.0),
#  legend_title = expression(lambda~"Transformation")
#)

## ------------------------------------------------------------
## Background / true absence plots
## ------------------------------------------------------------

back_simvstrue <- readRDS("~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/outfiles/grid_Background_simmedian_vs_true_pars.RDS")

background_pars <- paste0(c("True Absence", "Background"))

back_var_diff_df_list <- get_sim_true_vardifs(
  simvstrue = back_simvstrue,
  treatments = background_pars,
  par_names = par_names
)

plot_vardiff_boxplots(
  var_diff_df_list = back_var_diff_df_list,
  file_base = "~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/plots/Background_boxplots",
  label_name = label_name,
  cols = 6:7
)

#plot_sim_true_lines(
#  simvstrue = back_simvstrue,
#  file_base = "~/Downloads/grid_back_lines",
#  label_name = label_name,
#  cols = 6:7,
#  legend_labels = c("True Absence", "Background"),
#  legend_title = "Absence Data Type"
#)

cat("Wrote PNG and PDF versions of all boxplot and line-plot figures.\n")