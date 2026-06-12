#library(overlapping)
#require(LaplacesDemon)
#require(philentropy)
require(phytools)
require(coda)

library(BePhyNE)

sp_col    <- "species"
occ_col   <- "PA"
env_preds <- c("bio12", "bio1")
Npred <- length(env_preds)

phylo <- ENA_Pleth_Tree
tree  <- ENA_Pleth_Tree
Ntips <- length(tree$tip.label)

new_miss_logdf <- readRDS(
  file = file.path(
    "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_miss_compiled_logs",
    paste0("compiled_missing_species_logdf", ".RDS")
  )
)

new_full_logdf <- readRDS(
  "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_full_compiled_logs/compiled_full_species_logdf.RDS"
)

setwd("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp")

sets_full <- readRDS("data/sets_full.RDS")

pres_data_scaled <- readRDS("data/scaled_GBIF_clim_pres.RDS")

scale_atr <- list(
  scale = c(
    attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled) == "bio12"],
    attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled) == "bio1"]
  ),
  center = c(
    attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled) == "bio12"],
    attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled) == "bio1"]
  )
)

new_full_logdf <- readRDS(
  "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_full_compiled_logs_final/uninform_height_prior_news_compiled_full_species_logdf.RDS"
)

new_miss_logdf <- readRDS(
  "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_miss_compiled_logs_final/compiled_missing_species_logdf.RDS"
)

GLM_only_ml <- suppressWarnings(
  BePhyNE:::MLglmStartpars_general(
    species_data = sets_full$training,
    tree = tree
  )
)

GLM_only_ml$start_pars_bt[[1]][, 3][GLM_only_ml$start_pars_bt[[1]][, 3] > 0.99] <- 0.99
GLM_only_ml$start_pars_bt[[2]][, 3][GLM_only_ml$start_pars_bt[[2]][, 3] > 0.99] <- 0.99

full_log_summary <- summarize_logdf(new_full_logdf, scale_atr = NA)
miss_log_summary <- summarize_logdf(new_miss_logdf, scale_atr = NA)

full_log_summary$median_parlist$traits[[1]][[1]] <- full_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label, ]
full_log_summary$median_parlist$traits[[2]][[1]] <- full_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label, ]
miss_log_summary$median_parlist$traits[[1]][[1]] <- miss_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label, ]
miss_log_summary$median_parlist$traits[[2]][[1]] <- miss_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label, ]

AUC_full_list <- AUC_posterior_median(
  full_log_summary,
  sets_full$predicting
)

AUC_full <- unlist(lapply(AUC_full_list, function(i) i$auc))

AUC_miss_list <- AUC_posterior_median(
  miss_log_summary,
  sets_full$predicting
)

AUC_miss <- unlist(lapply(AUC_miss_list, function(i) i$auc))

AUC_glm_list <- BePhyNE:::predict_stats(
  traits = GLM_only_ml$start_pars_bt,
  pa_data = sets_full$predicting,
  plot = FALSE
)

AUC_glm <- unlist(lapply(AUC_glm_list, function(i) i$auc))

AUC_df <- cbind(
  full = AUC_full,
  miss = AUC_miss,
  glm  = AUC_glm
) * 100

rownames(AUC_df) <- tree$tip.label

if (is.data.frame(AUC_df)) {
  AUC_df <- as.matrix(AUC_df)
}

if (is.vector(AUC_df)) {
  AUC_df <- matrix(AUC_df, ncol = 1)
  rownames(AUC_df) <- names(AUC_df[, 1])
  colnames(AUC_df) <- "AUC"
}

if (is.null(rownames(AUC_df))) {
  stop("AUC_df must have rownames matching phylo$tip.label.")
}

if (!all(phylo$tip.label %in% rownames(AUC_df))) {
  stop("Some phylo$tip.label values are missing from rownames(AUC_df).")
}

AUC_df <- AUC_df[phylo$tip.label, , drop = FALSE]

if (max(AUC_df, na.rm = TRUE) <= 1) {
  AUC_df <- AUC_df * 100
}

if (is.null(colnames(AUC_df))) {
  colnames(AUC_df) <- paste0("AUC_", seq_len(ncol(AUC_df)))
}

## Shared AUC scale for all columns.
auc_range <- c(50, 100)
auc_palette <- colorRampPalette(c("white", "#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D"))(100)

param_table <- data.frame(
  id = 1:6,
  label = c("theta_P", "theta_T", "omega_P", "omega_T", "tolerance_P", "tolerance_T"),
  pred = c("pred_1", "pred_2", "pred_1", "pred_2", "pred_1", "pred_2"),
  trait = c("opt", "opt", "brdth", "brdth", "tol", "tol"),
  stringsAsFactors = FALSE
)

param_plot_labels <- list(
  expression(theta[P]),
  expression(theta[T]),
  expression(omega[P]),
  expression(omega[T]),
  expression(tolerance[P]),
  expression(tolerance[T])
)

get_species_param_trace <- function(logdf, param_id, sp, phylo, param_table) {
  if (is.numeric(sp)) {
    sp_name <- phylo$tip.label[sp]
  } else {
    sp_name <- sp
  }

  pred_i <- param_table$pred[param_table$id == param_id]
  trait_i <- param_table$trait[param_table$id == param_id]

  col_i <- paste0(pred_i, "_dat.", trait_i, "_", sp_name)

  if (!col_i %in% colnames(logdf)) {
    stop(paste("Missing column:", col_i))
  }

  x <- as.numeric(logdf[[col_i]])

  if (trait_i == "brdth") {
    x <- exp(x)
  }

  x
}

npars <- nrow(param_table)
nspp <- length(phylo$tip.label)

lower_cut <- 0.05
upper_cut <- 0.95

panel_traits <- list(
  all_traits = 1:4,
  precip = c(1, 3),
  temp = c(2, 4)
)

panel_xlabs <- list(
  all_traits = "Quantile of Missing vs. Data Distance against Null",
  precip = expression("Quantile of Missing vs. Data Distance against Null: " ~ theta[P] ~ " and " ~ omega[P]),
  temp = expression("Quantile of Missing vs. Data Distance against Null: " ~ theta[T] ~ " and " ~ omega[T])
)

panel_outfiles <- c(
  all_traits = "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/plots/quantile_contmap_all_traits_AUC_tiles",
  precip = "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/plots/quantile_contmap_precip_AUC_tiles",
  temp = "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/plots/quantile_contmap_temp_AUC_tiles"
)

for (panel_name in names(panel_traits)) {

  traits <- panel_traits[[panel_name]]

  dist_same <- list()
  dist_other <- list()

  for (j in 1:nspp) {

    pred_post <- do.call(
      cbind,
      lapply(traits, function(x) {
        get_species_param_trace(new_miss_logdf, x, j, phylo, param_table)
      })
    )

    pred_center <- apply(pred_post, 2, mean)
    Sigma <- var(pred_post)

    dist_other[[j]] <- list()

    for (i in 1:nspp) {

      post <- do.call(
        cbind,
        lapply(traits, function(x) {
          get_species_param_trace(new_full_logdf, x, i, phylo, param_table)
        })
      )

      if (i == j) {
        dist_same[[j]] <- mahalanobis(post, pred_center, Sigma)
      } else {
        dist_other[[j]][[i]] <- mahalanobis(post, pred_center, Sigma)
      }
    }
  }

  dist_other <- lapply(dist_other, unlist)

  pdens <- list()

  for (i in 1:nspp) {
    ecdens <- ecdf(dist_other[[i]])
    pdens[[i]] <- ecdens(dist_same[[i]])
  }

  psum <- sapply(pdens, median)

  psig <- phylosig(
    phylo,
    setNames(psum, phylo$tip.label),
    "K",
    test = TRUE
  )

  psumcM <- contMap(
    phylo,
    setNames(psum, phylo$tip.label),
    fsize = 0,
    plot = FALSE
  )

  psumcM <- setMap(psumcM, cm.colors(10))

  for (dev_type in c("pdf", "png")) {

    if (dev_type == "pdf") {
      pdf(
        paste0(panel_outfiles[[panel_name]], ".pdf"),
        height = 16.5,
        width = 8.5
      )
    }

    if (dev_type == "png") {
      png(
        paste0(panel_outfiles[[panel_name]], ".png"),
        height = 16.5,
        width = 8.5,
        units = "in",
        res = 300
      )
    }

    par(
      fig = c(0, 1, 0.805, 1),
      mar = c(3.2, 4.8, 0.05, 1.0),
      new = FALSE
    )

    hist(
      unlist(pdens),
      main = "",
      xlab = "",
      ylab = "Frequency",
      col = "white",
      border = "gray55",
      cex.axis = 0.85,
      cex.lab = 0.9,
      las = 1
    )

    mtext(
      panel_xlabs[[panel_name]],
      side = 1,
      line = 1.8,
      cex = 0.95
    )

    abline(
      v = median(sapply(pdens, median), na.rm = TRUE),
      lwd = 2,
      lty = 2
    )

    ## Large right margin is intentional: it gives space for AUC tiles and species names.
    par(
      fig = c(0, 1, 0.000, 0.805),
      mar = c(2.4, 0.8, 0, 12.2),
      new = TRUE
    )

    tree_height <- max(nodeHeights(phylo))

    plot(
      psumcM,
      fsize = 0.85,
      offset = 4.2,
      legend = TRUE,
      lwd = 2,
      xlim = c(0, tree_height * 2.18)
      #ylim = c(-5, nspp + 2)
    )

    lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

    tip_x <- lp$xx[1:nspp]
    tip_y <- lp$yy[1:nspp]

    auc_scaled <- (AUC_df - auc_range[1]) / diff(auc_range)
    auc_scaled[auc_scaled < 0] <- 0
    auc_scaled[auc_scaled > 1] <- 1

    auc_col_idx <- round(auc_scaled * (length(auc_palette) - 1)) + 1

    auc_cols <- matrix(
      auc_palette[auc_col_idx],
      nrow = nrow(AUC_df),
      ncol = ncol(AUC_df),
      dimnames = dimnames(AUC_df)
    )

    ## Keep tiles close to the tips, but leave clear space before labels.
    auc_x_start <- max(tip_x, na.rm = TRUE) + tree_height * 0.050
    auc_x_step  <- tree_height * 0.065
    species_label_x <- auc_x_start + (ncol(AUC_df) - 1) * auc_x_step + tree_height * 0.160

    for (j in seq_len(ncol(AUC_df))) {
      points(
        x = rep(auc_x_start + (j - 1) * auc_x_step, nspp),
        y = tip_y,
        pch = 22,
        bg = auc_cols[, j],
        col = "gray20",
        cex = 1.00,
        lwd = 0.40,
        xpd = NA
      )
    }

    ## Column names above the heat-tile strip.
    ## Rotated 90 degrees, anchored downward so they do not clip into the histogram panel.
    ## Species labels are left to plot.contMap()/plot.phylo through offset,
    ## not manually redrawn here.
    ## AUC column names below the heat-tile strip.
    ## Labels start just below the lowest AUC squares and extend downward.
    ## The side nearest the squares is aligned to the same y coordinate.
    auc_label_y <- min(tip_y, na.rm = TRUE) - 0.55

    text(
      x = auc_x_start + (seq_len(ncol(AUC_df)) - 1) * auc_x_step,
      y = rep(auc_label_y, ncol(AUC_df)),
      labels = c("BePhyNE", "BePhyNE no data", "GLM only"),
      srt = 270,
      adj = c(0, 0.5),
      cex = 0.50,
      font = 2,
      xpd = NA
    )

    ## Shared AUC heat scale.
    ## Draw this farther right, under the Blomberg's K label.
    auc_legend_x0 <- tree_height * 1.45
    auc_legend_x1 <- tree_height * 1.90 
    auc_legend_y0 <- min(tip_y, na.rm = TRUE) - 4.70 -2
    auc_legend_y1 <- min(tip_y, na.rm = TRUE) - 4.20-2

    auc_legend_n <- length(auc_palette)
    auc_legend_x <- seq(auc_legend_x0, auc_legend_x1, length.out = auc_legend_n + 1)

    for (auc_k in seq_len(auc_legend_n)) {
      rect(
        xleft = auc_legend_x[auc_k],
        ybottom = auc_legend_y0,
        xright = auc_legend_x[auc_k + 1],
        ytop = auc_legend_y1,
        col = auc_palette[auc_k],
        border = NA,
        xpd = NA
      )
    }

    rect(
      xleft = auc_legend_x0,
      ybottom = auc_legend_y0,
      xright = auc_legend_x1,
      ytop = auc_legend_y1,
      border = "gray30",
      lwd = 0.5,
      xpd = NA
    )

    text(
      x = c(auc_legend_x0, auc_legend_x1),
      y = auc_legend_y0 - 0.45,
      labels = auc_range,
      cex = 0.55,
      adj = c(0.5, 1),
      xpd = NA
    )

    text(
      x = mean(c(auc_legend_x0, auc_legend_x1)),
      y = auc_legend_y0 - 1.15,
      labels = "AUC",
      cex = 0.65,
      font = 2,
      xpd = NA
    )

    ## Move Blomberg K label to the lower right, away from the contMap legend.
    text(
      x = tree_height * 1.45,
      y = -1.0-2,
      labels = paste(
        "Blomberg's K = ",
        round(psig$K, 2),
        ", p = ",
        round(psig$P, 4),
        sep = ""
      ),
      adj = 0,
      cex = 0.95,
      xpd = NA
    )

    dev.off()
  }
}
