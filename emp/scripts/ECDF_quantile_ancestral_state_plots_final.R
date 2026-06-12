#library(overlapping)
require(LaplacesDemon)
require(philentropy)
require(phytools)
require(coda)

library(BePhyNE)
phylo <- ENA_Pleth_Tree

new_miss_logdf <- readRDS(
  file = file.path(
    "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_miss_compiled_logs_final",
    paste0("compiled_missing_species_logdf", ".RDS")
  )
)

new_full_logdf <- readRDS(
  "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_full_compiled_logs_final/uninform_height_prior_news_compiled_full_species_logdf.RDS"
)

plot_dir <- "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/plots"

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

  pred_i  <- param_table$pred[param_table$id == param_id]
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

save_ecdf_phylo_plot <- function(pdens, psumcM, psig, file_base, xlab, phylo) {

  for(dev_type in c("pdf", "png")) {

    if(dev_type == "pdf") {
      pdf(
        file.path(plot_dir, paste0(file_base, ".pdf")),
        width = 14,
        height = 10,
        pointsize = 8
      )
    }

    if(dev_type == "png") {
      png(
        file.path(plot_dir, paste0(file_base, ".png")),
        width = 14,
        height = 10,
        units = "in",
        res = 300,
        pointsize = 8
      )
    }

    pp <- par(no.readonly = TRUE)

    layout(matrix(c(1,2), nrow=1), widths=c(0.90,1.35))

    par(mar=c(5,5,2,2))
    hist(
      unlist(pdens),
      main = "",
      xlab = xlab
    )
    abline(v = median(sapply(pdens, median)))

    par(mar=c(5,1,2,7), xpd=NA)
    plot(psumcM)

    text(
      x = 0.83 * max(nodeHeights(phylo)),
      y = -4.8,
      labels = paste("Blomberg's K = ", round(psig$K, 2), ", p = ", round(psig$P, 4), sep = ""),
      adj = c(0,0.5),
      cex = 0.9
    )

    par(pp)
    dev.off()
  }
}

save_true_pred_plot <- function(pred_median, data_median, file_base) {

  for(dev_type in c("pdf", "png")) {

    if(dev_type == "pdf") {
      pdf(
        file.path(plot_dir, paste0(file_base, ".pdf")),
        width = 10,
        height = 10,
        pointsize = 11
      )
    }

    if(dev_type == "png") {
      png(
        file.path(plot_dir, paste0(file_base, ".png")),
        width = 10,
        height = 10,
        units = "in",
        res = 300,
        pointsize = 11
      )
    }

    pp <- par(no.readonly = TRUE)

    par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

    lims <- list(c(-2, 3), c(-2, 3), c(0, 4), c(0, 4))

    true_pred_labels <- list(
      expression(theta[P]),
      expression(theta[T]),
      expression(omega[P]),
      expression(omega[T])
    )

    for(i in 1:4) {
      lm1 <- lm(pred_median[[i]] ~ data_median[[i]])
      slm1 <- summary(lm1)

      plot(
        data_median[[i]],
        pred_median[[i]],
        main = true_pred_labels[[i]],
        xlim = lims[[i]],
        ylim = lims[[i]],
        xlab = bquote("Data-estimated " ~ .(true_pred_labels[[i]])),
        ylab = bquote("Missing-species predicted " ~ .(true_pred_labels[[i]])),
        pch = 21,
        bg = "gray80"
      )

      text(
        diff(lims[[i]]) * 0.08 + min(lims[[i]]),
        diff(lims[[i]]) * 0.92 + min(lims[[i]]),
        paste(
          "Slope = ",
          signif(slm1$coef[2, 1], 2),
          ", p = ",
          signif(slm1$coef[2, 4], 2),
          sep = ""
        ),
        pos = 4
      )

      abline(0, 1, lty = 2, col = "gray50")
      abline(lm1, col = "red", lwd = 2)
    }

    par(pp)
    dev.off()
  }
}

get_ecdf_phylo_objects <- function(traits, xlab_expr, file_base) {

  dist_same <- list()
  dist_other <- list()

  for(j in 1:nspp) {

    pred_post <- do.call(
      cbind,
      lapply(traits, function(x) get_species_param_trace(new_miss_logdf, x, j, phylo, param_table))
    )

    pred_center <- apply(pred_post, 2, mean)
    Sigma <- var(pred_post)

    dist_other[[j]] <- list()

    for(i in 1:nspp) {

      post <- do.call(
        cbind,
        lapply(traits, function(x) get_species_param_trace(new_full_logdf, x, i, phylo, param_table))
      )

      if(i == j) {
        dist_same[[j]] <- mahalanobis(post, pred_center, Sigma)
      } else {
        dist_other[[j]][[i]] <- mahalanobis(post, pred_center, Sigma)
      }
    }
  }

  dist_other <- lapply(dist_other, unlist)

  pdens <- list()

  for(i in 1:nspp) {
    ecdens <- ecdf(dist_other[[i]])
    pdens[[i]] <- ecdens(dist_same[[i]])
  }

  D_i <- sapply(dist_same, median)

  psum <- sapply(pdens, median)
  psig <- phylosig(phylo, setNames(psum, phylo$tip.label), "K", test = TRUE)

  psumcM <- contMap(
    phylo,
    setNames(psum, phylo$tip.label),
    fsize = 0.42,
    plot = FALSE
  )

  psumcM <- setMap(psumcM, cm.colors(10))

  save_ecdf_phylo_plot(
    pdens = pdens,
    psumcM = psumcM,
    psig = psig,
    file_base = file_base,
    xlab = xlab_expr,
    phylo = phylo
  )

  invisible(
    list(
      pdens = pdens,
      D_i = D_i,
      psum = psum,
      psig = psig,
      psumcM = psumcM,
      dist_same = dist_same,
      dist_other = dist_other
    )
  )
}

npars <- nrow(param_table)
nspp  <- length(phylo$tip.label)

PhyNE_miss_median_parlist <- lapply(1:npars, function(par) {
  unlist(lapply(1:nspp, function(sp) {
    median(get_species_param_trace(new_miss_logdf, par, sp, phylo, param_table), na.rm = TRUE)
  }))
})

PhyNE_miss_medians <- list(
  cbind(
    PhyNE_miss_median_parlist[[1]],
    PhyNE_miss_median_parlist[[2]],
    PhyNE_miss_median_parlist[[3]]
  ),
  cbind(
    PhyNE_miss_median_parlist[[4]],
    PhyNE_miss_median_parlist[[5]],
    PhyNE_miss_median_parlist[[6]]
  )
)

logitTransform <- function(p) {
  log(p / (1 - p))
}

ECDF_miss_sp <- list()
ECDF_fit_samesp <- list()

for(sp in 1:length(phylo$tip.label)) {

  reps <- length(get_species_param_trace(new_miss_logdf, 1, sp, phylo, param_table))

  print(paste("sp", sp))

  ECDF_miss_sp[[sp]] <- list()
  ECDF_fit_samesp[[sp]] <- list()

  if(reps == 0) {
    print(paste(sp, "no runs"))
    next
  }

  for(par in 1:npars) {

    print(paste("par", par))

    ECDF_miss_sp[[sp]][[par]] <- ecdf(
      get_species_param_trace(new_miss_logdf, par, sp, phylo, param_table)
    )

    ECDF_fit_samesp[[sp]][[par]] <- ECDF_miss_sp[[sp]][[par]](
      get_species_param_trace(new_full_logdf, par, sp, phylo, param_table)
    )
  }
}

lower_cut <- 0.05
upper_cut <- 0.95

props_samesp <- lapply(
  1:length(ECDF_fit_samesp),
  function(i) sapply(ECDF_fit_samesp[[i]], function(x) mean(x > lower_cut & x < upper_cut))
)

median_samesp <- lapply(
  1:length(ECDF_fit_samesp),
  function(i) sapply(ECDF_fit_samesp[[i]], function(x) median(x))
)

props_samesp <- do.call(rbind, props_samesp)
median_samesp <- do.call(rbind, median_samesp)

ECDF_fit_diffsp_many_green <- list()

for(sp in 1:length(phylo$tip.label)) {

  reps <- length(get_species_param_trace(new_miss_logdf, 1, sp, phylo, param_table))

  print(paste("sp", sp))

  ECDF_fit_diffsp_many_green[[sp]] <- list()

  if(reps == 0) {
    print(paste(sp, "no runs"))
    next
  }

  for(par in 1:npars) {

    print(paste("par", par))

    ECDF_fit_diffsp_many_green[[sp]][[par]] <- lapply(
      (1:length(phylo$tip.label))[(1:length(phylo$tip.label)) != sp],
      function(x) {
        ECDF_miss_sp[[sp]][[par]](
          get_species_param_trace(new_full_logdf, par, x, phylo, param_table)
        )
      }
    )
  }
}

props_diffsp <- lapply(
  1:length(ECDF_fit_diffsp_many_green),
  function(i) sapply(
    ECDF_fit_diffsp_many_green[[i]],
    function(x) sapply(x, function(j) mean(j > lower_cut & j < upper_cut))
  )
)

all_param_res <- get_ecdf_phylo_objects(
  traits = 1:4,
  xlab_expr = expression("Quantile of Missing vs. Data Distance against Null: " ~ theta[P] ~ "," ~ theta[T] ~ "," ~ omega[P] ~ "," ~ omega[T]),
  file_base = "ECDF_quantile_ancestral_state_all_params"
)

pred_P_res <- get_ecdf_phylo_objects(
  traits = c(1, 3),
  xlab_expr = expression("Quantile of Missing vs. Data Distance against Null: " ~ theta[P] ~ " and " ~ omega[P]),
  file_base = "ECDF_quantile_ancestral_state_predictor_P"
)

pred_T_res <- get_ecdf_phylo_objects(
  traits = c(2, 4),
  xlab_expr = expression("Quantile of Missing vs. Data Distance against Null: " ~ theta[T] ~ " and " ~ omega[T]),
  file_base = "ECDF_quantile_ancestral_state_predictor_T"
)

pred_median <- list(numeric(0))
data_median <- list(numeric(0))

for(par in 1:4) {

  pred_median[[par]] <- numeric(0)
  data_median[[par]] <- numeric(0)

  for(j in 1:nspp) {

    pred_median[[par]][j] <- median(
      get_species_param_trace(new_miss_logdf, par, j, phylo, param_table),
      na.rm = TRUE
    )

    data_median[[par]][j] <- median(
      get_species_param_trace(new_full_logdf, par, j, phylo, param_table),
      na.rm = TRUE
    )
  }
}

save_true_pred_plot(
  pred_median = pred_median,
  data_median = data_median,
  file_base = "true_vs_predicted_parameter_medians"
)

