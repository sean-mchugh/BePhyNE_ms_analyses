## Combined correlation histogram + representative posterior chain plots by treatment
## Writes PNG and PDF
## Reads RDS files from Neww_New_New_Height_runs
##
## Histograms show all within-chain species x predictor x replicate correlations.
## Line/scatter panels show the single representative species/predictor/replicate
## whose correlation is closest to the median correlation for that treatment.

height_vec <- c("0.5", "0.7", "0.95")
sd_vec     <- c("0.15", "0.5", "1")

indir  <- "~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/outfiles/"
outdir <- "~/Projects/BePhyNE/BePhyNE_ms_analyses/sim/plots/"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

set.seed(1)

label_to_num <- function(x) as.numeric(gsub("p", ".", x, fixed = TRUE))
num_to_label <- function(x) gsub(".", "p", as.character(x), fixed = TRUE)

open_plot_device <- function(filename_base, dev_type,
                             width, height, res = 300) {
  if (dev_type == "png") {
    png(filename = file.path(outdir, paste0(filename_base, ".png")),
        width = width,
        height = height,
        units = "in",
        res = res)
  } else if (dev_type == "pdf") {
    pdf(file = file.path(outdir, paste0(filename_base, ".pdf")),
        width = width,
        height = height)
  } else {
    stop("Unknown dev_type: ", dev_type)
  }
}

parse_treatment_from_file <- function(file) {
  b <- basename(file)
  b <- sub("_trait_height_posterior_representative_samples\\.rds$", "", b)
  b <- sub("_trait_height_posterior_samples_downsampled\\.rds$", "", b)
  b <- sub("_trait_height_posterior_samples\\.rds$", "", b)
  b <- sub("_trait_height_posterior_samples_corr_summary\\.rds$", "", b)
  b <- sub("_trait_height_posterior_correlations\\.rds$", "", b)
  b <- sub("_trait_height_posterior_correlation_summary\\.rds$", "", b)
  b <- sub("_trait_height_median_correlation_representatives\\.rds$", "", b)

  parts <- strsplit(b, "_", fixed = TRUE)[[1]]

  if (length(parts) < 2) {
    stop("Could not parse treatment from filename: ", file)
  }

  list(
    mean_label = parts[[1]],
    sd_label   = parts[[2]],
    mean       = label_to_num(parts[[1]]),
    sd         = label_to_num(parts[[2]])
  )
}

make_treatment_key <- function(mean_label, sd_label) {
  paste(mean_label, sd_label, sep = "__")
}

rds_files <- list.files(
  path = indir,
  pattern = "_trait_height_.*\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No trait-height RDS files found in: ", indir)
}

cor_files <- character(0)
sample_files <- character(0)
metadata_files <- character(0)

for (f in rds_files) {
  obj <- readRDS(f)

  if (is.data.frame(obj) && all(c("comparison", "cor") %in% names(obj))) {
    cor_files <- c(cor_files, f)
  }

  if (is.data.frame(obj) && all(c("comparison", "x", "y") %in% names(obj))) {
    sample_files <- c(sample_files, f)
  }

  if (is.data.frame(obj) && all(c("comparison", "median_cor", "chosen_cor",
                                  "rep", "pred", "species") %in% names(obj)) &&
      !all(c("x", "y") %in% names(obj))) {
    metadata_files <- c(metadata_files, f)
  }
}

if (length(cor_files) == 0) {
  stop("No correlation-summary RDS files found. Expected columns: comparison, cor")
}

if (length(sample_files) == 0) {
  stop("No representative-sample RDS files found. Expected columns: comparison, x, y")
}

cor_by_treatment <- list()
for (f in cor_files) {
  tr <- parse_treatment_from_file(f)
  key <- make_treatment_key(tr$mean_label, tr$sd_label)
  cor_by_treatment[[key]] <- readRDS(f)
}

sample_by_treatment <- list()
sample_file_by_treatment <- list()

## Prefer the new representative sample files over any older downsampled files.
sample_files <- sample_files[
  order(!grepl("_trait_height_posterior_representative_samples\\.rds$",
               sample_files))
]

for (f in sample_files) {
  tr <- parse_treatment_from_file(f)
  key <- make_treatment_key(tr$mean_label, tr$sd_label)

  if (!key %in% names(sample_by_treatment)) {
    sample_by_treatment[[key]] <- readRDS(f)
    sample_file_by_treatment[[key]] <- f
  }
}

metadata_by_treatment <- list()
for (f in metadata_files) {
  tr <- parse_treatment_from_file(f)
  key <- make_treatment_key(tr$mean_label, tr$sd_label)
  metadata_by_treatment[[key]] <- readRDS(f)
}

all_comparisons <- unique(unlist(lapply(cor_by_treatment, function(d) {
  as.character(d$comparison)
})))

all_comparisons <- all_comparisons[!grepl("_ft", all_comparisons, fixed = TRUE)]

preferred <- c(
  "optimum_vs_tolerance",
  "breadth_vs_tolerance",
  "optimum_vs_breadth"
)

if (all(preferred %in% all_comparisons)) {
  comparisons <- preferred
} else {
  comparisons <- unique(c(preferred[preferred %in% all_comparisons], all_comparisons))
}

comparison_labels <- list(
  optimum_vs_height     = expression(theta~"vs"~tau),
  breadth_vs_height     = expression(omega~"vs"~tau),
  optimum_vs_breadth    = expression(theta~"vs"~omega),
  optimum_vs_tolerance  = expression(theta~"vs"~tau),
  breadth_vs_tolerance  = expression(omega~"vs"~tau)
)

xlab_from_comparison <- function(comp) {
  if (grepl("^optimum", comp)) return(expression(theta))
  if (grepl("^breadth", comp)) return(expression(omega))
  sub("_vs_.*$", "", comp)
}

ylab_from_comparison <- function(comp) {
  if (grepl("height$", comp)) return(expression(tau))
  if (grepl("tolerance$", comp)) return(expression(tau))
  if (grepl("breadth$", comp)) return(expression(omega))
  sub("^.*_vs_", "", comp)
}

plot_cor_hist <- function(d, comp, treatment_title) {
  x <- d$cor[d$comparison == comp]
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    plot.new()
    title(main = treatment_title, cex.main = 0.9)
    text(0.5, 0.5, "no finite correlations")
    return(invisible(NULL))
  }

  med_x <- median(x, na.rm = TRUE)

  hist(
    x,
    breaks = 40,
    main = treatment_title,
    xlab = "correlation",
    col = "grey85",
    border = "grey40",
    cex.main = 0.9
  )

  abline(v = 0, lty = 2)
  abline(v = med_x, col = "red", lwd = 1.5)
}

plot_representative_line <- function(d, comp, treatment_title) {
  d <- d[d$comparison == comp, , drop = FALSE]
  good <- is.finite(d$x) & is.finite(d$y)
  d <- d[good, , drop = FALSE]

  if (nrow(d) == 0) {
    plot.new()
    title(main = treatment_title, cex.main = 0.9)
    text(0.5, 0.5, "no representative chain")
    return(invisible(NULL))
  }

  fit <- lm(y ~ x, data = d)

  plot(
    d$x,
    d$y,
    pch = 16,
    cex = 0.25,
    xlab = xlab_from_comparison(comp),
    ylab = ylab_from_comparison(comp),
    main = treatment_title,
    cex.main = 0.9
  )

  abline(fit, col = "red", lwd = 2)

  if (all(c("rep", "pred", "species", "median_cor", "chosen_cor") %in% names(d))) {
    usr <- par("usr")
    text(
      usr[1],
      usr[4],
      labels = paste0(
        "rep = ", d$rep[[1]],
        "\npred = ", d$pred[[1]],
        "\nspecies = ", d$species[[1]],
        "\nmedian r = ", round(d$median_cor[[1]], 3),
        "\nchosen r = ", round(d$chosen_cor[[1]], 3)
      ),
      adj = c(0, 1),
      cex = 0.65
    )
  }
}

## ------------------------------------------------------------
## Label render test
## ------------------------------------------------------------

for (dev_type in c("png", "pdf")) {

  open_plot_device("correlation_label_render_test",
                   dev_type = dev_type,
                   width = 8,
                   height = 3)

  oldpar <- par(no.readonly = TRUE)

  par(mfrow = c(1, 3),
      mar = c(1, 1, 3, 1))

  plot.new()
  title(main = expression(theta~"vs"~tau), cex.main = 2)

  plot.new()
  title(main = expression(omega~"vs"~tau), cex.main = 2)

  plot.new()
  title(main = expression(theta~"vs"~omega), cex.main = 2)

  dev.off()
  par(oldpar)
}

## ------------------------------------------------------------
## Main combined histogram + representative-chain panels
## ------------------------------------------------------------

for (comp in comparisons) {

  safe_comp <- gsub("[^A-Za-z0-9_]+", "_", comp)

  for (dev_type in c("png", "pdf")) {

    open_plot_device(
      filename_base = paste0(safe_comp, "_hist_representative_by_treatment"),
      dev_type = dev_type,
      width = 18,
      height = 10
    )

    oldpar <- par(no.readonly = TRUE)

    par(
      mfrow = c(length(sd_vec), length(height_vec) * 2),
      mar = c(4.0, 4.0, 2.4, 0.8),
      oma = c(0, 0, 2.8, 0)
    )

    for (height_sd in sd_vec) {
      for (height_mean in height_vec) {

        mean_label <- num_to_label(height_mean)
        sd_label   <- num_to_label(height_sd)
        key <- make_treatment_key(mean_label, sd_label)

        treatment_title_hist <- paste0(
          "mean = ", height_mean,
          ", sd = ", height_sd,
          "\nhist"
        )

        treatment_title_line <- paste0(
          "mean = ", height_mean,
          ", sd = ", height_sd,
          "\nmedian representative"
        )

        if (!key %in% names(cor_by_treatment)) {
          plot.new()
          title(main = treatment_title_hist, cex.main = 0.9)
          text(0.5, 0.5, paste0("missing cor RDS\n", key))
        } else {
          plot_cor_hist(cor_by_treatment[[key]], comp, treatment_title_hist)
        }

        if (!key %in% names(sample_by_treatment)) {
          plot.new()
          title(main = treatment_title_line, cex.main = 0.9)
          text(0.5, 0.5, paste0("missing representative RDS\n", key))
        } else {
          plot_representative_line(sample_by_treatment[[key]], comp, treatment_title_line)
        }
      }
    }

    if (!is.null(comparison_labels[[comp]])) {
      mtext(comparison_labels[[comp]],
            outer = TRUE,
            line = 0.7,
            cex = 1.5)
    } else {
      mtext(comp,
            outer = TRUE,
            line = 0.7,
            cex = 1.5)
    }

    dev.off()
    par(oldpar)
  }
}

cat("Wrote combined histogram/representative-chain treatment-panel PNGs and PDFs to: ",
    outdir, "\n", sep = "")

cat("Read RDS files from: ", indir, "\n", sep = "")
cat("Comparisons plotted: ", paste(comparisons, collapse = ", "), "\n", sep = "")

