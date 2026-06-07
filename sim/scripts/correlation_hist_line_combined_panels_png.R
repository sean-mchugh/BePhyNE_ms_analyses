## Combined correlation histogram + posterior sample line plots by treatment
##
## Reads already-processed RDS files from the current working directory.
## Treatment is parsed from filenames like:
##   0p5_0p15_trait_height_posterior_samples.rds
##   0p5_0p15_trait_height_posterior_samples_downsampled.rds
##
## The script auto-detects which files contain correlation summaries
## by checking for columns: comparison, cor
## and which files contain sampled x/y values by checking for: comparison, x, y.
##
## Output:
##   plots/<comparison>_hist_line_by_treatment.png
##
## Layout per PDF:
##   rows    = sd values
##   columns = mean values, with two columns per mean:
##             histogram | line/scatter plot
##
## Predictors are pooled by NOT filtering on pred.
## _ft comparisons are ignored.

height_vec <- c("0.5", "0.7", "0.95")
sd_vec     <- c("0.15", "0.5", "1")

outdir <- "~/Downloads/Height_plots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

set.seed(1)

label_to_num <- function(x) as.numeric(gsub("p", ".", x, fixed = TRUE))
num_to_label <- function(x) gsub(".", "p", as.character(x), fixed = TRUE)

parse_treatment_from_file <- function(file) {
  b <- basename(file)
  b <- sub("_trait_height_posterior_samples_downsampled\\.rds$", "", b)
  b <- sub("_trait_height_posterior_samples\\.rds$", "", b)
  b <- sub("_trait_height_posterior_correlations\\.rds$", "", b)
  b <- sub("_trait_height_posterior_correlation_summary\\.rds$", "", b)

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

rds_files <- list.files(
  path = ".",
  pattern = "_trait_height_posterior_.*\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No trait-height posterior RDS files found in current working directory: ", getwd())
}

cor_files <- character(0)
sample_files <- character(0)

for (f in rds_files) {
  obj <- readRDS(f)

  if (is.data.frame(obj) && all(c("comparison", "cor") %in% names(obj))) {
    cor_files <- c(cor_files, f)
  }

  if (is.data.frame(obj) && all(c("comparison", "x", "y") %in% names(obj))) {
    sample_files <- c(sample_files, f)
  }
}

if (length(cor_files) == 0) {
  stop("No correlation-summary RDS files found. Expected data.frames with columns: comparison, cor")
}

if (length(sample_files) == 0) {
  stop("No posterior-sample RDS files found. Expected data.frames with columns: comparison, x, y")
}

make_treatment_key <- function(mean_label, sd_label) {
  paste(mean_label, sd_label, sep = "__")
}

cor_by_treatment <- list()
for (f in cor_files) {
  tr <- parse_treatment_from_file(f)
  key <- make_treatment_key(tr$mean_label, tr$sd_label)
  cor_by_treatment[[key]] <- readRDS(f)
}

sample_by_treatment <- list()
for (f in sample_files) {
  tr <- parse_treatment_from_file(f)
  key <- make_treatment_key(tr$mean_label, tr$sd_label)
  sample_by_treatment[[key]] <- readRDS(f)
}

all_comparisons <- unique(unlist(lapply(cor_by_treatment, function(d) as.character(d$comparison))))
all_comparisons <- all_comparisons[!grepl("_ft", all_comparisons, fixed = TRUE)]

## Prefer the height names if present; otherwise use whatever non-_ft names exist.
preferred <- c("optimum_vs_height", "breadth_vs_height")
if (all(preferred %in% all_comparisons)) {
  comparisons <- preferred
} else {
  comparisons <- all_comparisons
}

comparison_labels <- list(
  optimum_vs_height     = expression(italic(theta)~"vs"~italic(h)),
  breadth_vs_height     = expression(italic(omega)~"vs"~italic(h)),
  optimum_vs_tolerance  = expression(italic(theta)~"vs"~italic(h)),
  breadth_vs_tolerance  = expression(italic(omega)~"vs"~italic(h))
)

xlab_from_comparison <- function(comp) {
  if (grepl("^optimum", comp)) return("optimum")
  if (grepl("^breadth", comp)) return("breadth")
  sub("_vs_.*$", "", comp)
}

ylab_from_comparison <- function(comp) {
  if (grepl("height$", comp)) return("height")
  if (grepl("tolerance$", comp)) return("tolerance")
  sub("^.*_vs_", "", comp)
}

plot_cor_hist <- function(d, comp, treatment_title) {
  x <- d$cor[d$comparison == comp]
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    plot.new()
    title(main = treatment_title, cex.main = 0.9)
    text(0.5, 0.5, "no finite correlations")
    mtext(comp, side = 3, line = -1.2, cex = 0.75)
    return(invisible(NULL))
  }

  hist(
    x,
    breaks = 40,
    main = treatment_title,
    xlab = "correlation",
    cex.main = 0.9
  )
  abline(v = 0, lty = 2)
}

plot_sample_line <- function(d, comp, treatment_title, point_fraction = 0.10) {
  d <- d[d$comparison == comp, , drop = FALSE]
  good <- is.finite(d$x) & is.finite(d$y)
  d <- d[good, , drop = FALSE]

  if (nrow(d) == 0) {
    plot.new()
    title(main = treatment_title, cex.main = 0.9)
    text(0.5, 0.5, "no finite samples")
    mtext(comp, side = 3, line = -1.2, cex = 0.75)
    return(invisible(NULL))
  }

  keep <- sample(seq_len(nrow(d)), size = max(1, ceiling(point_fraction * nrow(d))))
  fit <- lm(y ~ x, data = d)

  plot(
    d$x[keep],
    d$y[keep],
    pch = 16,
    cex = 0.25,
    xlab = xlab_from_comparison(comp),
    ylab = ylab_from_comparison(comp),
    main = treatment_title,
    cex.main = 0.9
  )

  abline(fit, col = "red", lwd = 2)
}

for (comp in comparisons) {

  safe_comp <- gsub("[^A-Za-z0-9_]+", "_", comp)

  png(
    filename = file.path(outdir, paste0(safe_comp, "_hist_line_by_treatment.png")),
    width = 18,
    height = 10,
    units = "in",
    res = 300
  )

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

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

      treatment_title_hist <- paste0("mean = ", height_mean, ", sd = ", height_sd, "\nhist")
      treatment_title_line <- paste0("mean = ", height_mean, ", sd = ", height_sd, "\nline")

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
        text(0.5, 0.5, paste0("missing sample RDS\n", key))
      } else {
        plot_sample_line(sample_by_treatment[[key]], comp, treatment_title_line)
      }
    }
  }

  if (!is.null(comparison_labels[[comp]])) {
    mtext(comparison_labels[[comp]], outer = TRUE, line = 0.7, cex = 1.5)
  } else {
    mtext(comp, outer = TRUE, line = 0.7, cex = 1.5)
  }

  dev.off()
}

cat("Wrote combined histogram/line treatment-panel PNGs to: ", outdir, "\n", sep = "")
cat("Comparisons plotted: ", paste(comparisons, collapse = ", "), "\n", sep = "")
