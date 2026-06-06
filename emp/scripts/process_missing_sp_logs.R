#args <- commandArgs(trailingOnly = TRUE)

#run_id <- as.integer(args[1])


library(MCMCpack, lib = "packages/")
library(robustbase, lib = "packages/")
library(caTools, lib = "packages/")
library(flux, lib = "packages/")
library(evd, lib = "packages/")
library(truncdist, lib = "packages/")
library(MultiRNG, lib = "packages/")
library(Rphylopars, lib = "packages/")
library(BePhyNE, lib = "packages/")

base_dir <- "pletho_miss_out_height_uninform_prior_new"
out_dir  <- "pletho_miss_compiled_logs_height_uninform_prior_new"

burnin   <- 3000
min_rows <- 5000
target_n <- 1000

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

thin_rows <- function(x, target_n = 1000) {
  
  n <- nrow(x)
  
  if (n <= target_n) {
    return(x)
  }
  
  idx <- round(seq(1, n, length.out = target_n))
  
  x[idx, , drop = FALSE]
}


thin_rows <- function(x, target_n = 1000) {
  
  n <- nrow(x)
  
  if (n < target_n) {
    return(NULL)
  }
  
  if (n == target_n) {
    return(x)
  }
  
  idx <- round(seq(1, n, length.out = target_n))
  
  x[idx, , drop = FALSE]
}

extract_cols_for_species <- function(logdf, species_name) {
  
  evo_pattern <- "^pred_[0-9]+_(R[0-9]+|R_cor[0-9]+|R_sd[0-9]+|A[0-9]+)$"
  
  sp_pattern <- paste0(
    "^pred_[0-9]+_dat\\.(opt|brdth|tol|tol_ft)_",
    species_name,
    "$"
  )
  
  keep_cols <- grep(
    paste0(evo_pattern, "|", sp_pattern),
    names(logdf),
    value = TRUE
  )
  
  keep_cols <- keep_cols[order(match(keep_cols, names(logdf)))]
  
  if (length(keep_cols) == 0) {
    return(NULL)
  }
  
  logdf[, keep_cols, drop = FALSE]
}

sp_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

species_logdf_list <- list()

summary_df <- data.frame(
  species = character(),
  file = character(),
  raw_rows = integer(),
  post_burnin_rows = integer(),
  kept = logical(),
  reason = character(),
  stringsAsFactors = FALSE
)

for (sp_dir_idx in 1:length(sp_dirs)) {
  
  sp_dir = sp_dirs[[sp_dir_idx]]
  species_name <- basename(sp_dir)
  message("\nProcessing species: ", species_name)
  
  sp_log_files <- list.files(
    sp_dir,
    pattern = "pars\\.log$",
    full.names = TRUE
  )
  
  if (length(sp_log_files) == 0) {
    summary_df <- rbind(summary_df, data.frame(
      species = species_name,
      file = NA_character_,
      raw_rows = NA_integer_,
      post_burnin_rows = NA_integer_,
      kept = FALSE,
      reason = "no pars.log files found"
    ))
    next
  }
  
  chunks <- list()
  
  for (f_idx in 1:length(sp_log_files)) {
    
    f = sp_log_files[[f_idx]]
    message("  reading: ", basename(f))
    
    logdf <- try(read_BePhyNE_log(file_name = f))
    
    if(class(logdf)=="try-error"){
      
      summary_df <- rbind(summary_df, data.frame(
        species = species_name,
        file = basename(f),
        raw_rows = raw_n,
        post_burnin_rows = NA_integer_,
        kept = FALSE,
        reason = "error reading in file"
      ))
      rm(logdf)
      gc()
      next
      
      
    }
    
    raw_n <- nrow(logdf)
    
    if (raw_n < min_rows) {
      summary_df <- rbind(summary_df, data.frame(
        species = species_name,
        file = basename(f),
        raw_rows = raw_n,
        post_burnin_rows = NA_integer_,
        kept = FALSE,
        reason = "raw rows < 1000"
      ))
      rm(logdf)
      gc()
      next
    }
    
    logdf <- logdf[(burnin + 1):raw_n, , drop = FALSE]
    
    subdf <- extract_cols_for_species(logdf, species_name)
    
    if (is.null(subdf)) {
      summary_df <- rbind(summary_df, data.frame(
        species = species_name,
        file = basename(f),
        raw_rows = raw_n,
        post_burnin_rows = nrow(logdf),
        kept = FALSE,
        reason = "no matching columns found"
      ))
      rm(logdf, subdf)
      gc()
      next
    }
    
    chunks[[length(chunks) + 1]] <- subdf
    
    summary_df <- rbind(summary_df, data.frame(
      species = species_name,
      file = basename(f),
      raw_rows = raw_n,
      post_burnin_rows = nrow(subdf),
      kept = TRUE,
      reason = "kept"
    ))
    
    rm(logdf, subdf)
    gc()
  }
  
  if (length(chunks) == 0) {
    message("  no usable files for ", species_name)
    next
  }
  
  sp_df <- do.call(rbind, chunks)
  
  species_logdf_list[[species_name]] <- sp_df
  
  rm(chunks, sp_df)
  gc()
}


saveRDS(
  species_logdf_list,
  file = file.path(out_dir, paste0("species_logdf_list", ".RDS"))
)

write.csv(
  summary_df,
  file = file.path(out_dir, paste0("compiled_missing_species_summary_", ".csv")),
  row.names = FALSE
)

if (length(species_logdf_list) == 0) {
  stop("No usable species log files were found.")
}

#species_logdf_list = readRDS(paste0(out_dir,"/",' species_logdf_list.RDS'))

#compiled_logdf <- species_logdf_list[[1]]
#
#if (length(species_logdf_list) > 1) {
#  
#  for (i in 2:length(species_logdf_list)) {
#    
#    x <- species_logdf_list[[i]]
#    
#    new_cols <- setdiff(names(x), names(compiled_logdf))
#    
#    compiled_logdf <- cbind(
#      compiled_logdf,
#      x[, new_cols, drop = FALSE]
#    )
#  }
#}
#
#compiled_logdf <- thin_rows(compiled_logdf, target_n = target_n)


compiled_logdf <- thin_rows(species_logdf_list[[1]], target_n = target_n)

if (is.null(compiled_logdf)) {
  stop("First species block has fewer than target_n rows after pooling.")
}

if (length(species_logdf_list) > 1) {
  
  for (i in 2:length(species_logdf_list)) {
    
    x <- thin_rows(species_logdf_list[[i]], target_n = target_n)
    
    if (is.null(x)) {
      stop(
        "Species block has fewer than target_n rows: ",
        names(species_logdf_list)[i]
      )
    }
    
    if (nrow(x) != nrow(compiled_logdf)) {
      stop(
        "Row mismatch before cbind: ",
        names(species_logdf_list)[i],
        " has ", nrow(x),
        " rows; compiled_logdf has ", nrow(compiled_logdf)
      )
    }
    
    new_cols <- setdiff(names(x), names(compiled_logdf))
    
    compiled_logdf <- cbind(
      compiled_logdf,
      x[, new_cols, drop = FALSE]
    )
  }
}

saveRDS(
  compiled_logdf,
  file = file.path(out_dir, paste0("compiled_missing_species_logdf", ".RDS"))
)

