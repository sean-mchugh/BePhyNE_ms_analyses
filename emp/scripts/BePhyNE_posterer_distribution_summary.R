
library(BePhyNE) 

setwd("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp")

sets_full = readRDS("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/data/sets_full.RDS")


pres_data_scaled<-readRDS("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/data/scaled_GBIF_clim_pres.RDS")

pres_data_scaled[,colnames(pres_data_scaled)=="bio1"]

attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]

tree= ENA_Pleth_Tree

scale_atr <- list(scale=c(attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio12"],
                          attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio1"]) ,
                  center= c(attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio12"],
                            attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]))


sp_col    = "species"
occ_col   = "PA"
env_preds = c("bio12", "bio1")
Npred = length(env_preds)
Ntips = length(tree$tip.label)


new_miss_logdf <- readRDS(file = file.path("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_miss_compiled_logs_final", paste0("compiled_missing_species_logdf", ".RDS")))
new_full_logdf <- readRDS("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_full_compiled_logs_final/uninform_height_prior_news_compiled_full_species_logdf.RDS")

GLM_only_ml <- suppressWarnings(BePhyNE:::MLglmStartpars_general(species_data = sets_full$training, 
                                                                 tree = tree))

#cbind(GLM_only_ml$start_pars_bt[[2]], traits_out[[2]][[2]])

scale_atr$scale = scale_atr$scale/10

scale_atr$center = scale_atr$center/10


full_log_summary <- summarize_logdf(new_full_logdf,scale_atr = scale_atr)



make_posterior_table <- function(log_summary){
  
  pred_names <- c(
    "mean annual precipitation",
    "mean annual temperature"
  )
  
  pred_units <- c(
    "mm",
    "°C"
  )
  
  out <- data.frame(
    parameter = character(),
    posterior_median = numeric(),
    HPD_lower = numeric(),
    HPD_upper = numeric(),
    HPD_95 = character(),
    stringsAsFactors = FALSE
  )
  
  for(pred in 1:2){
    
    pred_name <- pred_names[pred]
    pred_unit <- pred_units[pred]
    
    ## A / optimum terms
    A_med <- log_summary$median_parlist$A[[pred]][[1]]
    A_low <- log_summary$HPDlower_parlist$A[[pred]][[1]]
    A_up  <- log_summary$HPDupper_parlist$A[[pred]][[1]]
    
    out <- rbind(
      out,
      data.frame(
        parameter = paste0(pred_name, " A_theta (", pred_unit, ")"),
        posterior_median = A_med[1, 1],
        HPD_lower = A_low[1, 1],
        HPD_upper = A_up[1, 1],
        HPD_95 = paste0(round(A_low[1, 1], 2), "–", round(A_up[1, 1], 2))
      ),
      data.frame(
        parameter = paste0(pred_name, " A_omega (", pred_unit, ")"),
        posterior_median = A_med[1, 2],
        HPD_lower = A_low[1, 2],
        HPD_upper = A_up[1, 2],
        HPD_95 = paste0(round(A_low[1, 2], 2), "–", round(A_up[1, 2], 2))
      )
    )
    
    ## Rsd / phylogenetic SD terms
    Rsd_med <- log_summary$median_parlist$Rsd[[pred]][[1]]
    Rsd_low <- log_summary$HPDlower_parlist$Rsd[[pred]][[1]]
    Rsd_up  <- log_summary$HPDupper_parlist$Rsd[[pred]][[1]]
    
    out <- rbind(
      out,
      data.frame(
        parameter = paste0(pred_name, " sigma_theta"),
        posterior_median = Rsd_med[1],
        HPD_lower = Rsd_low[1],
        HPD_upper = Rsd_up[1],
        HPD_95 = paste0(round(Rsd_low[1], 2), "–", round(Rsd_up[1], 2))
      ),
      data.frame(
        parameter = paste0(pred_name, " sigma_omega"),
        posterior_median = Rsd_med[2],
        HPD_lower = Rsd_low[2],
        HPD_upper = Rsd_up[2],
        HPD_95 = paste0(round(Rsd_low[2], 2), "–", round(Rsd_up[2], 2))
      )
    )
    
    ## Rcor / correlation term
    Rcor_med <- log_summary$median_parlist$Rcor[[pred]][[1]]
    Rcor_low <- log_summary$HPDlower_parlist$Rcor[[pred]][[1]]
    Rcor_up  <- log_summary$HPDupper_parlist$Rcor[[pred]][[1]]
    
    out <- rbind(
      out,
      data.frame(
        parameter = paste0(pred_name, " R_COR"),
        posterior_median = Rcor_med[1, 2],
        HPD_lower = Rcor_low[1, 2],
        HPD_upper = Rcor_up[1, 2],
        HPD_95 = paste0(round(Rcor_low[1, 2], 2), "–", round(Rcor_up[1, 2], 2))
      )
    )
  }
  
  out$posterior_median <- round(out$posterior_median, 2)
  out$HPD_lower <- round(out$HPD_lower, 2)
  out$HPD_upper <- round(out$HPD_upper, 2)
  
  out
}

posterior_table <- make_posterior_table(full_log_summary)

posterior_table
write.csv(posterior_table, "posterior_parameter_table.csv", row.names = FALSE)


