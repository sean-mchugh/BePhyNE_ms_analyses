library(BePhyNE)


sets_full = readRDS("data/sets_full.RDS")


pres_data_scaled<-readRDS("data/scaled_GBIF_clim_pres.RDS")

pres_data_scaled[,colnames(pres_data_scaled)=="bio1"]

attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]

tree= ENA_Pleth_Tree


scale_atr <- list(scale=c(attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio12"],
                          attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio1"]) ,
                  center= c(attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio12"],
                            attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]))

#
#logdf <- read_BePhyNE_log(
#  file_name = log_files[[i]]
#)
#####summarize only thinned compiled_logs


sp_col    = "species"
occ_col   = "PA"
env_preds = c("bio12", "bio1")
Npred = length(env_preds)
Ntips = length(tree$tip.label)


new_full_logdf = readRDS( "pletho_full_compiled_logs/compiled_full_species_logdf.RDS")
new_miss_logdf = readRDS( "pletho_miss_compiled_logs/compiled_missing_species_logdf.RDS")


full_log_summary <- summarize_logdf(new_full_logdf,scale_atr = NA)
miss_log_summary <- summarize_logdf(new_miss_logdf,scale_atr = NA)


#full_log_summary$median_parlist$traits[[1]][[1]][,1:2] = full_log_summary$median_parlist$traits[[1]][[1]][,1:2]/10
#full_log_summary$median_parlist$traits[[2]][[1]][,1:2] = full_log_summary$median_parlist$traits[[2]][[1]][,1:2]/10
#miss_log_summary$median_parlist$traits[[1]][[1]][,1:2] = miss_log_summary$median_parlist$traits[[1]][[1]][,1:2]/10
#miss_log_summary$median_parlist$traits[[2]][[1]][,1:2] = miss_log_summary$median_parlist$traits[[2]][[1]][,1:2]/10

full_log_summary$median_parlist$traits[[1]][[1]] = full_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
full_log_summary$median_parlist$traits[[2]][[1]] = full_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]
miss_log_summary$median_parlist$traits[[1]][[1]] = miss_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
miss_log_summary$median_parlist$traits[[2]][[1]] = miss_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]



AUC_full_list <- AUC_posterior_median(
  full_log_summary,
  sets_full$predicting
)


AUC_full = unlist(lapply(AUC_full_list, function(i) i$auc))

AUC_miss_list <- AUC_posterior_median(
  miss_log_summary,
  sets_full$predicting
)

AUC_miss = unlist(lapply(AUC_miss_list, function(i) i$auc))




{
  pdf("plots/final_Median_New_AUC_Barplots.pdf", height=10, width=8)
  
  cols <- c("blue", "#F6C344")
  xlim <- c(50, 100)
  fsize <- 0.6
  mar <- c(5.1, 1, 1.1, 0.5)
  label.offset <- 1
  
  
  AUC_df <- cbind(
    full = AUC_full,
    miss = AUC_miss
  ) * 100
  
  
  rownames(AUC_df) <- tree$tip.label
  
  bp <- plotTree.barplot(
    tree,
    AUC_df,
    args.barplot = list(
      beside = TRUE,
      col = cols,
      border = cols,
      xlab = "AUC",
      xlim = xlim,
      mar = c(5.1, 0, 0.1, 4),
      
      # c(within-species gap, between-species gap)
      # 0 = blue/yellow bars touch within species
      # 0.5 = whitespace between species
      space = c(0, 1.5)
    ),
    args.plotTree = list(
      fsize = fsize,
      mar = mar,
      label.offset = label.offset
    )
  )
  
  abline(v = 90, lty = 2, lwd = 4, col=2)
  abline(v = 80, lty = 2, lwd = 4, col=2)
  abline(v = 70, lty = 2, lwd = 4, col=2)
  
  dev.off()
  
}

