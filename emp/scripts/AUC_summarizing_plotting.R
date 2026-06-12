library(BePhyNE)


setwd("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp")

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





new_full_logdf = readRDS( "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_full_compiled_logs_final/uninform_height_prior_news_compiled_full_species_logdf.RDS")
new_miss_logdf = readRDS( "~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/outfiles/pletho_miss_compiled_logs_final/compiled_missing_species_logdf.RDS")


#new_full_logdf = readRDS( "pletho_full_compiled_logs/uninform_height_prior_news_compiled_full_species_logdf.RDS")
#new_miss_logdf = readRDS( "pletho_miss_compiled_logs/uninform_height_prior_news_compiled_missing_species_logdf.RDS")

GLM_only_ml <- suppressWarnings(BePhyNE:::MLglmStartpars_general(species_data = sets_full$training, 
                                                       tree = tree))

GLM_only_ml$start_pars_bt[[1]][,3][GLM_only_ml$start_pars_bt[[1]][,3]>0.99]=0.99
GLM_only_ml$start_pars_bt[[2]][,3][GLM_only_ml$start_pars_bt[[2]][,3]>0.99]=0.99


full_log_summary <- summarize_logdf(new_full_logdf,scale_atr = NA)
miss_log_summary <- summarize_logdf(new_miss_logdf,scale_atr = NA)


#final_full_log_summary <- summarize_logdf(new_full_logdf_final,scale_atr = NA)
#final_miss_log_summary <- summarize_logdf(new_miss_logdf_final,scale_atr = NA)



#full_log_summary$median_parlist$traits[[1]][[1]][,1:2] = full_log_summary$median_parlist$traits[[1]][[1]][,1:2]/10
#full_log_summary$median_parlist$traits[[2]][[1]][,1:2] = full_log_summary$median_parlist$traits[[2]][[1]][,1:2]/10
#miss_log_summary$median_parlist$traits[[1]][[1]][,1:2] = miss_log_summary$median_parlist$traits[[1]][[1]][,1:2]/10
#miss_log_summary$median_parlist$traits[[2]][[1]][,1:2] = miss_log_summary$median_parlist$traits[[2]][[1]][,1:2]/10

full_log_summary$median_parlist$traits[[1]][[1]] = full_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
full_log_summary$median_parlist$traits[[2]][[1]] = full_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]
miss_log_summary$median_parlist$traits[[1]][[1]] = miss_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
miss_log_summary$median_parlist$traits[[2]][[1]] = miss_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]


full_log_summary$median_parlist$traits[[1]][[1]] - miss_log_summary$median_parlist$traits[[1]][[1]]

#final_miss_log_summary$median_parlist$traits[[1]][[1]] = final_miss_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
#final_miss_log_summary$median_parlist$traits[[2]][[1]] = final_miss_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]
#
#
#cbind(full_log_summary$median_parlist$traits[[1]][[1]][,1], miss_log_summary$median_parlist$traits[[1]][[1]][,1])
#cbind(full_log_summary$median_parlist$traits[[1]][[1]][,2], miss_log_summary$median_parlist$traits[[1]][[1]][,2])
#
#
#cbind(final_full_log_summary$median_parlist$traits[[2]][[1]][,1], full_log_summary$median_parlist$traits[[2]][[1]][,1])
#cbind(final_full_log_summary$median_parlist$traits[[2]][[1]][,2], full_log_summary$median_parlist$traits[[2]][[1]][,2])
#
#cbind(final_full_log_summary$median_parlist$traits[[1]][[1]][,1],
#      full_log_summary$median_parlist$traits[[1]][[1]][,1], 
#      final_miss_log_summary$median_parlist$traits[[1]][[1]][,1], 
#      miss_log_summary$median_parlist$traits[[1]][[1]][,1], 
#      GLM_only_ml$start_pars_bt[[1]][,1])[c(
#        which(grepl("folk", tree$tip.label)),
#        which(grepl("juna", tree$tip.label)),
#        which(grepl("gul", tree$tip.label))
#        
#      ),]
#
#cbind(final_full_log_summary$median_parlist$traits[[1]][[1]][,2],
#      full_log_summary$median_parlist$traits[[1]][[1]][,2], 
#  final_miss_log_summary$median_parlist$traits[[1]][[1]][,2], 
#      miss_log_summary$median_parlist$traits[[1]][[1]][,2], 
#      GLM_only_ml$start_pars_bt[[1]][,2])[c(
#        which(grepl("folk", tree$tip.label)),
#        which(grepl("juna", tree$tip.label)),
#        which(grepl("gul", tree$tip.label))
#        
#      ),]
#
#cbind(final_full_log_summary$median_parlist$traits[[2]][[1]][,1],
#      full_log_summary$median_parlist$traits[[2]][[1]][,1], 
#  final_miss_log_summary$median_parlist$traits[[2]][[1]][,1], 
#            miss_log_summary$median_parlist$traits[[2]][[1]][,1], 
#                         GLM_only_ml$start_pars_bt[[2]][,1])[c(
#                                          which(grepl("folk", tree$tip.label)),
#                                          which(grepl("juna", tree$tip.label)),
#                                          which(grepl("gul", tree$tip.label))
#
#                                          ),]
#
#cbind(final_full_log_summary$median_parlist$traits[[2]][[1]][,2],
#      full_log_summary$median_parlist$traits[[2]][[1]][,2], 
#      final_miss_log_summary$median_parlist$traits[[2]][[1]][,2], 
#            miss_log_summary$median_parlist$traits[[2]][[1]][,2], 
#                         GLM_only_ml$start_pars_bt[[2]][,2])[c(
#        which(grepl("folk", tree$tip.label)),
#        which(grepl("juna", tree$tip.label)),
#        which(grepl("gul", tree$tip.label))
#        
#      ),]
#
#
#
#cbind(final_miss_log_summary$median_parlist$traits[[2]][[1]][,2], miss_log_summary$median_parlist$traits[[2]][[1]][,2], GLM_only_ml$start_pars_bt[[2]][,2])


cbind(
      full_log_summary$median_parlist$traits[[1]][[1]][,1], 
      miss_log_summary$median_parlist$traits[[1]][[1]][,1], 
      GLM_only_ml$start_pars_bt[[1]][,1])[c(
        which(grepl("folk", tree$tip.label)),
        which(grepl("juna", tree$tip.label)),
        which(grepl("gul", tree$tip.label))
        
      ),]

cbind(
      full_log_summary$median_parlist$traits[[1]][[1]][,2], 
      miss_log_summary$median_parlist$traits[[1]][[1]][,2], 
      GLM_only_ml$start_pars_bt[[1]][,2])[c(
        which(grepl("folk", tree$tip.label)),
        which(grepl("juna", tree$tip.label)),
        which(grepl("gul", tree$tip.label))
        
      ),]

cbind(
      full_log_summary$median_parlist$traits[[2]][[1]][,1], 
            miss_log_summary$median_parlist$traits[[2]][[1]][,1], 
                         GLM_only_ml$start_pars_bt[[2]][,1])[c(
                                          which(grepl("folk", tree$tip.label)),
                                          which(grepl("juna", tree$tip.label)),
                                          which(grepl("gul", tree$tip.label))

                                          ),]

cbind(
      full_log_summary$median_parlist$traits[[2]][[1]][,2], 
            miss_log_summary$median_parlist$traits[[2]][[1]][,2], 
                         GLM_only_ml$start_pars_bt[[2]][,2])[c(
        which(grepl("folk", tree$tip.label)),
        which(grepl("juna", tree$tip.label)),
        which(grepl("gul", tree$tip.label))
        
      ),]



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


AUC_glm_list = BePhyNE:::predict_stats(traits = GLM_only_ml$start_pars_bt, pa_data = sets_full$predicting, plot=F)

AUC_glm = unlist(lapply(AUC_glm_list, function(i) i$auc))



AUC_df <- cbind(
  full = AUC_full,
  miss = AUC_miss,
  glm  = AUC_glm
) * 100

rownames(AUC_df) <- tree$tip.label

colMedians(AUC_df)

range(AUC_df[,1])
range(AUC_df[,2])
range(AUC_df[,3])


{

  
  {
    pdf("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/plots/Median_New_AUC_Barplots.pdf",
        height = 11,
        width  = 8.5)
    
    cols <- c("blue", "#F6C344", "black")
    xlim <- c(50, 100)
    
    plotTree.barplot(
      tree,
      AUC_df,
      args.barplot = list(
        beside = TRUE,
        col = cols,
        border = NA,
        xlab = "AUC",
        xlim = xlim,
        mar = c(4.6, 0, 0.5, 2.8),
        space = c(0.25, 1.25)
      ),
      args.plotTree = list(
        fsize = 0.80,
        mar = c(4.6, 0.8, 0.5, 0.35),
        label.offset = 0.55
      )
    )
    
    par(mfg = c(1,1), xpd = NA)
    
    legend(
      "bottomleft",
      inset = c(0.0, -0.06),
      legend = c("BePhyNE With Data",
                 "BePhyNE No Data",
                 "GLM"),
      fill = cols,
      border = NA,
      bty = "n",
      box.lwd=50,
      cex = 1.3,
      pt.cex =10.2,
      y.intersp = 1.1,
      x.intersp = 1.2
    )
    dev.off()
    
    
    
    ## PNG
    png("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/plots/Median_New_AUC_Barplots.png",
        height = 11,
        width = 8.5,
        units = "in",
        res = 300)
    
    plotTree.barplot(
      tree,
      AUC_df,
      args.barplot = list(
        beside = TRUE,
        col = cols,
        border = NA,
        xlab = "AUC",
        xlim = xlim,
        mar = c(4.6, 0, 0.5, 2.8),
        space = c(0.25, 1.25)
      ),
      args.plotTree = list(
        fsize = 0.80,
        mar = c(4.6, 0.8, 0.5, 0.35),
        label.offset = 0.55
      )
    )
    
    par(mfg = c(1, 1), xpd = NA)
    
    legend(
      "bottomleft",
      inset = c(0.0, -0.08),
      legend = c("BePhyNE With Data",
                 "BePhyNE No Data",
                 "GLM"),
      fill = cols,
      border = NA,
      bty = "n",
      box.lwd = 50,
      cex = 1.3,
      pt.cex = 10.2,
      y.intersp = 1.1,
      x.intersp = 1.2
    )
    
    dev.off()
  }  
  
}



########pred 1 only ##########


sets_precip = lapply(sets_full$predicting, function(i) i[1:3])
  
full_traits        = full_log_summary$median_parlist$traits[[1]]
full_predict_stats = predict_stats(full_traits, sets_precip)
AUC_full           = unlist(lapply(full_predict_stats, function(i) i$auc))

miss_traits        = miss_log_summary$median_parlist$traits[[1]]
miss_predict_stats = predict_stats(miss_traits, sets_precip)
AUC_miss           = unlist(lapply(miss_predict_stats, function(i) i$auc))

{
  pdf("plots/precip_final_Median_New_AUC_Barplots.pdf", height=10, width=8)
  
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
  
  abline(v = 90, lty = 2, lwd = 2, col=2)
  
  dev.off()
}


{
  pdf("~/Projects/BePhyNE/BePhyNE_ms_analyses/emp/plots/Median_New_AUC_Barplots.pdf",
      height = 15,
      width  = 8.5)
  
  cols <- c("blue", "#F6C344", "black")
  xlim <- c(50, 100)
  
  plotTree.barplot(
    tree,
    AUC_df,
    args.barplot = list(
      beside = TRUE,
      col = cols,
      border = NA,
      xlab = "AUC",
      xlim = xlim,
      mar = c(4.6, 0, 0.5, 2.8),
      space = c(0.25, 1.25)
    ),
    args.plotTree = list(
      fsize = 0.65,
      mar = c(4.6, 0.8, 0.5, 0.35),
      label.offset = 0.55
    )
  )
  
  usr <- par("usr")
  
  legend(
    x = usr[1] + 1,
    y = usr[3] + 2,
    legend = c("BePhyNE With Data", "BePhyNE No Data", "GLM"),
    fill = cols,
    border = NA,
    bty = "n",
    cex = 0.75,
    xjust = 0,
    yjust = 0,
    xpd = TRUE
  )
  
  dev.off()
}

########pred 2 only ##########


sets_temp = lapply(sets_full$predicting, function(i) i[c(1,2,4)])


full_traits        = full_log_summary$median_parlist$traits[[2]]
full_predict_stats = predict_stats(full_traits, sets_temp)
AUC_full           = unlist(lapply(full_predict_stats, function(i) i$auc))

miss_traits        = miss_log_summary$median_parlist$traits[[2]]
miss_predict_stats = predict_stats(miss_traits, sets_temp)
AUC_miss           = unlist(lapply(miss_predict_stats, function(i) i$auc))

{
  pdf("plots/temp_final_Median_New_AUC_Barplots.pdf", height=10, width=8)
  
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
  
  abline(v = 90, lty = 2, lwd = 2, col=2)
  dev.off()
}

