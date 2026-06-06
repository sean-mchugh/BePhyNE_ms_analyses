
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


sp_col    = "species"
occ_col   = "PA"
env_preds = c("bio12", "bio1")
Npred = length(env_preds)
Ntips = length(tree$tip.label)


new_full_logdf = readRDS( "pletho_full_compiled_logs/compiled_full_species_logdf.RDS")
new_miss_logdf = readRDS( "pletho_miss_compiled_logs_height_uninform_prior_new/compiled_missing_species_logdf.RDS")

GLM_only_ml <- suppressWarnings(MLglmStartpars_general(species_data = sets_full$training, 
                                                       tree = tree))

#cbind(GLM_only_ml$start_pars_bt[[2]], traits_out[[2]][[2]])



full_log_summary <- summarize_logdf(new_full_logdf,scale_atr = scale_atr)
miss_log_summary <- summarize_logdf(new_miss_logdf,scale_atr = scale_atr)
GLM_log_summary  <- full_log_summary

GLM_log_summary$median_parlist$traits[[1]][[1]] = GLM_only_ml$start_pars_bt[[1]]
GLM_log_summary$median_parlist$traits[[2]][[1]] = GLM_only_ml$start_pars_bt[[2]]

GLM_log_summary$median_parlist$traits[[1]][[1]][,1] = ((GLM_log_summary$median_parlist$traits[[1]][[1]][,1] * scale_atr$scale[[1]]) + scale_atr$center[[1]])/10 
GLM_log_summary$median_parlist$traits[[2]][[1]][,1] = ((GLM_log_summary$median_parlist$traits[[2]][[1]][,1] * scale_atr$scale[[2]]) + scale_atr$center[[2]])/10

GLM_log_summary$median_parlist$traits[[1]][[1]][,2] = (GLM_log_summary$median_parlist$traits[[1]][[1]][,2] * scale_atr$scale[[1]]) / 10 
GLM_log_summary$median_parlist$traits[[2]][[1]][,2] = (GLM_log_summary$median_parlist$traits[[2]][[1]][,2] * scale_atr$scale[[2]]) / 10



full_log_summary$median_parlist$traits[[1]][[1]][,1:2] = full_log_summary$median_parlist$traits[[1]][[1]][,1:2]/10
full_log_summary$median_parlist$traits[[2]][[1]][,1:2] = full_log_summary$median_parlist$traits[[2]][[1]][,1:2]/10
miss_log_summary$median_parlist$traits[[1]][[1]][,1:2] = miss_log_summary$median_parlist$traits[[1]][[1]][,1:2]/10
miss_log_summary$median_parlist$traits[[2]][[1]][,1:2] = miss_log_summary$median_parlist$traits[[2]][[1]][,1:2]/10


full_log_summary$median_parlist$traits[[1]][[1]] = full_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
full_log_summary$median_parlist$traits[[2]][[1]] = full_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]
miss_log_summary$median_parlist$traits[[1]][[1]] = miss_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
miss_log_summary$median_parlist$traits[[2]][[1]] = miss_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]
GLM_log_summary$median_parlist$traits[[1]][[1]]  = GLM_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,]
GLM_log_summary$median_parlist$traits[[2]][[1]]  = GLM_log_summary$median_parlist$traits[[2]][[1]][tree$tip.label,]



{
  pdf("plots/full_miss_GLM_ridge_pletho_full_new.pdf", height=12, width=10)
  
  plot_summarylist_ridgeplot(tree,
                             list(full_log_summary,miss_log_summary, GLM_log_summary)  ,  
                             model_names       = c("BePhyNE with data   ", "BePhyNE no data", "GLM only"), 
                             predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                             scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                             curve_colors      = c(make.transparent("black", 200/255) , make.transparent("red", 30/255 ), (make.transparent("gray40", 200/255) ) ) ,
                             curve_fill_colors = c(NA, make.transparent("red", 30/255 ), NA  ) ,
                             pt_colors         = c(make.transparent("black", 255/255) , make.transparent("red", 255/255) , (make.transparent("gray40", 255/255) ) ) ,
                             line_types        = c(1,1,2),
                             xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                             predictor_name_cex = 1.0,
                             xlabel_cex         = 0.8
                             
  )
  
  
  dev.off()
}


{
  pdf("plots/full_ridge_pletho_full_new.pdf", height=12, width=10)
  
  plot_summarylist_ridgeplot(tree,
                             list(full_log_summary)  ,  
                             model_names       = c("BePhyNE with data"), 
                             predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                             scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                             curve_colors      = c(make.transparent("black", 200/255))   ,
                             curve_fill_colors = c(NA )   ,
                             pt_colors         = c(make.transparent("black", 255/255))   ,
                             line_types        = c(1),
                             xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                             predictor_name_cex = 1.0,
                             xlabel_cex         = 0.8
                             
  )
  
  
  dev.off()
}




{
  pdf("plots/full_miss_ridge_pletho_full_new.pdf", height=12, width=10)
  
  plot_summarylist_ridgeplot(tree,
                             list(full_log_summary,miss_log_summary)  ,  
                             model_names       = c("BePhyNE with data   ", "BePhyNE no data"), 
                             predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                             scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                             curve_colors      = c(make.transparent("black", 200/255) , make.transparent("red", 30/255 )  ) ,
                             curve_fill_colors = c(NA , make.transparent("red", 30/255 )  ) ,
                             pt_colors         = c(make.transparent("black", 255/255) , make.transparent("red", 255/255)  ) ,
                             line_types        = c(1,1),
                             xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                             predictor_name_cex = 1.0,
                             xlabel_cex         = 0.8
                             
  )
  
  
  dev.off()
}


{
pdf("plots/full_GLM_ridge_pletho_full_new.pdf", height=12, width=10)

plot_summarylist_ridgeplot(tree,
                           list(full_log_summary, GLM_log_summary)  ,  
                           model_names       = c("BePhyNE with data   ", "GLM only"), 
                           predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                           scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                           curve_colors      = c(make.transparent("black", 200/255) , (make.transparent("gray40", 200/255) ) ) ,
                           curve_fill_colors = c(make.transparent("white",  0/255) , NA  ) ,
                           pt_colors         = c(make.transparent("black", 255/255) ,  (make.transparent("gray40", 255/255) ) ) ,
                           line_types        = c(1,2),
                           xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                           predictor_name_cex = 1.0,
                           xlabel_cex         = 0.8
                           
)


dev.off()
}








{
  png("plots/full_miss_GLM_ridge_pletho_full_new.png", height=12, width=10,units = "in", res = 300)
  
  plot_summarylist_ridgeplot(tree,
                             list(full_log_summary,miss_log_summary, GLM_log_summary)  ,  
                             model_names       = c("BePhyNE with data   ", "BePhyNE no data", "GLM only"), 
                             predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                             scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                             curve_colors      = c(make.transparent("black", 200/255) , make.transparent("red", 30/255 ), (make.transparent("gray40", 200/255) ) ) ,
                             curve_fill_colors = c(NA, make.transparent("red", 30/255 ), NA  ) ,
                             pt_colors         = c(make.transparent("black", 255/255) , make.transparent("red", 255/255) , (make.transparent("gray40", 255/255) ) ) ,
                             line_types        = c(1,1,2),
                             xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                             predictor_name_cex = 1.0,
                             xlabel_cex         = 0.8
                             
  )
  
  
  dev.off()
}


{
  png("plots/full_ridge_pletho_full_new.png", height=12, width=10, units = "in", res = 300)
  
  plot_summarylist_ridgeplot(tree,
                             list(full_log_summary)  ,  
                             model_names       = c("BePhyNE with data"), 
                             predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                             scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                             curve_colors      = c(make.transparent("black", 200/255))   ,
                             curve_fill_colors = c(NA )   ,
                             pt_colors         = c(make.transparent("black", 255/255))   ,
                             line_types        = c(1),
                             xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                             predictor_name_cex = 1.0,
                             xlabel_cex         = 0.8
                             
  )
  
  
  dev.off()
}




{
  png("plots/full_miss_ridge_pletho_full_new.png", height=12, width=10,units = "in", res = 300)
  
  plot_summarylist_ridgeplot(tree,
                             list(full_log_summary,miss_log_summary)  ,  
                             model_names       = c("BePhyNE with data   ", "BePhyNE no data"), 
                             predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                             scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                             curve_colors      = c(make.transparent("black", 200/255) , make.transparent("red", 30/255 )  ) ,
                             curve_fill_colors = c(NA , make.transparent("red", 30/255 )  ) ,
                             pt_colors         = c(make.transparent("black", 255/255) , make.transparent("red", 255/255)  ) ,
                             line_types        = c(1,1),
                             xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                             predictor_name_cex = 1.0,
                             xlabel_cex         = 0.8
                             
  )
  
  
  dev.off()
}


{
  png("plots/full_GLM_ridge_pletho_full_new.png", height=12, width=10, units = "in",res = 300)
  
  plot_summarylist_ridgeplot(tree,
                             list(full_log_summary, GLM_log_summary)  ,  
                             model_names       = c("BePhyNE with data   ", "GLM only"), 
                             predictor_names   = c("Average Annual \n Precipitation (mm)\n", "Average Annual \n Temperature (degrees celsius)\n"),   # character vector of predictor names
                             scale_atr         = NA,         # list with $center and $scale (same length as predictors)
                             curve_colors      = c(make.transparent("black", 200/255) , (make.transparent("gray40", 200/255) ) ) ,
                             curve_fill_colors = c(make.transparent("white",  0/255) , NA  ) ,
                             pt_colors         = c(make.transparent("black", 255/255) ,  (make.transparent("gray40", 255/255) ) ) ,
                             line_types        = c(1,2),
                             xlims = list(c(62.7,204.3), c(-1.1,22.4)),
                             predictor_name_cex = 1.0,
                             xlabel_cex         = 0.8,
                             tip_cex=0.65
                             
  )
  
  
  dev.off()
}


