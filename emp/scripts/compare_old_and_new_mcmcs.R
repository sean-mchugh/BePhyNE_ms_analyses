


## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

run_id  <- as.integer(args[1])
#miss    <- as.integer(args[2])

set.seed(run_id)

library(MCMCpack, lib="packages/")
library(robustbase, lib="packages/")
library(caTools, lib="packages/")

library(flux, lib="packages/")
library(evd, lib="packages/")
library(truncdist, lib="packages/")
library(MultiRNG, lib="packages/")
library(Rphylopars, lib="packages/")
library(BePhyNE, lib="packages/")


## -----------------------------------------------------------------------------
## Load and format data
## -----------------------------------------------------------------------------

tree <- ENA_Pleth_Tree
#pa_data <- read.csv("data/Pleth_data_vignette.csv")


miss=1


sp_col <- "species"
occ_col <- "PA"
env_preds <- c("bio12", "bio1")

Npred <- length(env_preds)
Ntips <- length(tree$tip.label)



sp_dirs = list.dirs("pletho_full_out",recursive = F)
#miss=82
#sp_dir = sp_dirs[[miss]]

sp_files = list.files("pletho_full_out", full.names = T)

sp_log_files = sp_files[grepl("pars.log",sp_files)]

logdf_list=list()

for(i in 1:length(sp_log_files)){
  
  logdf <- read_BePhyNE_log(
    file_name = paste0(sp_log_files[[i]])
  )
  print(nrow(logdf))
  logdf_list[[i]] = logdf
  
}


logdf_list_wburnin =  logdf_list

logdf_list = logdf_list_wburnin 


logdf_list= lapply(logdf_list, function(i) i[-(1:2500),])

logdf_list = logdf_list[unlist(lapply( logdf_list, nrow))>100]

full_logdf = do.call(rbind, logdf_list)


sp_startPar_files = sp_files[grepl("start",sp_files)]

startPars_list=list()

for(i in 1:length(sp_log_files)){
  
  #print(nrow(logdf))
  startPars_list[[i]] = readRDS(sp_startPar_files[[i]])
  
}

startPars_list[[1]]$sim_dat$sim_dat_bt[[1]][,3]

lapply( logdf_list, nrow)


logdf_list[[1]]$pred_1_dat.tol_Plethodon_savannah
logdf_list[[5]]$pred_1_dat.tol_Desmognathus_fuscus

logdf_list[[5]]$pred_1_dat.tol_Desmognathus_imitator
logdf_list[[5]]$pred_1_dat.tol_Desmognathus_imitator
logdf_list[[5]]$pred_1_dat.tol_Desmognathus_ocoee

logdf_list = logdf_list[-5]



#full_logdf =full_logdf[sample(1:nrow(full_logdf), size = 650,replace = F),]

old_missing_data_mcmc = readRDS("~/Manuscripts/BePhyNE/R_code/presub_HPC_analyses_pub_FullandMiss/pletho_analyses/Josef_posterior_ecdf_miss_dat/missing_data_mcmc.RDS")
old_data_mcmc         = readRDS("~/Manuscripts/BePhyNE/R_code/presub_HPC_analyses_pub_FullandMiss/pletho_analyses/Josef_posterior_ecdf_miss_dat/full_data_mcmc.RDS")


old_missing_data_mcmc = readRDS("~/Downloads/rstudio-export (33)/missing_data_mcmc.RDS")
old_data_mcmc         = readRDS("~/Downloads/rstudio-export (33)/full_data_mcmc.RDS")


traits_out = readRDS("~/Manuscripts/BePhyNE/R_code/traits_out.RDS")



old_data_mcmc = lapply(old_missing_data_mcmc, function(i) i[[1]])

plot(density(old_missing_data_mcmc[[1]][[1]][,1]))
plot(density(old_data_mcmc[[1]][,1]))

lines(density(full_logdf$pred_1_dat.opt_Stereochilus_marginatus))

x=1

{
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred1_opt.pdf")
    par(mfrow=c(3,3))
    
    for(i in 1:82){
      
      min= min( c( min(full_logdf[,paste0("pred_1_dat.opt_", tree$tip.label[[i]])]) , min(old_data_mcmc[[1]][,i]) ) )
      max= max( c( max(full_logdf[,paste0("pred_1_dat.opt_", tree$tip.label[[i]])]) , max(old_data_mcmc[[1]][,i])))
      
      plot(c(0,nrow(full_logdf)), c(min,max), type="n", main= tree$tip.label[[i]])
      lines(1:nrow(full_logdf), full_logdf[,paste0("pred_1_dat.opt_", tree$tip.label[[i]])], col=3)
      lines(1:length(old_data_mcmc[[1]][,i]), old_data_mcmc[[1]][,i])
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred1_brdth.pdf")
    par(mfrow=c(3,3))
    
    for(i in 1:82){
      
      brdth_new = exp(full_logdf[,paste0("pred_1_dat.brdth_", tree$tip.label[[i]])])
      
      min= min( c( min(brdth_new) , min(old_data_mcmc[[3]][,i]) ) )
      max= max( c( max(brdth_new), max(old_data_mcmc[[3]][,i])))
      
      plot(c(0,nrow(full_logdf)), c(min,max), type="n", main= tree$tip.label[[i]])
      lines(1:nrow(full_logdf), brdth_new, col=3)
      lines(1:length(old_data_mcmc[[3]][,i]), old_data_mcmc[[3]][,i])
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred1_tol.pdf")
    par(mfrow=c(3,3))
    
    for(i in 1:82){
      
      tol_ft = full_logdf[,paste0("pred_1_dat.tol_", tree$tip.label[[i]])]
      tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
      min= min( c( min(tol_new) , min(old_data_mcmc[[5]][,i]) ) )
      max= max(c(max(tol_new), max(old_data_mcmc[[5]][,i])))
      
      plot(c(0,nrow(full_logdf)), c(min,max), type="n", main = tree$tip.label[[i]])
      lines(1:nrow(full_logdf), tol_new, col=3)
      lines(1:length(old_data_mcmc[[5]][,i]), old_data_mcmc[[5]][,i])
      
    }
    dev.off()
  }
  
  
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred2_opt.pdf")
    par(mfrow=c(3,3))
    
    for(i in 1:82){
      
      min= min( c( min(full_logdf[,paste0("pred_2_dat.opt_", tree$tip.label[[i]])]) , min(old_data_mcmc[[2]][,i]) ) )
      max= max(c(max(full_logdf[,paste0("pred_2_dat.opt_", tree$tip.label[[i]])]), max(old_data_mcmc[[2]][,i])))
      
      plot(c(0,nrow(full_logdf)), c(min,max), type="n", main= tree$tip.label[[i]])
      lines(1:nrow(full_logdf), full_logdf[,paste0("pred_2_dat.opt_", tree$tip.label[[i]])], col=3)
      lines(1:length(old_data_mcmc[[2]][,i]), old_data_mcmc[[2]][,i])
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred2_brdth.pdf")
    par(mfrow=c(3,3))
    
    for(i in 1:82){
      
      brdth_new = exp(full_logdf[,paste0("pred_2_dat.brdth_", tree$tip.label[[i]])])
      
      min= min( c( min( brdth_new) , min(old_data_mcmc[[4]][,i]) ) )
      max= max(c(max( brdth_new), max(old_data_mcmc[[4]][,i])))
      
      plot(c(0,nrow(full_logdf)), c(min,max), type="n", main= tree$tip.label[[i]])
      lines(1:nrow(full_logdf),  brdth_new, col=3)
      lines(1:length(old_data_mcmc[[4]][,i]), old_data_mcmc[[4]][,i])
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred2_tol.pdf")
    par(mfrow=c(3,3))
    
    for(i in 1:82){
      
      tol_ft = full_logdf[,paste0("pred_2_dat.tol_", tree$tip.label[[i]])]
      tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
      min= min( c( min(tol_new) , min(old_data_mcmc[[6]][,i]) ) )
      max= max(c(max(tol_new), max(old_data_mcmc[[6]][,i])))
      
      plot(c(0,nrow(full_logdf)), c(min,max), type="n", main= tree$tip.label[[i]])
      lines(1:nrow(full_logdf), tol_new, col=3)
      lines(1:length(old_data_mcmc[[6]][,i]), old_data_mcmc[[6]][,i])
      
    }
    dev.off()
  }
  
}


{
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred1_opt_chainoverlay.pdf")
    par(mfrow=c(4,2))
    
    for(i in 1:82){
      
      min= min( c( min(full_logdf[,paste0("pred_1_dat.opt_", tree$tip.label[[i]])]) , min(old_data_mcmc[[1]][,i]) ) )
      max= max( c( max(full_logdf[,paste0("pred_1_dat.opt_", tree$tip.label[[i]])]) , max(old_data_mcmc[[1]][,i])))
      
      plot(c(0,nrow(logdf_list[[x]])), c(min,max), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        lines(1:nrow(logdf_list[[x]]), logdf_list[[x]][,paste0("pred_1_dat.opt_", tree$tip.label[[i]])], col=3)
      }
      lines(1:length(old_data_mcmc[[1]][,i]), old_data_mcmc[[1]][,i])
      
      
      dens = density( old_data_mcmc[[1]][,i])
      
      plot( c(min,max),c(0,max(dens$y)*1.5), type="n", main= tree$tip.label[[i]])
      
      #plot( c(min,max),c(0,5), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        lines(density(logdf_list[[x]][,paste0("pred_1_dat.opt_", tree$tip.label[[i]])]), col=3)
      }
      lines(density( old_data_mcmc[[1]][,i]))
      
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred1_brdth_chainoverlay.pdf")
    par(mfrow=c(4,2))
    
    for(i in 1:82){
      
      brdth_new = exp(full_logdf[,paste0("pred_1_dat.brdth_", tree$tip.label[[i]])])
      
      min= min( c( min(brdth_new) , min(old_data_mcmc[[3]][,i]) ) )
      max= max( c( max(brdth_new), max(old_data_mcmc[[3]][,i])))
      
      plot(c(0,nrow(logdf_list[[x]])), c(min,max), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        
        brdth_new = exp(logdf_list[[x]][,paste0("pred_1_dat.brdth_", tree$tip.label[[i]])])
        
        lines(1:nrow(logdf_list[[x]]),    brdth_new, col=3)
      }
      
      lines(1:length(old_data_mcmc[[3]][,i]), old_data_mcmc[[3]][,i])
      
      
      dens = density( old_data_mcmc[[3]][,i])
      
      plot( c(min,max),c(0,max(dens$y)*1.5), type="n", main= tree$tip.label[[i]])
      
      #plot( c(min,max),c(0,5), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        brdth_new = exp(logdf_list[[x]][,paste0("pred_1_dat.brdth_", tree$tip.label[[i]])])
        
        lines(density(brdth_new), col=3)
      }
      lines(density( old_data_mcmc[[3]][,i]))
      
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred1_tol_chainoverlay.pdf")
    par(mfrow=c(4,2))
    
    for(i in 1:82){
      
      tol_ft = full_logdf[,paste0("pred_1_dat.tol_", tree$tip.label[[i]])]
      tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
      min= min( c( min(tol_new) , min(old_data_mcmc[[5]][,i]) ) )
      max= max(c(max(tol_new), max(old_data_mcmc[[5]][,i])))
      
      plot(c(0,nrow(logdf_list[[x]])), c(min,max), type="n", main = tree$tip.label[[i]])
      
      for(x in 1:length(logdf_list)){
        tol_ft = logdf_list[[x]][,paste0("pred_1_dat.tol_", tree$tip.label[[i]])]
        tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
        
        lines(1:nrow(logdf_list[[x]]),  tol_new, col=3)
      }
      
      #lines(1:nrow(full_logdf), tol_new, col=3)
      lines(1:length(old_data_mcmc[[5]][,i]), old_data_mcmc[[5]][,i])
      
      
      dens = density( old_data_mcmc[[5]][,i])
      
      plot( c(min,max),c(0,max(dens$y)*1.5), type="n", main= tree$tip.label[[i]])
      
      #plot( c(min,max),c(0,5), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        tol_ft = logdf_list[[x]][,paste0("pred_1_dat.tol_", tree$tip.label[[i]])]
        tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
        
        lines(density(tol_new ), col=3)
      }
      lines(density( old_data_mcmc[[5]][,i]))
      
    }
    dev.off()
  }
  
  
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred2_opt_chainoverlay.pdf")
    par(mfrow=c(4,2))
    
    for(i in 1:82){
      
      min= min( c( min(full_logdf[,paste0("pred_2_dat.opt_", tree$tip.label[[i]])]) , min(old_data_mcmc[[2]][,i]) ) )
      max= max(c(max(full_logdf[,paste0("pred_2_dat.opt_", tree$tip.label[[i]])]), max(old_data_mcmc[[2]][,i])))
      
      plot(c(0,nrow(logdf_list[[x]])), c(min,max), type="n", main= tree$tip.label[[i]])
      
      for(x in 1:length(logdf_list)){
        lines(1:nrow(logdf_list[[x]]), logdf_list[[x]][,paste0("pred_2_dat.opt_", tree$tip.label[[i]])], col=3)
      }
      
      #lines(1:nrow(full_logdf), full_logdf[,paste0("pred_2_dat.opt_", tree$tip.label[[i]])], col=3)
      lines(1:length(old_data_mcmc[[2]][,i]), old_data_mcmc[[2]][,i])
      
      
      dens = density( old_data_mcmc[[2]][,i])
      
      plot( c(min,max),c(0,max(dens$y)*1.5), type="n", main= tree$tip.label[[i]])
      
      #plot( c(min,max),c(0,5), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        lines(density(logdf_list[[x]][,paste0("pred_2_dat.opt_", tree$tip.label[[i]])]), col=3)
      }
      lines(density( old_data_mcmc[[2]][,i]))
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred2_brdth_chainoverlay.pdf")
    par(mfrow=c(4,2))
    
    for(i in 1:82){
      
      brdth_new = exp(full_logdf[,paste0("pred_2_dat.brdth_", tree$tip.label[[i]])])
      
      min= min( c( min( brdth_new) , min(old_data_mcmc[[4]][,i]) ) )
      max= max(c(max( brdth_new), max(old_data_mcmc[[4]][,i])))
      
      plot(c(0,nrow(logdf_list[[x]])), c(min,max), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        
        brdth_new = exp(logdf_list[[x]][,paste0("pred_2_dat.brdth_", tree$tip.label[[i]])])
        
        lines(1:nrow(logdf_list[[x]]),    brdth_new, col=3)
      }
      lines(1:length(old_data_mcmc[[4]][,i]), old_data_mcmc[[4]][,i])
      
      
      dens = density( old_data_mcmc[[4]][,i])
      
      plot( c(min,max),c(0,max(dens$y)*1.5), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        brdth_new = exp(logdf_list[[x]][,paste0("pred_2_dat.brdth_", tree$tip.label[[i]])])
        
        lines(density(brdth_new), col=3)
      }
      lines(density( old_data_mcmc[[4]][,i]))
      
    }
    dev.off()
  }
  
  {
    
    pdf("~/Downloads/new_Trace_comp_pletho_pred2_tol_chainoverlay.pdf")
    par(mfrow=c(4,2))
    
    for(i in 1:82){
      
      tol_ft = full_logdf[,paste0("pred_2_dat.tol_", tree$tip.label[[i]])]
      tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
      min= min( c( min(tol_new) , min(old_data_mcmc[[6]][,i]) ) )
      max= max(c(max(tol_new), max(old_data_mcmc[[6]][,i])))
      
      plot(c(0,nrow(logdf_list[[x]])), c(min,max), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        tol_ft = logdf_list[[x]][,paste0("pred_2_dat.tol_", tree$tip.label[[i]])]
        tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
        
        lines(1:nrow(logdf_list[[x]]),  tol_new, col=3)
      }
      lines(1:length(old_data_mcmc[[6]][,i]), old_data_mcmc[[6]][,i])
      
      
      dens = density( old_data_mcmc[[6]][,i])
      max(dens$y)*1.5
      plot( c(min,max),c(0,max(dens$y)*1.5), type="n", main= tree$tip.label[[i]])
      for(x in 1:length(logdf_list)){
        tol_ft = logdf_list[[x]][,paste0("pred_2_dat.tol_", tree$tip.label[[i]])]
        tol_new = 0.05 + 0.95 * (exp(-1 * tol_ft)/(1 + exp(-1 * tol_ft)))
        
        lines(density(tol_new ), col=3)
      }
      lines(density( old_data_mcmc[[6]][,i]))
      
    }
    dev.off()
  }
  
}



plot(1:length(old_missing_data_mcmc[[1]][[82]][,82]), old_missing_data_mcmc[[1]][[82]][,82])
lines(1:length(old_data_mcmc[[1]][,82]), old_data_mcmc[[1]][,82])


log_summary <- summarize_logdf(logdf,scale_atr = scale_atr)

log_summary$median_parlist$Rsd

cbind(colMedians(old_data_mcmc[[1]]), log_summary$median_parlist$traits[[1]][[1]][,1])


log_summary$median_parlist$traits[[1]][[1]][,1]
log_summary$median_parlist$Rsd


#we take the scaling attributes from ALL presence data prior to downsampling 

pres_data_scaled<-readRDS("data/scaled_GBIF_clim_pres.RDS")

pres_data_scaled[,colnames(pres_data_scaled)=="bio1"]

attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]


scale_atr <- list(scale=c(attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio12"],
                          attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio1"]) ,
                  center= c(attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio12"],
                            attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]))


#(-3 * scale_atr$scale[[2]] + scale_atr$center[[2]]) /10

#(1.2 * scale_atr$scale[[1]] ) /10


{
  
  trait_info <- data.frame(
    pred     = c(1, 2, 1, 2, 1, 2),
    trait    = c("opt", "opt", "brdth", "brdth", "tol", "tol"),
    old_idx  = c(1, 2, 3, 4, 5, 6),
    main_lab = c("C 1", "C 2", "W 1", "W 2", "H 1", "H 2"),
    stringsAsFactors = FALSE
  )
  
  {
    
    pdf(
      file = "~/Downloads/miss_1_old_new_Trace_comp_pletho_all_species_traces_1.pdf",
      width = 10,
      height = 9
    )
    
    for (i in 1:82) {
      
      sp <- tree$tip.label[[i]]
      
      par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
      
      for (j in seq_len(nrow(trait_info))) {
        
        pred    <- trait_info$pred[j]
        trait   <- trait_info$trait[j]
        old_idx <- trait_info$old_idx[j]
        
        mainlab <- paste(trait_info$main_lab[j], sp)
        
        col_name <- paste0("pred_", pred, "_dat.", trait, "_", sp)
        
        old_vals <- old_data_mcmc[[old_idx]][, i]
        full_vals <- full_logdf[, col_name]
        
        if (trait == "brdth") {
          full_vals <- exp(full_vals)
        }
        
        if (trait == "tol") {
          full_vals <- 0.05 + 0.95 * (exp(-1 * full_vals) / (1 + exp(-1 * full_vals)))
        }
        
        ymin <- min(c(full_vals, old_vals), na.rm = TRUE)
        ymax <- max(c(full_vals, old_vals), na.rm = TRUE)
        
        xmax <- max(sapply(logdf_list, nrow), length(old_vals))
        
        plot(
          c(0, xmax), c(ymin, ymax),
          type = "n",
          main = mainlab,
          xlab = "Index",
          ylab = ""
        )
        
        for (x in seq_along(logdf_list)) {
          
          new_vals <- logdf_list[[x]][, col_name]
          
          if (trait == "brdth") {
            new_vals <- exp(new_vals)
          }
          
          if (trait == "tol") {
            new_vals <- 0.05 + 0.95 * (exp(-1 * new_vals) / (1 + exp(-1 * new_vals)))
          }
          
          lines(seq_along(new_vals), new_vals, col = 3)
        }
        
        lines(seq_along(old_vals), old_vals)
      }
    }
    
    dev.off()
  }
  
  
  
  pdf(
    file = "~/Downloads/miss_1_old_new_Trace_comp_pletho_all_species_densities_1.pdf",
    width = 10,
    height = 9
  )
  
  for (i in 1:82) {
    
    sp <- tree$tip.label[[i]]
    
    par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
    
    for (j in seq_len(nrow(trait_info))) {
      
      pred    <- trait_info$pred[j]
      trait   <- trait_info$trait[j]
      old_idx <- trait_info$old_idx[j]
      
      mainlab <- paste(trait_info$main_lab[j], sp)
      
      col_name <- paste0("pred_", pred, "_dat.", trait, "_", sp)
      
      old_vals <- old_data_mcmc[[old_idx]][, i]
      full_vals <- full_logdf[, col_name]
      
      if (trait == "brdth") {
        full_vals <- exp(full_vals)
      }
      
      if (trait == "tol") {
        full_vals <- 0.05 + 0.95 * (exp(-1 * full_vals) / (1 + exp(-1 * full_vals)))
      }
      
      xmin <- min(c(full_vals, old_vals), na.rm = TRUE)
      xmax <- max(c(full_vals, old_vals), na.rm = TRUE)
      
      dens_old <- density(old_vals)
      ymax <- max(dens_old$y, na.rm = TRUE)
      
      dens_new_list <- vector("list", length(logdf_list))
      
      for (x in seq_along(logdf_list)) {
        
        new_vals <- logdf_list[[x]][, col_name]
        
        if (trait == "brdth") {
          new_vals <- exp(new_vals)
        }
        
        if (trait == "tol") {
          new_vals <- 0.05 + 0.95 * (exp(-1 * new_vals) / (1 + exp(-1 * new_vals)))
        }
        
        dens_new_list[[x]] <- density(new_vals)
        ymax <- max(ymax, dens_new_list[[x]]$y, na.rm = TRUE)
      }
      
      plot(
        c(xmin, xmax), c(0, ymax * 1.5),
        type = "n",
        main = mainlab,
        xlab = "",
        ylab = "Density"
      )
      
      for (x in seq_along(dens_new_list)) {
        lines(dens_new_list[[x]], col = 3)
      }
      
      lines(dens_old)
    }
  }
  
  dev.off()
  
}





{
  
  trait_info <- data.frame(
    pred     = c(1, 2, 1, 2, 1, 2),
    trait    = c("opt", "opt", "brdth", "brdth", "tol", "tol"),
    old_idx  = c(1, 2, 3, 4, 5, 6),
    main_lab = c("C 1", "C 2", "W 1", "W 2", "H 1", "H 2"),
    stringsAsFactors = FALSE
  )
  
  {
    
    pdf(
      file = "~/Downloads/newprior__Trace_comp_pletho_all_species_traces_single_post.pdf",
      width = 10,
      height = 9
    )
    
    for (i in 1:82) {
      
      sp <- tree$tip.label[[i]]
      
      par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
      
      for (j in seq_len(nrow(trait_info))) {
        
        pred    <- trait_info$pred[j]
        trait   <- trait_info$trait[j]
        old_idx <- trait_info$old_idx[j]
        
        mainlab <- paste(trait_info$main_lab[j], sp)
        
        col_name <- paste0("pred_", pred, "_dat.", trait, "_", sp)
        
        old_vals <- old_data_mcmc[[old_idx]][, i]
        full_vals <- full_logdf[, col_name]
        
        if (trait == "brdth") {
          full_vals <- exp(full_vals)
        }
        
        if (trait == "tol") {
          full_vals <- 0.05 + 0.95 * (exp(-1 * full_vals) / (1 + exp(-1 * full_vals)))
        }
        
        ymin <- min(c(full_vals, old_vals), na.rm = TRUE)
        ymax <- max(c(full_vals, old_vals), na.rm = TRUE)
        
        xmax <- max(sapply(logdf_list, nrow), length(old_vals))
        
        plot(
          c(0, xmax), c(ymin, ymax),
          type = "n",
          main = mainlab,
          xlab = "Index",
          ylab = ""
        )
        
        for (x in 1) {
          
          new_vals <- full_logdf[, col_name]
          
          if (trait == "brdth") {
            new_vals <- exp(new_vals)
          }
          
          if (trait == "tol") {
            new_vals <- 0.05 + 0.95 * (exp(-1 * new_vals) / (1 + exp(-1 * new_vals)))
          }
          
          lines(seq_along(new_vals), new_vals, col = 3)
        }
        
        lines(seq_along(old_vals), old_vals)
      }
    }
    
    dev.off()
  }
  
  
  
  pdf(
    file = "~/Downloads/newprior_Trace_comp_pletho_all_species_densities_single_post.pdf",
    width = 10,
    height = 9
  )
  
  for (i in 1:82) {
    
    sp <- tree$tip.label[[i]]
    
    par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
    
    for (j in seq_len(nrow(trait_info))) {
      
      pred    <- trait_info$pred[j]
      trait   <- trait_info$trait[j]
      old_idx <- trait_info$old_idx[j]
      
      mainlab <- paste(trait_info$main_lab[j], sp)
      
      col_name <- paste0("pred_", pred, "_dat.", trait, "_", sp)
      
      old_vals <- old_data_mcmc[[old_idx]][, i]
      full_vals <- full_logdf[, col_name]
      
      if (trait == "brdth") {
        full_vals <- exp(full_vals)
      }
      
      if (trait == "tol") {
        full_vals <- 0.05 + 0.95 * (exp(-1 * full_vals) / (1 + exp(-1 * full_vals)))
      }
      
      xmin <- min(c(full_vals, old_vals), na.rm = TRUE)
      xmax <- max(c(full_vals, old_vals), na.rm = TRUE)
      
      dens_old <- density(old_vals)
      ymax <- max(dens_old$y, na.rm = TRUE)
      
      dens_new_list <- vector("list", length(logdf_list))
      
      for (x in 1) {
        
        new_vals <- full_logdf[, col_name]
        
        if (trait == "brdth") {
          new_vals <- exp(new_vals)
        }
        
        if (trait == "tol") {
          new_vals <- 0.05 + 0.95 * (exp(-1 * new_vals) / (1 + exp(-1 * new_vals)))
        }
        
        dens_new_list[[x]] <- density(new_vals)
        ymax <- max(ymax, dens_new_list[[x]]$y, na.rm = TRUE)
      }
      
      plot(
        c(xmin, xmax), c(0, ymax * 1.5),
        type = "n",
        main = mainlab,
        xlab = "",
        ylab = "Density"
      )
      
      for (x in seq_along(dens_new_list)) {
        lines(dens_new_list[[x]], col = 3)
      }
      
      lines(dens_old)
    }
  }
  
  dev.off()
  
}



#####compare old and new missing /missing and full mcmcs




tree <- ENA_Pleth_Tree
#pa_data <- read.csv("data/Pleth_data_vignette.csv")


miss=1


sp_col <- "species"
occ_col <- "PA"
env_preds <- c("bio12", "bio1")

Npred <- length(env_preds)
Ntips <- length(tree$tip.label)



sp_dirs = list.dirs("pletho_full_out_uninform_height_prior_news",recursive = F)
#miss=82
#sp_dir = sp_dirs[[miss]]



sp_files = list.files("pletho_full_out_uninform_height_prior_news", full.names = T)

sp_log_files = sp_files[grepl("pars.log",sp_files)]

logdf_list=list()

for(i in 1:length(sp_log_files)){
  
  logdf <- read_BePhyNE_log(
    file_name = paste0(sp_log_files[[i]])
  )
  print(nrow(logdf))
  logdf_list[[i]] = logdf
  
}


logdf_list_wburnin =  logdf_list

logdf_list = logdf_list_wburnin 


logdf_list= lapply(logdf_list, function(i) i[-(1:2500),])

logdf_list = logdf_list[unlist(lapply( logdf_list, nrow))>100]

new_full_logdf = do.call(rbind, logdf_list)

thin_rows <- function(x, target_n = 1000) {
  
  n <- nrow(x)
  
  if (n <= target_n) {
    return(x)
  }
  
  idx <- unique(round(seq(1, n, length.out = target_n)))
  
  x[idx, , drop = FALSE]
}

new_full_logdf = thin_rows(new_full_logdf)


saveRDS(new_full_logdf , file.path("pletho_full_compiled_logs", paste0("uninform_height_prior_news_compiled_full_species_logdf", ".RDS")))

rownames(new_full_logdf)

compiled_logdf = readRDS(file = file.path("pletho_miss_compiled_logs_height_uninform_prior_new", paste0("compiled_missing_species_logdf", ".RDS")))


new_miss_logdf = compiled_logdf

#saveRDS(new_full_logdf , "~/Downloads/full_logdf_new.RDS")
#saveRDS(new_miss_logdf , "~/Downloads/miss_logdf_new.RDS")

new_full_logdf = readRDS(file = file.path("pletho_full_compiled_logs/uninform_height_prior_news_compiled_full_species_logdf.RDS"))


#new_full_logdf = readRDS(new_full_logdf , "~/Downloads/full_logdf_new.RDS")
#new_miss_logdf = readRDS(new_miss_logdf , "~/Downloads/miss_logdf_new.RDS")


old_missing_data_mcmc = readRDS("~/Manuscripts/BePhyNE/R_code/presub_HPC_analyses_pub_FullandMiss/pletho_analyses/Josef_posterior_ecdf_miss_dat/missing_data_mcmc.RDS")
old_data_mcmc         = readRDS("~/Manuscripts/BePhyNE/R_code/presub_HPC_analyses_pub_FullandMiss/pletho_analyses/Josef_posterior_ecdf_miss_dat/full_data_mcmc.RDS")


#old_missing_data_mcmc = readRDS("~/Downloads/rstudio-export (33)/missing_data_mcmc.RDS")
#old_data_mcmc         = readRDS("~/Downloads/rstudio-export (33)/full_data_mcmc.RDS")

traits_out = readRDS("~/Manuscripts/BePhyNE/R_code/traits_out.RDS")


medians_full = traits_out[[1]]
medians_miss = traits_out[[3]]

new_median_full = lapply(medians_full, function(i) i*0)
new_median_miss = lapply(medians_full, function(i) i*0)


for(i in 1:2){
rownames(medians_miss[[i]]) = rownames(medians_full[[i]])
}

{
  
  trait_info <- data.frame(
    pred     = c(1, 2, 1, 2, 1, 2),
    trait    = c("opt", "opt", "brdth", "brdth", "tol", "tol"),
    old_idx  = c(1, 2, 3, 4, 5, 6),
    main_lab = c("C 1", "C 2", "W 1", "W 2", "H 1", "H 2"),
    stringsAsFactors = FALSE
  )
  
  {
    
    pdf(
      file = "~/Downloads/uninform_height_prior__miss_full_newprior_Trace_comp_pletho_all_species_traces_single_post.pdf",
      width = 10,
      height = 9
    )
    
    for (i in 1:82) {
      
      sp <- tree$tip.label[[i]]
      
      par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
      
      for (j in seq_len(nrow(trait_info))) {
        
        pred    <- trait_info$pred[j]
        trait   <- trait_info$trait[j]
        old_idx <- trait_info$old_idx[j]
        
        mainlab <- paste(trait_info$main_lab[j], sp)
        
        col_name <- paste0("pred_", pred, "_dat.", trait, "_", sp)
        
        old_full_vals <- old_data_mcmc[[old_idx]][, i]
        new_full_vals <- new_full_logdf[, col_name]
        old_miss_vals <- old_missing_data_mcmc[[old_idx]][[i]][, i]
        new_miss_vals <- new_miss_logdf[, col_name]
        
        
        if (trait == "brdth") {
          new_full_vals <- exp(new_full_vals)
          new_miss_vals <- exp(new_miss_vals)
          
        }
        
        if (trait == "tol") {
          new_full_vals <- 0.05 + 0.95 * (exp(-1 * new_full_vals) / (1 + exp(-1 * new_full_vals)))
          new_miss_vals <- 0.05 + 0.95 * (exp(-1 * new_miss_vals) / (1 + exp(-1 * new_miss_vals)))
          
        }
        
        ymin <- min(c( old_full_vals
                       ,new_full_vals
                       ,old_miss_vals
                       ,new_miss_vals), na.rm = TRUE)
        ymax <- max(c(old_full_vals
                      ,new_full_vals
                      ,old_miss_vals
                      ,new_miss_vals), na.rm = TRUE)
        
        xmax <- max(c(length(old_full_vals)
                      ,length(new_full_vals)
                      ,length(old_miss_vals)
                      ,length(new_miss_vals))
        )
        
        plot(
          c(0, xmax), c(ymin, ymax),
          type = "n",
          main = mainlab,
          xlab = "Index",
          ylab = ""
        )
        
        
        lines(seq_along(old_full_vals), old_full_vals, col = 1)
        lines(seq_along(new_full_vals), new_full_vals, col = 3)
        lines(seq_along(old_miss_vals), old_miss_vals, col = 2)
        lines(seq_along(new_miss_vals), new_miss_vals, col = 4)
        
        abline(h= median(old_full_vals), col = 1)
        abline(h= median(new_full_vals), col = 3)
        abline(h= median(old_miss_vals), col = 2)
        abline(h= median(new_miss_vals), col = 4)
        
        
      }
    }
    
    dev.off()
  }
  
  
  {
  pdf(
    file = "~/Downloads/new_only_uninform_height_prior__miss_full_newprior_Trace_comp_pletho_all_species_densities_single_post.pdf",
    width = 10,
    height = 9
  )
  
  for (i in 1:82) {
    
    sp <- tree$tip.label[[i]]
    
    par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
    
    for (j in seq_len(nrow(trait_info))) {
      
      pred    <- trait_info$pred[j]
      trait   <- trait_info$trait[j]
      old_idx <- trait_info$old_idx[j]
      
      mainlab <- paste(trait_info$main_lab[j], sp)
      
      col_name <- paste0("pred_", pred, "_dat.", trait, "_", sp)
      
      old_full_vals <- old_data_mcmc[[old_idx]][, i]
      new_full_vals <- new_full_logdf[, col_name]
      old_miss_vals <- old_missing_data_mcmc[[old_idx]][[i]][, i]
      new_miss_vals <- new_miss_logdf[, col_name]
      
      
      if (trait == "opt") {
        col_id=1
        new_full_vals <- (new_full_vals)
        new_miss_vals <- (new_miss_vals)
        
        old_median_full_trait = medians_full[[pred]][sp,1] 
        old_median_miss_trait = medians_miss[[pred]][sp,1] 
        
      }
      
      
      if (trait == "brdth") {
        col_id=2
        new_full_vals <- exp(new_full_vals)
        new_miss_vals <- exp(new_miss_vals)
        
        old_median_full_trait = medians_full[[pred]][sp,2] 
        old_median_miss_trait = medians_miss[[pred]][sp,2] 
        
      }
      
      if (trait == "tol") {
        col_id=3
        
        new_full_vals <- 0.05 + 0.95 * (exp(-1 * new_full_vals) / (1 + exp(-1 * new_full_vals)))
        new_miss_vals <- 0.05 + 0.95 * (exp(-1 * new_miss_vals) / (1 + exp(-1 * new_miss_vals)))
        
        old_median_full_trait = medians_full[[pred]][sp,3] 
        old_median_miss_trait = medians_miss[[pred]][sp,3] 
        
      }
      
      xmin <- min(c( old_full_vals
                     ,new_full_vals
                     ,old_miss_vals
                     ,new_miss_vals), na.rm = TRUE)
      xmax <- max(c(old_full_vals
                    ,new_full_vals
                    ,old_miss_vals
                    ,new_miss_vals), na.rm = TRUE)
      

      
      old_full_dens = density(old_full_vals)
      new_full_dens = density(new_full_vals)
      old_miss_dens = density(old_miss_vals)
      new_miss_dens = density(new_miss_vals)
      
      ymax= max(c(old_full_dens$y,
                  new_full_dens$y,
                  old_miss_dens$y,
                  new_miss_dens$y), na.rm = TRUE)
      

      xmax= max(c(old_full_dens$x,
                  new_full_dens$x,
                  old_miss_dens$x,
                  new_miss_dens$x), na.rm = TRUE)
      
      xmin= min(c(old_full_dens$x,
                  new_full_dens$x,
                  old_miss_dens$x,
                  new_miss_dens$x), na.rm = TRUE)
      

      
      plot(
        c(xmin, xmax), c(0, ymax * 1.5),
        type = "n",
        main = mainlab,
        xlab = "",
        ylab = "Density"
      )
      
      

      
      lines(old_full_dens, col = 1)
      lines(new_full_dens, col = 3)
      lines(old_miss_dens, col = 2)
      lines(new_miss_dens, col = 4)
      
      #abline(v= old_median_full_trait, col = 1, lwd=4.5)
      #abline(v= median(new_full_vals), col = 3)
      #abline(v=old_median_miss_trait, col = "red4", lwd=4.5)
      
      
      abline(v= median(old_full_vals), col = "black", lwd=3)
      abline(v= median(new_full_vals), col = 3, lwd=3)
      abline(v= median(old_miss_vals), col = 2, lwd=3)
      abline(v= median(new_miss_vals), col = 4, lwd=3)
      

      new_median_full[[pred]][sp,col_id] =  median(new_full_vals)
      new_median_miss[[pred]][sp,col_id] =  median(new_miss_vals)
      
    }
  }
  
  
  dev.off()
  }
}


new_median_full
new_median_miss[[1]][tree$tip.label,] -

miss_log_summary$median_parlist$traits[[1]][[1]][tree$tip.label,] 
#plot response curve medians relative to occurrence points


new_full_logdf






