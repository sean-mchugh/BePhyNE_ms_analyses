library(stringr)

library(coda)

#setwd("/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims")
collect_sim_medians_truevals = function(dirs, file_grep, names){
  
  job_id=NA
  percent_error_w_1=list()
  percent_error_c_1=list()
  percent_error_w_2=list()
  percent_error_c_2=list()
  percent_error_A_C_1=list()
  percent_error_A_C_2=list()
  percent_error_Rsd_C_1=list()
  percent_error_Rsd_C_2=list()
  percent_error_A_W_1=list()
  percent_error_A_W_2=list()
  percent_error_Rsd_W_1=list()
  percent_error_Rsd_W_2=list()
  percent_error_Rcor_1=list()
  percent_error_Rcor_2=list()
  
  MCMC_summary_list=list()
  
  full_traits_w     = list()
  full_traits_c     = list()
  full_traits_A_C   = list()
  full_traits_Rsd_C = list()
  full_traits_A_W   = list()
  full_traits_Rsd_W = list()
  full_traits_Rcor  = list()
  
  
  for (k in 1:length(file_grep)){
    
    
    
    traits_w_1=list()
    traits_c_1=list()
    traits_w_2=list()
    traits_c_2=list()
    traits_A_C_1=list()
    traits_A_C_2=list()
    traits_Rsd_C_1=list()
    traits_Rsd_C_2=list()
    traits_A_W_1=list()
    traits_A_W_2=list()
    traits_Rsd_W_1=list()
    traits_Rsd_W_2=list()
    traits_Rcor_1=list()
    traits_Rcor_2=list()
    
    traits_w    = list()
    traits_c    = list()
    traits_A_C  = list()
    traits_Rsd_C= list()
    traits_A_W  = list()
    traits_Rsd_W= list()
    traits_Rcor = list()
    
    mcmc_files= list.files(dirs[[k]])
    mcmc_files = mcmc_files[grep(file_grep[[k]],   mcmc_files)]
    
    MCMC_runs=list()
    for(i in 1:length( mcmc_files)){
    #for(i in 1:1){
      
      
      #load(paste(dirs[[k]],mcmc_files[[i]], sep=""))
      MCMC_run = readRDS(paste(dirs[[k]],mcmc_files[[i]], sep=""))
      MCMC_runs=MCMC_run[[1]]
      #names(MCMC_runs)[[i]]=str_split(mcmc_files[[i]],pattern ="_" )[[1]][[1]]
      
      print(i)
      #}
      #
      #
      #
      #
      #
      #MCMC_run=MCMC_runs
      x=i
      ### pul par estimates and true values from saved MCMC results ######
      #for (x in 1:length(MCMC_run)){
      
      #load(MCMC_files[[x]]) # load file
      MCMC_run_single=MCMC_runs
      
      mymcmc_trait_full     = MCMC_run_single$chains$mymcmc_trait_full
      mymcmc_trait_full_2   = MCMC_run_single$chains$mymcmc_trait_full_2
      mymcmc_A              = MCMC_run_single$chains$mymcmc_A
      mymcmc_A_2            = MCMC_run_single$chains$mymcmc_A_2
      mymcmc_R_sd           = MCMC_run_single$chains$mymcmc_R_sd
      mymcmc_R_sd_2         = MCMC_run_single$chains$mymcmc_R_sd_2
      mymcmc_R_corr         = MCMC_run_single$chains$mymcmc_R_corr
      mymcmc_R_corr_2       = MCMC_run_single$chains$mymcmc_R_corr_2
      
      TruePars_scale<-MCMC_run_single$TruePars_scale
      
      
      
      
      
      
      # setwd("/home/seanwm/GLM_MCMC/Plots/line_plots")
      
      traits_w_1[[x]]    <-do.call(rbind, lapply(1:ncol(mymcmc_trait_full$mymcmc.center), function(sp)  c(median(mymcmc_trait_full$mymcmc.width[,sp]),
                                                                                                          TruePars_scale$sim_dat$sim_dat_bt[[1]][[sp,2]])
      )
      )
      
      
      traits_c_1[[x]]    <- do.call(rbind, lapply(1:ncol(mymcmc_trait_full$mymcmc.center), function(sp)  c(median(mymcmc_trait_full$mymcmc.center[,sp]),
                                                                                                           TruePars_scale$sim_dat$sim_dat_bt[[1]][[sp,1]] )
      )
      )
      
      traits_w_2[[x]]    <-     do.call(rbind, lapply(1:ncol(mymcmc_trait_full$mymcmc.center), function(sp)  c(median(mymcmc_trait_full_2$mymcmc.width[,sp]),
                                                                                                               TruePars_scale$sim_dat$sim_dat_bt[[2]][[sp,2]])
      )
      )
      
      traits_c_2[[x]]    <-      do.call(rbind, lapply(1:ncol(mymcmc_trait_full$mymcmc.center), function(sp)  c(median(mymcmc_trait_full_2$mymcmc.center[,sp]),
                                                                                                                TruePars_scale$sim_dat$sim_dat_bt[[2]][[sp,1]] )
      )
      )
      
      traits_A_C_1[[x]]  <-        c(median(mymcmc_A$mymcmc_A.center),
                                     TruePars_scale$A$A_bt[[1]][[1]], HPDinterval(mymcmc_A$mymcmc_A.center)[1:2]  )
      
      
      traits_A_C_2[[x]]   <-        c(median(mymcmc_A_2$mymcmc_A.center),
                                      TruePars_scale$A$A_bt[[2]][[1]] , HPDinterval(mymcmc_A_2$mymcmc_A.center)[1:2] )
      
      
      traits_Rsd_C_1[[x]] <-        c(median(mymcmc_R_sd$mymcmc_sd.center),
                                      TruePars_scale$R$R_sd[[1]][[1]]  , HPDinterval(mymcmc_R_sd$mymcmc_sd.center)[1:2] )
      
      
      traits_Rsd_C_2[[x]] <-        c(median(mymcmc_R_sd_2$mymcmc_sd.center),
                                      TruePars_scale$R$R_sd[[2]][[1]] , HPDinterval(mymcmc_R_sd_2$mymcmc_sd.center)[1:2] )
      
      
      traits_A_W_1[[x]]   <-        c(median(mymcmc_A$mymcmc_A.width),
                                      TruePars_scale$A$A_bt[[1]][[2]] ,  HPDinterval(mymcmc_A$mymcmc_A.width)[1:2])
      
      
      traits_A_W_2[[x]]   <-        c(median(mymcmc_A_2$mymcmc_A.width),
                                      TruePars_scale$A$A_bt[[2]][[2]] , HPDinterval(mymcmc_A_2$mymcmc_A.width)[1:2] )
      
      
      traits_Rsd_W_1[[x]] <-        c(median(mymcmc_R_sd$mymcmc_sd.width),
                                      TruePars_scale$R$R_sd[[1]][[2]] , HPDinterval(mymcmc_R_sd$mymcmc_sd.width)[1:2])
      
      traits_Rsd_W_2[[x]] <-        c(median(mymcmc_R_sd_2$mymcmc_sd.width),
                                      TruePars_scale$R$R_sd[[2]][[2]] , HPDinterval(mymcmc_R_sd_2$mymcmc_sd.width)[1:2])
      
      traits_Rcor_1[[x]]   <-     c(median(mymcmc_R_corr$mymcmc_c.w),
                                    TruePars_scale$R$R_cor[[1]][1,2], HPDinterval(mymcmc_R_corr$mymcmc_c.w)[1:2])
      
      
      traits_Rcor_2[[x]]   <-     c(median(mymcmc_R_corr_2$mymcmc_c.w),
                                    TruePars_scale$R$R_cor[[2]][1,2], HPDinterval(mymcmc_R_corr_2$mymcmc_c.w)[1:2])
      
      
      traits_w    [[x]] = rbind(traits_w_1[[x]]   , traits_w_2[[x]] )
      traits_c    [[x]] = rbind(traits_c_1[[x]]   , traits_c_2[[x]] )
      traits_A_C  [[x]] = rbind(traits_A_C_1[[x]], traits_A_C_2[[x]])
      traits_Rsd_C[[x]] = rbind(traits_Rsd_C_1[[x]], traits_Rsd_C_2[[x]])
      traits_A_W  [[x]] = rbind(traits_A_W_1[[x]], traits_A_W_2[[x]])
      traits_Rsd_W[[x]] = rbind(traits_Rsd_W_1[[x]], traits_Rsd_W_2[[x]])
      traits_Rcor [[x]] = rbind(traits_Rcor_1[[x]], traits_Rcor_2[[x]])
      
      
      
    }
    
    full_traits_w    [[ k ]]= do.call(rbind,traits_w    )[,1:2]
    full_traits_c    [[ k ]]= do.call(rbind,traits_c    )[,1:2]
    full_traits_A_C  [[ k ]]= do.call(rbind,traits_A_C  )[,1:2]
    full_traits_Rsd_C[[ k ]]= do.call(rbind,traits_Rsd_C)[,1:2]
    full_traits_A_W  [[ k ]]= do.call(rbind,traits_A_W  )[,1:2]
    full_traits_Rsd_W[[ k ]]= do.call(rbind,traits_Rsd_W)[,1:2]
    full_traits_Rcor [[ k ]]= do.call(rbind,traits_Rcor )[,1:2]
    
    names(full_traits_w    )[[k]] = names[[k]]
    names(full_traits_c    )[[k]] = names[[k]]
    names(full_traits_A_C  )[[k]] = names[[k]]
    names(full_traits_Rsd_C)[[k]] = names[[k]]
    names(full_traits_A_W  )[[k]] = names[[k]]
    names(full_traits_Rsd_W)[[k]] = names[[k]]
    names(full_traits_Rcor )[[k]] = names[[k]]
    
  }
  
  simmedian_vs_true_pars = list( full_traits_w    = full_traits_w    
                                 ,full_traits_c    = full_traits_c    
                                 ,full_traits_A_C  = full_traits_A_C  
                                 ,full_traits_Rsd_C= full_traits_Rsd_C
                                 ,full_traits_A_W  = full_traits_A_W  
                                 ,full_traits_Rsd_W= full_traits_Rsd_W
                                 ,full_traits_Rcor = full_traits_Rcor )
  
  return(simmedian_vs_true_pars)
}  

#
dirs     = c("Background_Abs_runs_grid/", "Background_Abs_runs_grid/")
file_grep=c ("FALSE", "TRUE")
names    = c("Only_True_Abs", "Background_Abs")

simmedian_vs_true_pars = lapply(1:length(dirs ), function(treatment) collect_sim_medians_truevals(dirs[[treatment]], file_grep[[treatment]], names[[treatment]]))
saveRDS(simmedian_vs_true_pars,"grid_Background_simmedian_vs_true_pars.RDS")
  
  
dirs     = rep("Kappa_runs/", 5)
file_grep     =paste0("kappa_", c("0_", "0.25", "0.5", "0.75", "1"))
names         = paste0("kappa_", c("0", "0.25", "0.5", "0.75", "1"))



simmedian_vs_true_pars = lapply(1:length(dirs ), function(treatment) collect_sim_medians_truevals(dirs[[treatment]], file_grep[[treatment]], names[[treatment]]))


saveRDS(simmedian_vs_true_pars,"Kappa_simmedian_vs_true_pars.RDS")


dirs          = c("Kappa_runs/", "Lambda_runs/")
file_grep     = c("kappa_1", "lambda_0" )
names         = c("lambda_1","lambda_0" )

simmedian_vs_true_pars = lapply(1:length(dirs ), function(treatment) collect_sim_medians_truevals(dirs[[treatment]], file_grep[[treatment]], names[[treatment]]))
saveRDS(simmedian_vs_true_pars,"Lambda_simmedian_vs_true_pars.RDS")



dirs     = c("Kappa_runs_new_oldPA/")
file_grep=c ("iter")
names    = c("Only_True_Abs")

simmedian_vs_true_pars = lapply(1:length(dirs ), function(treatment) collect_sim_medians_truevals(dirs[[treatment]], file_grep[[treatment]], names[[treatment]]))
saveRDS(simmedian_vs_true_pars,"Kappa_oldpa.RDS")

  
