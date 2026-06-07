
#devtools::install_github( "sean-mchugh/BePhyNE", lib = "packages", upgrade = "never", dependencies = FALSE )

#devtools::install_github( "sean-mchugh/BePhyNE", upgrade = "never", force = T, dependencies = FALSE )


#library(BePhyNE)
#
##other packages
#library(treeplyr)
#require(readr)
#require(tibble)
#require(MCMCpack)
#require(devtools)
#require(coda)
#require(devtools)
#require(coda)
#require(ape)
#require(truncdist)
#require(geiger)
#require(phytools)
#require(tidyr)
#require(ratematrix)
#require(mvMORPH)
#require(dplyr)
#require(mvtnorm)
#require(MultiRNG)
#require(Rphylopars)
#require(Rcpp)
#require(doParallel)
#require(corpcor)
#require(Matrix)
#require(treeplyr)
#require(bayou)
#require(crayon)
#require(shape)
#require(scales)
#require(phytools)
#require(robustbase)

#setwd("/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims")

library(MCMCpack, lib="packages/")
library(robustbase, lib="packages/")
library(caTools, lib="packages/")

library(flux, lib="packages/")
library(evd, lib="packages/")
library(truncdist, lib="packages/")
library(MultiRNG, lib="packages/")
library(Rphylopars, lib="packages/")
library(BePhyNE, lib="packages/")

#setwd("/storage1/fs1/michael.landis/Active/BePhyNE_sims/Height_runs/")
#run_id=1
reps=1
#height_mean = 0.7
#height_sd   = 1.00


args <- commandArgs(trailingOnly = TRUE)

run_id      <- as.integer(args[1])
height_mean <- as.numeric(args[2])
height_sd   <- as.numeric(args[3])

cat("run_id:"     , run_id, "\n"     )
cat("height_mean:", height_mean, "\n")
cat("height_sd:"  , height_sd, "\n"  )


height_mean_label = gsub(pattern =".", replacement = "p", height_mean, fixed = T )
height_sd_label   = gsub(pattern =".", replacement = "p", height_sd  , fixed = T )

outdir = paste0("Neww_New_New_Height_runs/Height_runs_mean_",height_mean_label, "_sd_",height_sd_label)
#dir.create("New_New_Height_runs")
dir.create(outdir )

{
  
  {
    
    {
      #setwd(" ~/../../media/HDSSD1/shared/swm_outfiles")
      
      getwd()
      #dir<- "/home/seanwm/GLM_MCMC/out_files"
      #outname="single_R_fast"
      
      ID <- run_id #paste( sample(x=1:9, size=9, replace=TRUE), collapse="")
      
      
      set.seed(ID)
      
      H_fixed=F
      
      sparse_sp=F
      
      
      simMiss=F
      nMiss=0
      tips=200
      #n=tips/10
      
      norm = T
      unif = F
      
      
      pred=2
      Npred=pred
      batches=1
      
      set_prevalence=F
      prev=.1
      grid_size<-1500
      #N<-1000
      span=1.5
      
      PA_plot=F
      
      
      v=c(.05, .05)
      k=lapply(1:length(v), function(x) log(v[[x]]/(1-v[[x]])))
      
      
      kappa_scale=0
      
      
      phylo_unscaled <-pbtree(n=tips,scale=1) 
      
      phylo <-  phylo_unscaled # lapply(1:reps, function(x) rescale(  phylo_unscaled[[x]] , "lambda", kappa_scale))
      
      
      #tree<-geiger::rescale(ENA_Pleth_Tree, "kappa",0)
      
      
      #max(branching.times(phylo[[1]]))
      
      ####simulate environmental data######
      
      
      xmin_1=10
      xmax_1=30
      
      xmin_2=10
      xmax_2=30
      
      xdataamount_1=50
      xdataamount_2=50
      
      #X1 =sort(runif(xdataamount_1, min=xmin_1, max=xmax_1))
      #X2 =sort(runif(xdataamount_2, min=xmin_2, max=xmax_2))
      
      #simuklate real evirnmental axes
      
      X1 <- seq(xmin_1,xmax_1, length.out=1000)
      X2 <- seq(xmin_2,xmax_2, length.out=1000)
      
      
      #scale axes
      scaled_X1<-scale(X1)
      scaled_X2<-scale(X2)
      
      
      
      scale_atr <- list(scale=list(attributes(scaled_X1)$`scaled:scale`,
                                   attributes(scaled_X2)$`scaled:scale`) ,
                        center= list(attributes(scaled_X1)$`scaled:center`,
                                     attributes(scaled_X2)$`scaled:center`))
      
      
      #full_X1<-seq(-4,4 , length.out=xdataamount_1)
      #full_X2<-seq(-4,4 , length.out=xdataamount_1)
      
      #simulate data on this axes (normally distributed )
      
      full_X1<-sort(rnorm(xdataamount_1, mean=0, sd=2))
      full_X2<-full_X1
      
      (full_X1*scale_atr$scale[[1]])+scale_atr$center[[1]]
      (full_X2*scale_atr$scale[[2]])+scale_atr$center[[2]]
      
      
      log(range(scaled_X1)[[2]])
      
      .5^2
      
      
    }
    ###set scaled prior#####
    
      #norm priors
      
      
     #par.mu.norm.H.fixed<-list( matrix(c(  0,  .2, #center(norm)
     #                                      log(0.3), .1), #width (lnorm)
     #                                  nrow = 2, ncol = 2, byrow = TRUE) #pred 1 root prior (center mean, center sd, width mean, width sd)
     #                           , matrix(c(0, .2,
     #                                      log(0.3), .1),
     #                                    nrow = 2, ncol = 2, byrow = TRUE) #pred 2 root prior (center mean, center sd, width mean, width sd)
     #)
     #par.sd.lnorm.meanlog<-list(c(log(.1), log(.1)) #pred 1 (center, width, height)
     #                           ,c(log(.1), log(.1)) #pred 2
     #)
     #par.sd.lnorm.sdlog<-list(c(.5, .5) #pred 1
     #                         ,c(.5, .5) #pred 2
     #)
     #
     #par.sd.lnorm <-lapply(1:length(par.sd.lnorm.meanlog), function(x) cbind(par.sd.lnorm.meanlog[[x]], par.sd.lnorm.sdlog[[x]]))
     #par.sd.lnorm <-lapply(1:length(par.sd.lnorm.meanlog), function(x) cbind(par.sd.lnorm.meanlog[[x]], par.sd.lnorm.sdlog[[x]]))
     #
     #Prior_scale <-lapply(1:length(par.sd.lnorm.meanlog), function(x) makePrior_ENE(r = 2, p = 1, den.mu = "norm", heights_by_sp = sample(height_mean, size = tips, replace=T), heights_sd = height_sd, par.mu = par.mu.norm.H.fixed[[x]], den.sd = "lnorm" , par.sd = par.sd.lnorm[[x]], plot=T))
      
      
      
      Prior_scale <- make_all_priors(
        Npred,
        length(phylo$tip.label),
        bd_range = NA,
        
        ## root priors
        root_opt_mean      = rep(0, Npred),
        root_opt_sd        = rep(0.2, Npred),
        root_brdth_meanlog = rep(log(0.3), Npred),
        root_brdth_sdlog   = rep(0.1, Npred),
        
        ## log-normal hyper-priors
        sigsq_opt_meanlog   = rep(log(.1), Npred),
        sigsq_opt_sdlog     = rep(0.5, Npred),
        sigsq_brdth_meanlog = rep(log(.1), Npred),
        sigsq_brdth_sdlog   = rep(0.5, Npred),
        
        ## heights
        heights_sd_by_sp = height_sd,
        heights_mean_by_sp = height_mean,
        use_glm_height_mean = F,
        #species_data = data_final,
        #uninf_height_sd = 1.00,
        tree = phylo,#sample(0.95, size = Ntips, replace = TRUE),
        
        
        ## constants forwarded to makePrior_ENE
        r = Npred,
        p = 1,
        plot = TRUE,
        miss_data_uninf = F
      )
      
      
      Prior_scale[[1]]$pars$heights
      
      
      TruePars_scale <- priorSim_pars(Prior=Prior_scale, hard_coded_heights = NULL, phylo = phylo, dist="norm")
      #startPars_scaled <-lapply(1:reps, function(x) priorSim_pars(Prior_scale, hard_coded_heights = c(0.95, 0.95), phylo[[x]], dist="norm") )
      TruePars_scale$sim_dat$sim_dat_bt
      #TruePars_scale$sim_dat$sim_dat_ft
      
      
      job_id = as.numeric(Sys.getenv("SLURM_JOB_ID"))
      
      cat("Starting run:",job_id, ID, "\n")
      
      
      #dir<-getwd()
      
      
      
      #dir.create(dir)
      
    
    
    
    
    #presence absences data then simulated from scaled traits and beta coefficients
    
    pa_data_scaled<-simPA( res=TruePars_scale$sim_dat$sim_dat_bt,
                                                      tree=phylo, span=span,
                                                      grid_size=grid_size,
                                                      simMiss=simMiss,
                                                      nMiss=nMiss
    
    )
    
    
    startPars_scaled=list()

    ###set tunning parameters and Run#####
    center_slide=.18
    center_mult=.12
    width_slide=.15
    width_mult=.2
    height_slide=.5
    height_mult=.3
    
    move_details = make_tuning(tree = phylo, pred = 2)
    
    sparse_sp=F
    iterations=1000000
    trim_freq=iterations/10000
    chain_end=(iterations-(iterations/10))/trim_freq
    
    job_id = as.numeric(Sys.getenv("SLURM_JOB_ID"))
    #job_rep = as.numeric(Sys.getenv("PARALLEL_JOBSLOT"))
    
    #ID <- paste( sample(x=1:9, size=9, replace=TRUE), collapse="")
    cat("Starting run:",job_id, ID, "\n")
    
    
    i=reps
    
    
    for (i in 1:reps){
      
      
      rep=0
      #repeat{
        print("rep")
        print(rep)
        rep=rep+1
        #print
        #startPars_scaled[[i]] <- priorSim_pars(Prior_scale, phylo_unscaled[[i]], dist="norm",hard_coded_heights = NULL)
        #break when no tip has a -inf likelihood
        startPars = get_starting_values(Prior_scale = Prior_scale, tree=phylo,  data =   pa_data_scaled[[1]])
        
        #startPars$R
        #
        #cbind( startPars$sim_dat$sim_dat_bt[[1]][,3]
        #       , TruePars_scale$sim_dat$sim_dat_bt[[1]][,3])
        #
        #test=try(  metro_haste_full_MV(       R_corr_start = startPars_scaled[[i]]$R$R_cor
        #                                      , R_sd_start   = startPars_scaled[[i]]$R$R_sd
        #                                      , A_start      = startPars_scaled[[i]]$A$A_bt
        #                                      , Prior_scale
        #                                      , phylo_unscaled[[i]] #testing how well a kappa dataset handles under our model when we can observe the speciation driven patterns
        #                                      , tibble_data = startPars_scaled[[i]]$sim_dat$sim_td_bt
        #                                      , pa_data     = pa_data_scaled[[i]]$species_data
        #                                      , iterations=100
        #                                      , burnin=1
        #                                      , move_prob=move_prob
        #                                      , n=2
        #                                      , print.i.freq=10
        #                                      , print.ac.freq=10
        #                                      , printing=TRUE
        #                                      , trim=T
        #                                      , trim_freq=trim_freq
        #                                      , H_fixed=F
        #                                      , tuning=tuning
        #                                      , center_fixed=F
        #                                      , write_file=F
        #                                      #, IDlen=5
        #                                      #, dir=NA
        #                                      #, outname=NA
        #                                      , prior_only= F
        #                                      , glm_only  = F
        #                                      , plot=F
        #                                      , plot_freq=iterations/5
        #                                      , plot_file= NA
        #                                      , True_pars = NA
        #                                      , k=k
        #)
        #)
        
        saveRDS(TruePars_scale, file = paste0(outdir,"/True_",ID,".rds"))
        
        
        results_set_test=BePhyNE::BePhyNE_MCMC(tree = phylo_unscaled
                                              ,pa_data = pa_data_scaled$species_data
                                              ,Prior_scale = Prior_scale
                                              ,startPars = startPars
                                              ,move_details = move_details
                                              ,iterations = 10000
                                              ,trim_freq = 1
                                              ,write2file = F
                                              ,filename = paste0(outdir,ID,".log")
                                          
        )
        
        results_set_test$accept_ratios
        
        if(all(results_set_test$chain[[1]][[1]][[1]]$dat[,3] == results_set_test$chain[[1000]][[1]][[1]]$dat[,3])){
          
          
          xxxxxxxxx
        }
        
        #results_set_test$chain[[1]][[1]][[2]]$dat == results_set_test$chain[[1000]][[1]][[2]]$dat
        
        #if(length(findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data=pa_data_scaled[[i]]$species_data, plot=F))==0 && length(test)>1){
        #  print("good start")
        #  break}
        
        
     # }
      #findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data= pa_data_scaled[[i]]$species_data, plot=F)
      
      print("done")
    }
    
    
    #for (i in 1:reps){
    print("initializing MCMC")
    
    #startPars = get_starting_values(Prior_scale = Prior_scale, tree=phylo,  data =   pa_data_scaled[[1]])
    
   # registerDoParallel(cores=reps)
    #results_set<- foreach(i=1:reps) %dopar% {
    i=1  
      #```{r,message=FALSE, warning=F, results= "hide" }
      results_set=BePhyNE::BePhyNE_MCMC(tree          = phylo_unscaled
                                        ,pa_data      = pa_data_scaled$species_data
                                        ,Prior_scale  = Prior_scale
                                        ,startPars    = startPars
                                        ,move_details = move_details
                                        ,iterations   = iterations
                                        ,trim_freq    = trim_freq
                                        ,write2file   = T
                                        ,filename     =  paste0(outdir,"/Inf_",ID,".log")
                                  
                                          )
      
      
      
      
      
      #metro_haste_full_MV(          R_corr_start = startPars_scaled[[i]]$R$R_cor
      #                              , R_sd_start   = startPars_scaled[[i]]$R$R_sd
      #                              , A_start      = startPars_scaled[[i]]$A$A_bt
      #                              , Prior_scale
      #                              , phylo_unscaled[[i]] #testing how well a kappa dataset handles under our model when we can observe the speciation driven patterns
      #                              , tibble_data = startPars_scaled[[i]]$sim_dat$sim_td_bt
      #                              , pa_data=pa_data_scaled[[i]]$species_data
      #                              , iterations=iterations
      #                              , burnin=iterations/10
      #                              , move_prob=move_prob
      #                              , n=2
      #                              , print.i.freq=1000
      #                              , print.ac.freq=100
      #                              , printing=TRUE
      #                              , trim=T
      #                              , trim_freq=trim_freq
      #                              , H_fixed=T
      #                              , tuning=tuning
      #                              , center_fixed=F
      #                              , write_file=T
      #                              ,filename = paste0("Height_p95_runs/",ID,".log")
      #                              #, IDlen=5
      #                              #, dir=NA
      #                              #, outname=NA
      #                              , prior_only= F
      #                              , glm_only  = F
      #                              , plot=F
      #                              , plot_freq=iterations/5
      #                              , plot_file= NA
      #                              , True_pars = NA
      #                              , k=k
      #)
      #```
    
    
    #saveRDS(results, file=paste("BePhyNE_Kappa_run",toString(sample(1:1000,1)), sep="_" ))
    
    ###check likelihoods and acceptance ratios#######
    
    
    #dir<-paste("~/R/",job_id,sep="")
    
    #dir.create(dir)
    
    #dir.create(paste(dir,"/AR",sep=""))
    
    
   # MCMC_run=list()
   # 
   # 
   # 
   # 
   # for (i in 1:reps){
   #   
   #   results=results_set[[i]]
   #   
   #   mymcmc_trait_full     =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="trait", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   mymcmc_trait_full_2   =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="trait", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   mymcmc_A              =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="theta", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   mymcmc_A_2            =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="theta", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   mymcmc_R_sd           =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="R_sd", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   mymcmc_R_sd_2         =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="R_sd", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   mymcmc_R_corr         =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="R_corr", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   mymcmc_R_corr_2       =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="R_corr", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
   #   
   #   
   #   
   #   
   #   
   #   
   #   
   #   MCMC_run[[i]] <- list(chains=list(mymcmc_trait_full  =mymcmc_trait_full
   #                                     ,mymcmc_trait_full_2=mymcmc_trait_full_2
   #                                     ,mymcmc_A           =mymcmc_A
   #                                     ,mymcmc_A_2         =mymcmc_A_2
   #                                     ,mymcmc_R_sd        =mymcmc_R_sd
   #                                     ,mymcmc_R_sd_2      =mymcmc_R_sd_2
   #                                     ,mymcmc_R_corr      =mymcmc_R_corr
   #                                     ,mymcmc_R_corr_2    =mymcmc_R_corr_2
   #   ),
   #   TruePars_scale=TruePars_scale[[i]])
   #   
   #   
   #   
   #   
   #   
   #   #dir.create(paste(dir,"/MCMC_out",sep=""))
   #   
   #   
   #   
   # }
   # 
   # outputName=paste(ID, "Height_iter", iterations, "MCMC_run",".RData",sep="_")
   # outputPath=file.path(dir,outputName)
   # save("MCMC_run",file=outputPath)
      
      
      logdf = read_BePhyNE_log(paste0(paste0(outdir,"/Inf_",ID,".log.pars.log") ))
      logdf_summary = summarize_logdf(logdf)
      
      #TruePars_scale = readRDS(paste0(outdir,"/", true_files_list[[k]] ))
      
      
      
      
      
      traits_h_1[[k]] <-do.call(rbind, lapply(1:nrow(logdf_summary$median_parlist$traits[[1]][[1]]), function(sp)  c(logdf_summary$median_parlist$traits[[1]][[1]][sp,3],
                                                                                                                     TruePars_scale$sim_dat$sim_dat_bt[[1]][[sp,3]],logdf_summary$HPDlower_parlist$traits[[1]][[1]][sp,3], logdf_summary$HPDupper_parlist$traits[[1]][[1]][sp,3])
      )
      )
      
      
      
      traits_h_2[[k]] <-do.call(rbind, lapply(1:nrow(logdf_summary$median_parlist$traits[[2]][[1]]), function(sp)  c(logdf_summary$median_parlist$traits[[2]][[1]][sp,3],
                                                                                                                     TruePars_scale$sim_dat$sim_dat_bt[[2]][[sp,3]],logdf_summary$HPDlower_parlist$traits[[2]][[1]][sp,3], logdf_summary$HPDupper_parlist$traits[[2]][[1]][sp,3])
      )
      )
      
      
      
      
      
      traits_c_1[[k]] <-do.call(rbind, lapply(1:nrow(logdf_summary$median_parlist$traits[[1]][[1]]), function(sp)  c(logdf_summary$median_parlist$traits[[1]][[1]][sp,1],
                                                                                                                     TruePars_scale$sim_dat$sim_dat_bt[[1]][[sp,1]],logdf_summary$HPDlower_parlist$traits[[1]][[1]][sp,1], logdf_summary$HPDupper_parlist$traits[[1]][[1]][sp,1])
      )
      )
      
      
      
      traits_c_2[[k]] <-do.call(rbind, lapply(1:nrow(logdf_summary$median_parlist$traits[[2]][[1]]), function(sp)  c(logdf_summary$median_parlist$traits[[2]][[1]][sp,1],
                                                                                                                     TruePars_scale$sim_dat$sim_dat_bt[[2]][[sp,1]],logdf_summary$HPDlower_parlist$traits[[2]][[1]][sp,1], logdf_summary$HPDupper_parlist$traits[[2]][[1]][sp,1])
      )
      )
      
      traits_A_C_1[[k]]   <-        c(logdf_summary$median_parlist$A[[1]][[1]][[1]],
                                      TruePars_scale$A$A_bt[[1]][[1]], logdf_summary$HPDlower_parlist$A[[1]][[1]][[1]], logdf_summary$HPDupper_parlist$A[[1]][[1]][[1]] )
      
      
      traits_A_C_2[[k]]   <-        c(logdf_summary$median_parlist$A[[2]][[1]][[1]],
                                      TruePars_scale$A$A_bt[[2]][[1]], logdf_summary$HPDlower_parlist$A[[2]][[1]][[1]], logdf_summary$HPDupper_parlist$A[[2]][[1]][[1]] )
      
      
      traits_Rsd_C_1[[k]] <-        c(logdf_summary$median_parlist$Rsd[[1]][[1]][[1]],
                                      TruePars_scale$R$R_sd[[1]][[1]], logdf_summary$HPDlower_parlist$A[[1]][[1]][[1]], logdf_summary$HPDupper_parlist$A[[1]][[1]][[1]] )
      
      
      traits_Rsd_C_2[[k]] <-        c(logdf_summary$median_parlist$Rsd[[2]][[1]][[1]],
                                      TruePars_scale$R$R_sd[[2]][[1]], logdf_summary$HPDlower_parlist$Rsd[[2]][[1]][[1]], logdf_summary$HPDupper_parlist$Rsd[[2]][[1]][[1]] )
      
      
      traits_w_1[[k]] <-do.call(rbind, lapply(1:nrow(logdf_summary$median_parlist$traits[[1]][[1]]), function(sp)  c(logdf_summary$median_parlist$traits[[1]][[1]][sp,2],
                                                                                                                     TruePars_scale$sim_dat$sim_dat_bt[[1]][[sp,2]],logdf_summary$HPDlower_parlist$traits[[1]][[1]][sp,2], logdf_summary$HPDupper_parlist$traits[[1]][[1]][sp,2])
      )
      )
      
      
      traits_w_2[[k]] <-do.call(rbind, lapply(1:nrow(logdf_summary$median_parlist$traits[[2]][[1]]), function(sp)  c(logdf_summary$median_parlist$traits[[2]][[1]][sp,2],
                                                                                                                     TruePars_scale$sim_dat$sim_dat_bt[[2]][[sp,2]],logdf_summary$HPDlower_parlist$traits[[2]][[1]][sp,2], logdf_summary$HPDupper_parlist$traits[[2]][[1]][sp,2])
      )
      )
      
      traits_A_W_1[[k]]   <-        c(logdf_summary$median_parlist$A[[1]][[1]][[2]],
                                      TruePars_scale$A$A_bt[[1]][[2]], logdf_summary$HPDlower_parlist$A[[1]][[1]][[2]], logdf_summary$HPDupper_parlist$A[[1]][[1]][[2]] )
      
      
      traits_A_W_2[[k]]   <-        c(logdf_summary$median_parlist$A[[2]][[1]][[2]],
                                      TruePars_scale$A$A_bt[[2]][[2]], logdf_summary$HPDlower_parlist$A[[2]][[1]][[2]], logdf_summary$HPDupper_parlist$A[[2]][[1]][[2]] )
      
      
      traits_Rsd_W_1[[k]] <-        c(logdf_summary$median_parlist$Rsd[[1]][[1]][[2]],
                                      TruePars_scale$R$R_sd[[1]][[2]], logdf_summary$HPDlower_parlist$A[[1]][[1]][[2]], logdf_summary$HPDupper_parlist$A[[1]][[1]][[2]] )
      
      
      traits_Rsd_W_2[[k]] <-        c(logdf_summary$median_parlist$Rsd[[2]][[1]][[2]],
                                      TruePars_scale$R$R_sd[[2]][[2]], logdf_summary$HPDlower_parlist$Rsd[[2]][[1]][[2]], logdf_summary$HPDupper_parlist$Rsd[[2]][[1]][[2]] )
      
      
      
      traits_Rcor_1[[k]] <-     c(logdf_summary$median_parlist$Rcor[[1]][[1]][1,2],
                                  TruePars_scale$R$R_cor[[1]][1,2], logdf_summary$HPDlower_parlist$Rcor[[1]][[1]][1,2],logdf_summary$HPDupper_parlist$Rcor[[1]][[1]][1,2])
      
      
      traits_Rcor_2[[k]] <-     c(logdf_summary$median_parlist$Rcor[[2]][[1]][1,2],
                                  TruePars_scale$R$R_cor[[2]][1,2], logdf_summary$HPDlower_parlist$Rcor[[2]][[1]][1,2],logdf_summary$HPDupper_parlist$Rcor[[2]][[1]][1,2])
      
      
      
      
      
      
      
  }
  
  
  
  ## make full par lists#####
  full_traits_h_1    = do.call(rbind, lapply(1:length(traits_c_1), function(x)     traits_h_1[[x]]     ))
  full_traits_h_2    = do.call(rbind, lapply(1:length(traits_c_2), function(x)     traits_h_2[[x]]     ))
  
  full_traits_w_1    = do.call(rbind, lapply(1:length(traits_w_1), function(x)     traits_w_1[[x]]     ))
  full_traits_c_1    = do.call(rbind, lapply(1:length(traits_c_1), function(x)     traits_c_1[[x]]     ))
  full_traits_w_2    = do.call(rbind, lapply(1:length(traits_w_2), function(x)     traits_w_2[[x]]     ))
  full_traits_c_2    = do.call(rbind, lapply(1:length(traits_c_2), function(x)     traits_c_2[[x]]     ))
  full_traits_A_C_1  = do.call(rbind, lapply(1:length(traits_A_C_1), function(x)   traits_A_C_1[[x]]   ))
  full_traits_A_C_2  = do.call(rbind, lapply(1:length(traits_A_C_2), function(x)   traits_A_C_2[[x]]   ))
  full_traits_Rsd_C_1= do.call(rbind, lapply(1:length(traits_Rsd_C_1), function(x) traits_Rsd_C_1[[x]] ))
  full_traits_Rsd_C_2= do.call(rbind, lapply(1:length(traits_Rsd_C_2), function(x) traits_Rsd_C_2[[x]] ))
  full_traits_A_W_1  = do.call(rbind, lapply(1:length(traits_A_W_1), function(x)   traits_A_W_1[[x]]   ))
  full_traits_A_W_2  = do.call(rbind, lapply(1:length(traits_A_W_2), function(x)   traits_A_W_2[[x]]   ))
  full_traits_Rsd_W_1= do.call(rbind, lapply(1:length(traits_Rsd_W_1), function(x) traits_Rsd_W_1[[x]] ))
  full_traits_Rsd_W_2= do.call(rbind, lapply(1:length(traits_Rsd_W_2), function(x) traits_Rsd_W_2[[x]] ))
  full_traits_Rcor_1 = do.call(rbind, lapply(1:length(traits_Rcor_1), function(x)  traits_Rcor_1[[x]]  ))
  full_traits_Rcor_2 = do.call(rbind, lapply(1:length(traits_Rcor_2), function(x)  traits_Rcor_2[[x]]  ))
  
  
  MCMC_summary_list=list( 
    "full_traits_h_1"    =full_traits_h_1    ,
    "full_traits_w_1"    = full_traits_w_1    ,
    "full_traits_c_1"    =full_traits_c_1    ,
    "full_traits_h_2"    =full_traits_h_2    ,
    "full_traits_w_2"    =full_traits_w_2    ,
    "full_traits_c_2"    =full_traits_c_2    ,
    "full_traits_A_C_1"  =full_traits_A_C_1  ,
    "full_traits_A_C_2"  =full_traits_A_C_2  ,
    "full_traits_Rsd_C_1"=full_traits_Rsd_C_1,
    "full_traits_Rsd_C_2"=full_traits_Rsd_C_2,
    "full_traits_A_W_1"  =full_traits_A_W_1  ,
    "full_traits_A_W_2"  =full_traits_A_W_2  ,
    "full_traits_Rsd_W_1"=full_traits_Rsd_W_1,
    "full_traits_Rsd_W_2"=full_traits_Rsd_W_2,
    "full_traits_Rcor_1"  =full_traits_Rcor_1 ,
    "full_traits_Rcor_2 "=full_traits_Rcor_2 )
  
  
  saveRDS(MCMC_summary_list, paste(outdir,"/", ID, "_MCMC_summary_list.rds",  sep=""))
  
    
  
  
}
