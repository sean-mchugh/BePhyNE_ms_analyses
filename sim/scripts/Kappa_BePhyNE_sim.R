#library(BePhyNE)
#library(devtools)
##devtools::install_github("sean-mchugh/BePhyNE")
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

setwd("/storage1/fs1/michael.landis/Active/Sean/BePhyNE_ms_analyses/sim")

library(MCMCpack, lib="packages/")
library(robustbase, lib="packages/")

library(evd, lib="packages/")
library(truncdist, lib="packages/")
library(MultiRNG, lib="packages/")
library(Rphylopars, lib="packages/")
library(BePhyNE, lib="packages/")

#widths under this PA function underestimated, 
#setwd("/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/Kappa_runs_new_oldPA/")

{
  
  {
    
    {
      #setwd(" ~/../../media/HDSSD1/shared/swm_outfiles")
      
      getwd()
      #dir<- "/home/seanwm/GLM_MCMC/out_files"
      #outname="single_R_fast"
      
      ID <- paste( sample(x=1:9, size=9, replace=TRUE), collapse="")
      set.seed(ID)
      H_fixed=T
      sparse_sp=F
      simMiss=F
      nMiss=0
      tips=200
      #n=tips/10
      norm = T
      unif = F
      reps=1
      pred=2
      batches=1
      set_prevalence=F
      prev=.1
      grid_size<-1500
      #N<-1000
      span=1.5
      PA_plot=F
      
      
      v=c(.05, .05)
      k=lapply(1:length(v), function(x) log(v[[x]]/(1-v[[x]])))
      
      
      kappa_scale=1.0
      
      
      phylo_unscaled <- lapply(1:reps, function(x) pbtree(n=tips,scale=1) )
      
      phylo_kappa <- lapply(1:reps, function(x) rescale(  phylo_unscaled[[x]] , "kappa", kappa_scale))
      
      
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
    
    
    if(unif==T){
      
      
      par.mu.unif.H.fixed<-list( matrix(c(-.25, .25, 0.2, 0.4), nrow = 2, ncol = 2, byrow = TRUE) #pred 1 root prior
                                 , matrix(c(-.25, .25, 0.2, 0.4), nrow = 2, ncol = 2, byrow = TRUE) #pred 2 root prior
      )
      max_bd<-1.5
      min_bd<-.1
      bd<-((log(max_bd)-log(min_bd))/4)^2
      backTransform1(forwardTransform1(c(.5, bd, .06)))
      par.sd.unif.max<-list(c(.2, .2) #pred 1
                            ,c(.2, .2) #pred 2
      )
      par.sd.unif.min<-list(c(.05, .05) #pred 1
                            ,c(.05, .05) #pre 2
      )
      par.sd.unif <-lapply(1:length(par.sd.unif.max), function(x) cbind(par.sd.unif.min[[x]],par.sd.unif.max[[x]]))
      Prior_scale <-lapply(1:length(par.sd.unif.max), function(x) makePrior_ENE(r = 2, p = 1, den.mu = "unif", par.mu = par.mu.unif.H.fixed[[x]], den.sd = "unif" , par.sd = par.sd.unif[[x]]))
      
      TruePars_scale<-lapply(1:reps, function(x) priorSim_pars(Prior_scale, phylo[[x]], dist="unif"))
      startPars_scaled<-lapply(1:reps, function(x) priorSim_pars(Prior_scale, phylo[[x]], dist="unif") )
      
      
      
    }else if(norm==T){
      #norm priors
      
      
      par.mu.norm.H.fixed<-list( matrix(c(  0,  .2, #center(norm)
                                            log(0.3), .1), #width (lnorm)
                                        nrow = 2, ncol = 2, byrow = TRUE) #pred 1 root prior (center mean, center sd, width mean, width sd)
                                 , matrix(c(0, .2,
                                            log(0.3), .1),
                                          nrow = 2, ncol = 2, byrow = TRUE) #pred 2 root prior (center mean, center sd, width mean, width sd)
      )
      par.sd.lnorm.meanlog<-list(c(log(.1), log(.1)) #pred 1 (center, width, height)
                                 ,c(log(.1), log(.1)) #pred 2
      )
      par.sd.lnorm.sdlog<-list(c(.5, .5) #pred 1
                               ,c(.5, .5) #pred 2
      )
      
      par.sd.lnorm <-lapply(1:length(par.sd.lnorm.meanlog), function(x) cbind(par.sd.lnorm.meanlog[[x]], par.sd.lnorm.sdlog[[x]]))
      par.sd.lnorm <-lapply(1:length(par.sd.lnorm.meanlog), function(x) cbind(par.sd.lnorm.meanlog[[x]], par.sd.lnorm.sdlog[[x]]))
      
      #Prior_scale <-lapply(1:length(par.sd.lnorm.meanlog), function(x) makePrior_ENE(r = 2, p = 1, den.mu = "norm", heights_by_sp = sample(.95,size = tips, replace=T), par.mu = par.mu.norm.H.fixed[[x]], den.sd = "lnorm" , par.sd = par.sd.lnorm[[x]], plot=T))
      #Prior_scale <-lapply(1:length(par.sd.lnorm.meanlog), function(x) makePrior_ENE(r = 2, p = 1, den.mu = "norm", heights_by_sp = sample(.95,size = tips, replace=T), par.mu = par.mu.norm.H.fixed[[x]], den.sd = "lnorm" , par.sd = par.sd.lnorm[[x]], plot=T))
      #TruePars_scale <-lapply(1:reps, function(x) priorSim_pars(Prior_scale, hard_coded_heights = c(0.95, 0.95), phylo_kappa[[x]], dist="norm"))
      
      
      Prior_scale <-lapply(1:length(par.sd.lnorm.meanlog), function(x) makePrior_ENE(r = 2, p = 1, den.mu = "norm", heights_mean_by_sp = sample(.95,size = tips, replace=T), heights_sd_by_sp = sample(.15,size = tips, replace=T), par.mu = par.mu.norm.H.fixed[[x]], den.sd = "lnorm" , par.sd = par.sd.lnorm[[x]], plot=T))
      TruePars_scale <-lapply(1:reps, function(x) priorSim_pars(Prior_scale, hard_coded_heights = c(0.95, 0.95), phylo_kappa[[x]], dist="norm"))
      
      #startPars_scaled <-lapply(1:reps, function(x) priorSim_pars(Prior_scale, hard_coded_heights = c(0.95, 0.95), phylo[[x]], dist="norm") )
      #TruePars_scale[[1]]$sim_dat$sim_dat_bt
      job_id = as.numeric(Sys.getenv("SLURM_JOB_ID"))
      
      cat("Starting run:",job_id, ID, "\n")
      
      
      dir<-("outfiles/")
      
      
      
      dir.create(dir)
      
    }
    
    
    
    
    #presence absences data then simulated from scaled traits and beta coefficients
    
    #pa_data_scaled<-lapply(1:reps, function(x) simPA( res=TruePars_scale[[x]]$sim_dat$sim_dat_bt,
    #                                                  tree=phylo_kappa[[x]], span=span,
    #                                                  grid_size=grid_size,
    #                                                  simMiss=simMiss,
    #                                                  nMiss=nMiss,addbackground = T, A2P_ratio = 10
    #)
    #)
    
    
    pa_data_scaled<-lapply(1:reps, function(x) BePhyNE::simPA( res=TruePars_scale[[x]]$sim_dat$sim_dat_bt,
                                                      tree=phylo_kappa[[x]], span=span,
                                                      grid_size=grid_size,
                                                      simMiss=simMiss,
                                                      nMiss=nMiss
    )
    )
    
    
    #i=170
    #i=i+1
    #plot(pa_data_scaled[[1]]$species_data[[i]]$X1,  pa_data_scaled[[1]]$species_data[[i]]$X2, col=c(1,2)[pa_data_scaled[[1]]$species_data[[i]]$y+1])
    #pdf("~/Downloads/kappa_pa_test.pdf")
    ##plot(pa_data_scaled[[1]]$species_data[[i]]$X1,  pa_data_scaled[[1]]$species_data[[i]]$X2, col=c(1,2)[pa_data_scaled[[1]]$species_data[[i]]$y+1])
    #par(mfrow = c(4, 2),          # Adjust layout (rows, columns)
    #    mar = c(2.5, 2.5, 1, 1),  # Reduce space around each plot
    #    oma = c(0, 0, 0, 0),      # Outer margins
    #    mgp = c(1.5, 0.5, 0))     # Move axis labels closer
    #i=1
    ##plot(pa_data_scaled[[1]]$species_data[[i]]$X1[pa_data_scaled[[1]]$species_data[[i]]$y==1],  pa_data_scaled[[1]]$species_data[[i]]$X2[pa_data_scaled[[1]]$species_data[[i]]$y==1], col=c(sample(1:10)), ylim =c(-2,2), xlim=c(-2,2), pch=16)
    #
    #for(i in 1:200){
    #  plot(pa_data_scaled[[1]]$species_data[[i]]$X1,  pa_data_scaled[[1]]$species_data[[i]]$X2, col=c(0,2)[pa_data_scaled[[1]]$species_data[[i]]$y+1], pch=19, cex=0.3)
    #  points(pa_data_scaled[[1]]$species_data[[i]]$X1,  pa_data_scaled[[1]]$species_data[[i]]$X2, col=c(1,0)[pa_data_scaled[[1]]$species_data[[i]]$y+1], pch=19, cex=0.1)
    #  
    ##plot(pa_data_scaled[[1]]$species_data[[i]]$X1[pa_data_scaled[[1]]$species_data[[i]]$y==1],  pa_data_scaled[[1]]$species_data[[i]]$X2[pa_data_scaled[[1]]$species_data[[i]]$y==1], col=2)
    #}
    #
    #dev.off()
    
    startPars_scaled=list()
    #BePhyNE::lnL_ratematrix()
    
    #startPars_scaled=GLM_only_ml$
    
    
    #lapply(1:length(TruePars_scale), function(x) unlist(c(TruePars_scale[[x]]$A$A_bt, TruePars_scale[[x]]$R$R_sd, TruePars_scale[[x]]$R$R_cor[[1]][1,2], TruePars_scale[[x]]$R$R_cor[[2]][1,2])))
    
    
    ###set tunning parameters and Run#####
    center_slide=.18
    center_mult=.12
    width_slide=.15
    width_mult=.2
    height_slide=.5
    height_mult=.3
    
    tuning<-  list(
      niche_prop= lapply(1:pred, function(pred) list(slide=tibble(center=sample(center_slide, length(phylo_unscaled[[1]]$tip.label), replace=T),
                                                                  width= sample(width_slide,  length(phylo_unscaled[[1]]$tip.label), replace=T),
                                                                  height= sample(height_slide,length(phylo_unscaled[[1]]$tip.label), replace=T) ),
                                                     multi=tibble(center=sample(center_mult , length(phylo_unscaled[[1]]$tip.label), replace=T),
                                                                  width= sample(width_mult ,  length(phylo_unscaled[[1]]$tip.label), replace=T),
                                                                  height= sample(height_mult, length(phylo_unscaled[[1]]$tip.label), replace=T)  )
      )),
      w_mu =  lapply(1:pred, function(pred) list(slide=c(.6,.6),
                                                 multi=c(.6,.6)))
      ,
      w_sd =  lapply(1:pred, function(pred)  list(slide=c(.15,.12),
                                                  multi=c(.15,.12))
      ),
      v_cor       = lapply(1:pred, function(pred) 100)
    )
    
    true = lapply(1:reps, function(i) true=list(TruePars_scale[[i]]$sim_dat$sim_dat_bt[[1]],
                                                TruePars_scale[[i]]$sim_dat$sim_dat_bt[[2]])
    )
    
    moves_wieghts=c("height" =2,
                    "center" =3,
                    "width"  =3,
                    "theta"  =1,
                    "R_corr" =1,
                    "R_sd"   =1)
    
    move_prob=c("height"  = moves_wieghts[[1]]/sum(moves_wieghts),
                "center" = moves_wieghts[[2]]/sum(moves_wieghts),
                "width"  = moves_wieghts[[3]]/sum(moves_wieghts),
                "theta"  = moves_wieghts[[4]]/sum(moves_wieghts),
                "R_corr" = moves_wieghts[[5]]/sum(moves_wieghts),
                "R_sd"   = moves_wieghts[[6]]/sum(moves_wieghts))
    
    sparse_sp=F
    iterations=1000000
    trim_freq=iterations/10000
    chain_end=(iterations-(iterations/10))/trim_freq
    
    job_id = as.numeric(Sys.getenv("SLURM_JOB_ID"))
    #job_rep = as.numeric(Sys.getenv("PARALLEL_JOBSLOT"))
    
    ID <- paste( sample(x=1:9, size=9, replace=TRUE), collapse="")
    cat("Starting run:",job_id, ID, "\n")
    
    
    i=reps
    
    
    for (i in 1:reps){
      
      
      rep=0
      repeat{
        print("rep")
        print(rep)
        rep=rep+1
        #print
        startPars_scaled[[i]] <- priorSim_pars(Prior_scale, phylo_unscaled[[i]], dist="norm",hard_coded_heights = c(0.95, 0.95))
        #break when no tip has a -inf likelihood
        
        
        test=try(  metro_haste_full_MV(       R_corr_start = startPars_scaled[[i]]$R$R_cor
                                              , R_sd_start   = startPars_scaled[[i]]$R$R_sd
                                              , A_start      = startPars_scaled[[i]]$A$A_bt
                                              , Prior_scale
                                              , phylo_unscaled[[i]] #testing how well a kappa dataset handles under our model when we can observe the speciation driven patterns
                                              , tibble_data = startPars_scaled[[i]]$sim_dat$sim_td_bt
                                              , pa_data     = pa_data_scaled[[i]]$species_data
                                              , iterations=100
                                              , burnin=1
                                              , move_prob=move_prob
                                              , n=2
                                              , print.i.freq=10
                                              , print.ac.freq=10
                                              , printing=TRUE
                                              , trim=T
                                              , trim_freq=trim_freq
                                              , H_fixed=T
                                              , tuning=tuning
                                              , center_fixed=F
                                              , write_file=F
                                              , IDlen=5
                                              , dir=NA
                                              , outname=NA
                                              , prior_only= F
                                              , glm_only  = F
                                              , plot=F
                                              , plot_freq=iterations/5
                                              , plot_file= NA
                                              , True_pars = NA
                                              , k=k
        )
        )
        
        if(length(findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data=pa_data_scaled[[i]]$species_data, plot=F))==0 && length(test)>1){
          print("good start")
          break}
        
        
      }
      #findBadStart(res=startPars_scaled[[i]]$sim_dat$sim_dat_bt, pa_data= pa_data_scaled[[i]]$species_data, plot=F)
      
      print("done")
    }
    
    
    #for (i in 1:reps){
    print("initializing MCMC")
    
    registerDoParallel(cores=reps)
    results_set<- foreach(i=1:reps) %dopar% {
      
      #```{r,message=FALSE, warning=F, results= "hide" }
      metro_haste_full_MV(          R_corr_start = startPars_scaled[[i]]$R$R_cor
                                    , R_sd_start   = startPars_scaled[[i]]$R$R_sd
                                    , A_start      = startPars_scaled[[i]]$A$A_bt
                                    , Prior_scale
                                    , phylo_unscaled[[i]] #testing how well a kappa dataset handles under our model when we can observe the speciation driven patterns
                                    , tibble_data = startPars_scaled[[i]]$sim_dat$sim_td_bt
                                    , pa_data=pa_data_scaled[[i]]$species_data
                                    , iterations=iterations
                                    , burnin=iterations/10
                                    , move_prob=move_prob
                                    , n=2
                                    , print.i.freq=1000
                                    , print.ac.freq=100
                                    , printing=TRUE
                                    , trim=T
                                    , trim_freq=trim_freq
                                    , H_fixed=T
                                    , tuning=tuning
                                    , center_fixed=F
                                    , write_file=F
                                    , IDlen=5
                                    , dir=NA
                                    , outname=NA
                                    , prior_only= F
                                    , glm_only  = F
                                    , plot=F
                                    , plot_freq=iterations/5
                                    , plot_file= NA
                                    , True_pars = NA
                                    , k=k
      )
      #```
    }
    
    #saveRDS(results, file=paste("BePhyNE_Kappa_run",toString(sample(1:1000,1)), sep="_" ))
    
    ###check likelihoods and acceptance ratios#######
    
    
    #dir<-paste("~/R/",job_id,sep="")
    
    #dir.create(dir)
    
    #dir.create(paste(dir,"/AR",sep=""))
    
    
    MCMC_run=list()
    
    
    
    
    for (i in 1:reps){
      
      results=results_set[[i]]
      
      mymcmc_trait_full     =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="trait", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      mymcmc_trait_full_2   =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="trait", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      mymcmc_A              =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="theta", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      mymcmc_A_2            =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="theta", zeroed=F, startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      mymcmc_R_sd           =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="R_sd", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      mymcmc_R_sd_2         =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="R_sd", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      mymcmc_R_corr         =    chain2mcmcobj_full_H_fixed_MV(1, results, phylo_unscaled[[i]], object="R_corr", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      mymcmc_R_corr_2       =    chain2mcmcobj_full_H_fixed_MV(2, results, phylo_unscaled[[i]], object="R_corr", zeroed=F,  startPars_scaled[[i]]$sim_dat$sim_td, startPars_scaled[[i]]$A$A_ft, startPars_scaled[[i]]$R$R_sd, startPars_scaled[[i]]$R$R_cor)
      
      
      
      
      
      
      
      MCMC_run[[i]] <- list(chains=list(mymcmc_trait_full  =mymcmc_trait_full
                                        ,mymcmc_trait_full_2=mymcmc_trait_full_2
                                        ,mymcmc_A           =mymcmc_A
                                        ,mymcmc_A_2         =mymcmc_A_2
                                        ,mymcmc_R_sd        =mymcmc_R_sd
                                        ,mymcmc_R_sd_2      =mymcmc_R_sd_2
                                        ,mymcmc_R_corr      =mymcmc_R_corr
                                        ,mymcmc_R_corr_2    =mymcmc_R_corr_2
      ),
      TruePars_scale=TruePars_scale[[i]])
      
      
      
      
      
      #dir.create(paste(dir,"/MCMC_out",sep=""))
      
      
      
    }
    
    outputName=paste(ID,  kappa_scale,"iter", iterations, "MCMC_run",".RData",sep="_")
    
    #outputName=paste(ID, "kappa", kappa_scale,"iter", iterations, "MCMC_run",".RData",sep="_")
    outputPath=file.path(dir,outputName)
    save("MCMC_run",file=outputPath)
    
    
    
  }
  
}
