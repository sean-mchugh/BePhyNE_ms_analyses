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

args <- commandArgs(trailingOnly = TRUE)

#run_number <- as.integer(args[1])
addbackground  <- as.logical(args[1])

setwd("/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims")

library(MCMCpack, lib="packages/")
library(robustbase, lib="packages/")

library(evd, lib="packages/")
library(truncdist, lib="packages/")
library(MultiRNG, lib="packages/")
library(Rphylopars, lib="packages/")
library(BePhyNE, lib="packages/")

library(MCMCpack)
library(robustbase)
library(evd)
library(truncdist)
library(MultiRNG)
library(Rphylopars)
library(BePhyNE)

{
simPA<-function(res, tree, span, grid_size, simMiss=F, nMiss=1, addbackground=F,A2P_ratio=4){
  #A2P is absence to presence ratio  how many absences for every presence?
  res = TruePars_scale[[1]]$sim_dat$sim_dat_bt
  tree = phylo_kappa[[1]]
  span = span
  grid_size=grid_size
  simMiss=simMiss
  nMiss=nMiss
  #addbackground = T
  
  if(simMiss==T){
    miss<-sample(1:length(tree$tip.label), nMiss)
  }
  
  
  #generate bound for PA/ sim for each species
  max=lapply(1:length(res), function(pred) lapply( 1:nrow(res[[1]]), function(sp) round(res[[pred]][[sp,1]]+(res[[pred]][[sp,2]]), 3)))
  min=lapply(1:length(res), function(pred) lapply( 1:nrow(res[[1]]), function(sp) round(res[[1]][[sp,1]]-(res[[1]][[sp,2]]), 3)))
  
  #generate predictor vectors using seq based on what the expanded.grid output length will be (grid size being nth rooted based on # of pred)
  Pred_full<- lapply( 1:nrow(res[[1]]), function(sp) lapply(1:length(res), function(pred) seq(min[[pred]][[sp]]-span,max[[pred]][[sp]]+span, length.out = (grid_size^ (1 / length(res))) )))
  #Y_full<-lapply( 1:nrow(res[[1]]), function(sp) seq(ymin[[sp]]-span,ymax[[sp]]+span, length.out = sqrt(grid_size)))
  
  #generate predictor grid for each sp
  full_grid=lapply( 1:nrow(res[[1]]), function(sp) as.matrix(expand.grid( Pred_full[[sp]] )) )
  
  betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]])) # Convert to beta coefficients
  
  #yy_prelogit  = lapply( 1:nrow(res[[1]]), function(i) (betas[[1]][i,1] + betas[[1]][i,2]*full_grid[[i]][,1] + betas[[1]][i,3]*(full_grid[[i]][,1] ^2) + (betas[[2]][i,1] + betas[[2]][i,2]*full_grid[[i]][,2]+ betas[[2]][i,3]*(full_grid[[i]][,2] ^2) )) )
  
  #made the betas iterative so that we can flexibly simulate multiple pred
  
  yy_prelogit_sep  = lapply( 1:nrow(res[[1]]), function(sp)  lapply(1:length(res), function(pred) (betas[[pred]][sp,1] + betas[[pred]][sp,2]*full_grid[[sp]][,pred]+ betas[[pred]][sp,3]*(full_grid[[sp]][,pred]^2)  ) ) )
  
  yy_prelogit= lapply(1:length(yy_prelogit_sep), function(sp) rowSums(matrix(unlist(yy_prelogit_sep[[sp]]), ncol=length(yy_prelogit_sep[[sp]]), byrow=F)))
  
  presProb <-lapply(1:nrow(res[[1]]), function(i) 1/(1+exp(-1*yy_prelogit[[i]])) ) #convert to probability using logit link
  
  y <- lapply(1:nrow(res[[1]]), function(i) rbinom(n=length(presProb[[i]]), size=1, prob=presProb[[i]]) )
  
  names(y)<- lapply(1:nrow(res[[1]]), function(i)  "y")
  
  species_data=list()
  
  if(addbackground==T){
    
    pres=lapply(y, function(dat) (dat==1))
    
  }
  
  for (sp in  1:nrow(res[[1]])){
    
    if(addbackground==T){
      pres=lapply(y, function(dat) (dat==1))
      
      species_data[[sp]]<-list(
        species= tree$tip.label[[sp]],
        y=  c(y[[sp]], rep(0,sum(pres[[sp]])) )
      )
      
    }else{
      
      species_data[[sp]]<-list(
        species= tree$tip.label[[sp]],
        y= y[[sp]]
      )
    }    
    
    
    for(pred in 1:ncol(full_grid[[sp]])){
      
      if(addbackground==T){
        
        species_data[[sp]][[paste("X", pred, sep="") ]] <-c(full_grid[[sp]][,pred], full_grid[[sp]][pres[[sp]],pred])
        #View(species_data)
        
        species_data_df=as.data.frame(do.call(cbind, species_data[[sp]]) )
        
        
      }else{
        
        species_data[[sp]][[paste("X", pred, sep="") ]] <-full_grid[[sp]][,pred]
        
        
      }
      #species_data[[sp]][[length(species_data[[sp]])+1]] <- list((paste("X", pred, sep="") )=full_grid[[sp]][,pred])
      
      #names(species_data[[sp]][length(species_data[[sp]])])<- paste("X", pred, sep="")
    }
    
    #prune down data
    
    species_data_df=as.data.frame(do.call(cbind, species_data[[sp]]) )
    
    pres_df=species_data_df[species_data_df$y==1,]
    
    abs_df=species_data_df[species_data_df$y==0,][sample(1:sum(species_data_df$y==0), size = nrow( pres_df)*A2P_ratio, replace = F), ]
    
    new_species_data_df=rbind(pres_df, abs_df)
    
    species_data[[sp]]$species= new_species_data_df$species
    species_data[[sp]]$y= as.integer(new_species_data_df$y)
    for(i in 1:pred){
      species_data[[sp]][[paste0("X", i)]]= as.numeric(new_species_data_df[[paste0("X", i)]])

    }
  }
  
  
  if(simMiss==T){
    
    for (sp in miss){
      
      species_data[[sp]] <- list(species= tree$tip.label[[sp]],
                                 y= NA
      )
      for(pred in 1:ncol(full_grid[[sp]])){
        
        species_data[[sp]][[paste("X", pred, sep="") ]] <- NA
      }
      
    }
    return(list(species_data=species_data, missing_sp= sort(miss)))
    
  }
  
  
  
  return(list(species_data=species_data))
  
}

  
simPA = function (res, tree, span, grid_size, simMiss = F, nMiss = 1, addbackground=F) 
  {
    if (simMiss == T) {
      miss <- sample(1:length(tree$tip.label), nMiss)
    }
    max = lapply(1:length(res), function(pred) lapply(1:nrow(res[[1]]), 
                                                      function(sp) round(res[[pred]][[sp, 1]] + (res[[pred]][[sp, 
                                                                                                              2]]), 3)))
    min = lapply(1:length(res), function(pred) lapply(1:nrow(res[[1]]), 
                                                      function(sp) round(res[[1]][[sp, 1]] - (res[[1]][[sp, 
                                                                                                        2]]), 3)))
    Pred_full <- lapply(1:nrow(res[[1]]), function(sp) lapply(1:length(res), 
                                                              function(pred) seq(min[[pred]][[sp]] - span, max[[pred]][[sp]] + 
                                                                                   span, length.out = (grid_size^(1/length(res))))))
    full_grid = lapply(1:nrow(res[[1]]), function(sp) as.matrix(expand.grid(Pred_full[[sp]])))
    betas <- lapply(1:length(res), function(x) traits2coefs(res[[x]]))
    yy_prelogit_sep = lapply(1:nrow(res[[1]]), function(sp) lapply(1:length(res), 
                                                                   function(pred) (betas[[pred]][sp, 1] + betas[[pred]][sp, 
                                                                                                                        2] * full_grid[[sp]][, pred] + betas[[pred]][sp, 
                                                                                                                                                                     3] * (full_grid[[sp]][, pred]^2))))
    yy_prelogit = lapply(1:length(yy_prelogit_sep), function(sp) rowSums(matrix(unlist(yy_prelogit_sep[[sp]]), 
                                                                                ncol = length(yy_prelogit_sep[[sp]]), byrow = F)))
    presProb <- lapply(1:nrow(res[[1]]), function(i) 1/(1 + exp(-1 * 
                                                                  yy_prelogit[[i]])))
    y <- lapply(1:nrow(res[[1]]), function(i) rbinom(n = length(presProb[[i]]), 
                                                     size = 1, prob = presProb[[i]]))
    names(y) <- lapply(1:nrow(res[[1]]), function(i) "y")
    species_data = list()
    
    if(addbackground==T){
      for (sp in 1:nrow(res[[1]])) {
        
        pres = y[[sp]]==1
        
        y[[sp]] = c(y[[sp]], rep(0, sum(pres) ))
        
        full_grid[[sp]] = rbind( full_grid[[sp]], full_grid[[sp]][pres,])
      }
      
    }
    
    
    for (sp in 1:nrow(res[[1]])) {
      
      species_data[[sp]] <- list(species = tree$tip.label[[sp]], 
                                 y = y[[sp]])
      for (pred in 1:ncol(full_grid[[sp]])) {
        species_data[[sp]][[paste("X", pred, sep = "")]] <- full_grid[[sp]][, 
                                                                            pred]
      }
    }
    if (simMiss == T) {
      for (sp in miss) {
        species_data[[sp]] <- list(species = tree$tip.label[[sp]], 
                                   y = NA)
        for (pred in 1:ncol(full_grid[[sp]])) {
          species_data[[sp]][[paste("X", pred, sep = "")]] <- NA
        }
      }
      return(list(species_data = species_data, missing_sp = sort(miss)))
    }
    return(list(species_data = species_data))
  }

  
make_all_priors <- function(
    N,
    tips,
    # ── μ (root) priors ──────────────────────────────────────────
    root_opt_mean         = rep(0,        N),
    root_opt_sd           = rep(0.2,      N),
    root_brdth_meanlog    = rep(log(.3),  N),
    root_brdth_sdlog      = rep(0.1,      N),
    # ── σ (log‑normal) hyper‑priors ──────────────────────────────
    sigsq_opt_meanlog     = rep(log(.1),  N),
    sigsq_opt_sdlog       = rep(0.5,      N),
    sigsq_brdth_meanlog   = rep(log(.1),  N),
    sigsq_brdth_sdlog     = rep(0.5,      N),
    # ── heights (default hard‑coded) ─────────────────────────────
    heights_by_sp         = sample(.95, size = tips, replace = TRUE),
    # ── constants forwarded to makePrior_ENE ─────────────────────
    r = 2, p = 1,
    plot = TRUE
){
  ## 0. Coerce heights into a list of length N ------------------
  height_list <- if (is.list(heights_by_sp)){
    stopifnot(length(heights_by_sp) == N)
    heights_by_sp
  } else {
    replicate(N, heights_by_sp, simplify = FALSE)
  }
  
  ## 1. Build priors for each character -------------------------
  prior_scale <- lapply(seq_len(N), function(i){
    
    ## μ matrix (2 × 2)
    par_mu <- matrix(c(root_opt_mean[i],      root_opt_sd[i],
                       root_brdth_meanlog[i], root_brdth_sdlog[i]),
                     nrow = 2, byrow = TRUE)
    
    ## σ matrix (2 × 2)
    par_sigsq <- matrix(c(sigsq_opt_meanlog[i],  sigsq_opt_sdlog[i],
                          sigsq_brdth_meanlog[i], sigsq_brdth_sdlog[i]),
                        nrow = 2, byrow = TRUE)
    
    makePrior_ENE(
      r   = r,  p = p,
      den.mu = "norm",
      heights_by_sp = height_list[[i]],
      par.mu = par_mu,
      den.sd = "lnorm",
      par.sd = par_sigsq,
      plot   = plot
    )
  })
  
  invisible(prior_scale)
}

}


#setwd("/storage1/fs1/michael.landis/Active/BePhyNE_sims/Background_Abs_runs")

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
      
      phylo_kappa <- phylo_unscaled
      
      tree = phylo_kappa[[1]]
      
      #lapply(1:reps, function(x) geiger::rescale(  phylo_unscaled[[x]] , "kappa", kappa_scale))
      
      
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
      Prior_scale   = make_all_priors(pred,tips)
      TruePars_scale <-lapply(1:reps, function(x) priorSim_pars(Prior_scale, hard_coded_heights = rep(0.95, pred), phylo_kappa[[x]], dist="norm"))
      TruePars_scale
      #startPars_scaled <-lapply(1:reps, function(x) priorSim_pars(Prior_scale, hard_coded_heights = c(0.95, 0.95), phylo[[x]], dist="norm") )
      #TruePars_scale[[1]]$sim_dat$sim_dat_bt
      job_id = as.numeric(Sys.getenv("SLURM_JOB_ID"))
      
      cat("Starting run:",job_id, ID, "\n")
      
      
      #dir<-getwd()
      
      #dir = paste0(dir, "/Background_Abs_runs_grid")
      
      
      #dir.create(dir)
      
      dir<-("outfiles/")
      
      
    }
    
    
    
    
    #presence absences data then simulated from scaled traits and beta coefficients
    
    pa_data_scaled<-lapply(1:reps, function(x) simPA( res=TruePars_scale[[x]]$sim_dat$sim_dat_bt,
                                                      tree=phylo_kappa[[x]], span=span,
                                                      grid_size=grid_size,
                                                      simMiss=simMiss,
                                                      nMiss=nMiss,addbackground = addbackground 
    )
    )
    
    #pdf("~/Downloads/grid_back_pa_test.pdf")
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
    #  #plot(pa_data_scaled[[1]]$species_data[[i]]$X1[pa_data_scaled[[1]]$species_data[[i]]$y==1],  pa_data_scaled[[1]]$species_data[[i]]$X2[pa_data_scaled[[1]]$species_data[[i]]$y==1], col=2)
    #}
    
    ##i=170
    #pdf("~/Downloads/pa_test.pdf")
    ##plot(pa_data_scaled[[1]]$species_data[[i]]$X1,  pa_data_scaled[[1]]$species_data[[i]]$X2, col=c(1,2)[pa_data_scaled[[1]]$species_data[[i]]$y+1])
    #par(mfrow = c(4, 2),          # Adjust layout (rows, columns)
    #    mar = c(2.5, 2.5, 1, 1),  # Reduce space around each plot
    #    oma = c(0, 0, 0, 0),      # Outer margins
    #    mgp = c(1.5, 0.5, 0))     # Move axis labels closer
    #i=1
    ##plot(pa_data_scaled[[1]]$species_data[[i]]$X1[pa_data_scaled[[1]]$species_data[[i]]$y==1],  pa_data_scaled[[1]]$species_data[[i]]$X2[pa_data_scaled[[1]]$species_data[[i]]$y==1], col=c(sample(1:10)), ylim =c(-2,2), xlim=c(-2,2), pch=16)
    #for(i in 1:200){
    #  
    # plot(pa_data_scaled[[1]]$species_data[[i]]$X1[pa_data_scaled[[1]]$species_data[[i]]$y==1],  pa_data_scaled[[1]]$species_data[[i]]$X2[pa_data_scaled[[1]]$species_data[[i]]$y==1], col=c(2), ylim =c(-2,2), xlim=c(-2,2), pch=16)
    # points(pa_data_scaled[[1]]$species_data[[i]]$X1[pa_data_scaled[[1]]$species_data[[i]]$y==0],  pa_data_scaled[[1]]$species_data[[i]]$X2[pa_data_scaled[[1]]$species_data[[i]]$y==0], col=1, pch=16, cex=0.5)
    #}
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
    
    true = lapply(1:reps, function(i) true=lapply(TruePars_scale[[i]]$sim_dat$sim_dat_bt, function(x) x)
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
        startPars_scaled[[i]] <- priorSim_pars(Prior_scale, phylo_unscaled[[i]], dist="norm",hard_coded_heights = rep(0.95, pred))
        #break when no tip has a -inf likelihood
        
        
        test=try(  metro_haste_full_MV(         R_corr_start = startPars_scaled[[i]]$R$R_cor
                                              , R_sd_start   = startPars_scaled[[i]]$R$R_sd
                                              , A_start      = startPars_scaled[[i]]$A$A_bt
                                              , Prior_scale
                                              , phylo_unscaled[[i]] #testing how well a kappa dataset handles under our model when we can observe the speciation driven patterns
                                              , tibble_data = startPars_scaled[[i]]$sim_dat$sim_td_bt
                                              , pa_data     = pa_data_scaled[[i]]$species_data
                                              , iterations=100
                                              , burnin=1
                                              , move_prob=move_prob
                                              , n=3
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
      TruePars_scale=TruePars_scale[[i]]
      )
      
      
      
      
      
      #dir.create(paste(dir,"Background_Abs_runs",sep=""))
      
      
      
    }
    
    outputName=paste(ID, "background_clim_space", addbackground ,"iter", iterations, "MCMC_run",".RData",sep="_")
    
    #outputName=paste(ID, "kappa", kappa_scale,"iter", iterations, "MCMC_run",".RData",sep="_")
    outputPath=file.path(dir,outputName)
    saveRDS(MCMC_run, file=outputPath)
    
    
    
  }
  

