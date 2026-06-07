## BePhyNE Guide: fitting data to BePhyNE
## Converted from R Markdown vignette to a normal R script.

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


#devtools::install_github( "sean-mchugh/BePhyNE", lib = "packages", upgrade = "never", dependencies = FALSE )

#devtools::install_github( "sean-mchugh/BePhyNE", upgrade = "never", dependencies = FALSE )


## -----------------------------------------------------------------------------
## Load and format data
## -----------------------------------------------------------------------------

tree <- ENA_Pleth_Tree
pa_data <- read.csv("data/Pleth_data_vignette.csv")




miss=sample(1:length(tree$tip.label),1)


#pa_data= pa_data[pa_data$species!="Aneides_aeneus",]

sp_col <- "species"
occ_col <- "PA"
env_preds <- c("bio12", "bio1")

Npred <- length(env_preds)
Ntips <- length(tree$tip.label)

#we take the scaling attributes from ALL presence data prior to downsampling 

pres_data_scaled<-readRDS("data/scaled_GBIF_clim_pres.RDS")

pres_data_scaled[,colnames(pres_data_scaled)=="bio1"]

attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]


scale_atr <- list(scale=c(attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio12"],
                             attr(pres_data_scaled, "scaled:scale")[colnames(pres_data_scaled)=="bio1"]) ,
                  center= c(attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio12"],
                               attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]))

#scale_atr = list(scale = c(unlist(scale_atr$scale)), center = unlist(scale_atr$center))
#
#data_obj <- format_BePhyNE_data(
#  pa_data,
#  tree,
#  sp_col,
#  occ_col,
#  env_preds,
#  scale_atr = scale_atr ,
#  normalize_data = T
#    
#    
#)
#
#data_final <- data_obj$data
#

#data_final_new= rev(data_final_new)

#all(sort(data_final_new[[3]]$X2) == sort(data_final[[3]]$X2))
#data_final[[1]]

data_obj <- format_BePhyNE_data(
  pa_data,
  tree,
  sp_col,
  occ_col,
  env_preds,
  scale_atr = scale_atr ,
  normalize_data = T
  
)

data_final <- data_obj$data


#sets_full <- suppressWarnings(separate.data(data_final, ratio = 0.5))
#scale_atr <- data_obj$scale

#saveRDS(sets_full,file = "data/sets_full.RDS")
sets_full = readRDS("data/sets_full.RDS")

data_final <-sets_full$training

## -----------------------------------------------------------------------------
## Drop species data
## -----------------------------------------------------------------------------

data_final_miss<-sets_full$training

data_final_miss[[miss]]$y<-NA
data_final_miss[[miss]]$X1<-NA
data_final_miss[[miss]]$X2<-NA


miss_sp = data_final_miss[[miss]]$species

## -----------------------------------------------------------------------------
## Priors
## -----------------------------------------------------------------------------



bd_range <- c(0.1, 1.5)
#bd <- bd_range[[2]] - bd_range[[1]]

#bd_range <- c(0.1, 1.5)
bd<-((log(bd_range[[2]])-log(bd_range[[1]]))/4)^2



#GLM_only_ml<-MLglmStartpars_general(species_data = data_final_miss,tree = tree, height = NULL)
#heights_glm<- lapply(GLM_only_ml$start_pars_bt, function(pred) pred[,3])
#
#for(pred in 1:length(heights_glm)){
#  
#  heights_glm[[pred]][is.na(heights_glm[[pred]])] = 0.5
#  
#  less_than_tol_range    = heights_glm[[pred]] < 0.05
#  greater_than_tol_range = heights_glm[[pred]] > 0.95
#  
#  heights_glm[[pred]][less_than_tol_range   ] = 0.051
#  heights_glm[[pred]][greater_than_tol_range] = 0.949
#}



#Prior_scale <- make_all_priors(
#  Npred,
#  length(tree$tip.label),
#  bd_range = bd_range,
#
#  ## root priors
#  root_opt_mean      = rep(0, Npred),
#  root_opt_sd        = rep(0.2, Npred),
#  root_brdth_meanlog = rep(log(0.3), Npred),
#  root_brdth_sdlog   = rep(0.1, Npred),
#
#  ## log-normal hyper-priors
#  sigsq_opt_meanlog   = rep(log(bd), Npred),
#  sigsq_opt_sdlog     = rep(0.5, Npred),
#  sigsq_brdth_meanlog = rep(log(bd), Npred),
#  sigsq_brdth_sdlog   = rep(0.5, Npred),
#
#  ## heights
#  heights_by_sp = heights_glm ,#sample(0.95, size = Ntips, replace = TRUE),
#
#  ## constants forwarded to makePrior_ENE
#  r = Npred,
#  p = 1,
#  plot = TRUE
#)


#N=Npred
#tips= length(tree$tip.label)
#root_opt_mean      = rep(0, Npred)
#root_opt_sd        = rep(0.5, Npred)
#root_brdth_meanlog = rep(log(0.3), Npred)
#root_brdth_sdlog   = rep(0.2, Npred)
#
### log-normal hyper-priors
#sigsq_opt_meanlog   = rep(log(.2), Npred)
#sigsq_opt_sdlog     = rep(0.1, Npred)
#sigsq_brdth_meanlog = rep(log(bd), Npred)
#sigsq_brdth_sdlog   = rep(0.1, Npred)
#
### heights
#heights_sd_by_sp = 0.15
#use_glm_height_mean = T
#species_data = data_final_miss
#tree = tree#sample(0.95, size = Ntips, replace = TRUE),
#
#
### constants forwarded to makePrior_ENE
#r = Npred
#p = 1
#plot = TRUE


Prior_scale <- make_all_priors(
  Npred,
  length(tree$tip.label),
  bd_range = bd_range,
  
  ## root priors
  root_opt_mean      = rep(0, Npred),
  root_opt_sd        = rep(0.5, Npred),
  root_brdth_meanlog = rep(log(0.3), Npred),
  root_brdth_sdlog   = rep(0.2, Npred),
  
  ## log-normal hyper-priors
  sigsq_opt_meanlog   = rep(log(.2), Npred),
  sigsq_opt_sdlog     = rep(0.1, Npred),
  sigsq_brdth_meanlog = rep(log(bd), Npred),
  sigsq_brdth_sdlog   = rep(0.1, Npred),
  
  ## heights
  heights_sd_by_sp = 1.35,
  heights_mean_by_sp = 0.53,
  use_glm_height_mean = F,
  species_data = data_final_miss,
  uninf_height_sd = 1.35,
  tree = tree,#sample(0.95, size = Ntips, replace = TRUE),
    
  
  ## constants forwarded to makePrior_ENE
  r = Npred,
  p = 1,
  plot = TRUE
)

#Prior_scale[[1]]$pars

## -----------------------------------------------------------------------------
## Starting values
## -----------------------------------------------------------------------------

startPars <- get_starting_values(
  Prior_scale = Prior_scale,
  tree = tree,
  data = data_final_miss,
  reps_before_POE = 1000
)


## -----------------------------------------------------------------------------
## Tuning parameters
## -----------------------------------------------------------------------------

move_details <- make_tuning(
  tree = tree,
  pred = Npred
)


dir.create("pletho_miss_out_height_uninform_prior_new" )

outdir = paste0("pletho_miss_out_height_uninform_prior_new/", miss_sp,"/")

dir.create(outdir )


## -----------------------------------------------------------------------------
## Run MCMC
## -----------------------------------------------------------------------------


#sparse_sp <- FALSE
iterations=5000000
trim_freq= iterations/5000
burnin <- 0
chain_end <- (iterations - burnin) / trim_freq


int=1

filename <- paste0(outdir,"pars_", int)

new_filename = filename


while(file.exists(paste0(new_filename,".pars.log"))){
  int=int+1
  new_filename = paste0(outdir,"pars_",int)
  
}

filename = new_filename 

results <- BePhyNE_MCMC(
  tree,
  pa_data = data_final_miss,
  Prior_scale,
  startPars,
  move_details,
  iterations = iterations,
  trim_freq = trim_freq,
  write2file = TRUE,
  append2existingfile = F,
  filename
)


## -----------------------------------------------------------------------------
## Read and summarize MCMC output
## -----------------------------------------------------------------------------

logdf <- read_BePhyNE_log(
  file_name = paste0(filename,".pars.log")
)

logdf = logdf[-(1:500),]
log_summary <- summarize_logdf(logdf)


saveRDS(log_summary, file = paste0(outdir,"log_summary_",int,".RDS") )

saveRDS(sets_full, file = paste0(outdir,"sets_full_data_",int,".RDS") )


## -----------------------------------------------------------------------------
## Plot posterior median response curves
## -----------------------------------------------------------------------------

#plot_summary_ridgeplot(
#  tree,
#  log_summary,
#  model_names = paste0("model", 1),
#  predictor_names = env_preds,
#  scale_atr = scale_atr,
#  curve_colors = rep(make.transparent("blue", 255 / 255), 1),
#  curve_fill_colors = rep(make.transparent("blue", 30 / 255), 1),
#  line_types = rep(1, 1),
#  xlims = NULL
#)


## -----------------------------------------------------------------------------
## Continuous stochastic character mapping
## -----------------------------------------------------------------------------

## Install these if needed:
## devtools::install_github("https://github.com/bstaggmartin/evorates")
## devtools::install_github("https://github.com/bstaggmartin/contsimmap")

#library(evorates)
#library(contsimmap)
#
#char_names <- c(
#  "bio12_optimum",
#  "bio12_breadth",
#  "bio1_optimum",
#  "bio1_breadth"
#)
#
#contsimmaps <- make_simmaps_BePhyNE(
#  tree = tree,
#  logdf = logdf,
#  nsims = 10,
#  char_names = char_names
#)
#
#plot_BePhyNE_simmap(
#  contsimmaps,
#  scale_atr = scale_atr
#)


## -----------------------------------------------------------------------------
## Model adequacy / AUC
## -----------------------------------------------------------------------------

predict_stats_list <- AUC_posterior_median(
  log_summary,
  sets_full$predicting
)

saveRDS(predict_stats_list, file = paste0(outdir,"AUC_predict_stats_list_",int,".RDS") )


#plot_AUC_treebarplot(
#  tree,
#  predict_stats_list
#)
