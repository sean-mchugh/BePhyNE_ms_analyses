## BePhyNE Guide: fitting data to BePhyNE
## Converted from R Markdown vignette to a normal R script.

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

run_id  <- as.integer(args[1])
#miss    <- as.integer(args[2])
int=run_id
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
pa_data <- read.csv("data/Pleth_data_vignette.csv")

class(pa_data)
#pa_data[which(pa_data$species=="Eurycea_bislineata"),]

#pa_data[which(pa_data$bio1==(min(pa_data$bio1))),]

#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   986,  -11.1))
#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   987,  -12))
#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   987,  -13))
#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   985,  -13))
#
#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   986,  -11.1))
#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   987,  -12))
#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   987,  -13))
#pa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   985,  -13))


#pa_data$PA = as.numeric(pa_data$PA)
#
#pa_data$bio12 = as.numeric(pa_data$bio12)
#pa_data$bio1 = as.numeric(pa_data$bio1)
#
##mipa_data = rbind(pa_data, c( "Eurycea_bislineata",  0,   987,  -13))
#misss=sample(1:length(tree$tip.label),1)
#miss_sp = data_final_miss[[miss]]$species

outdir = paste0("pletho_full_out_uninform_height_prior_new")

dir.create(outdir )


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
                  center=c(attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio12"],
                           attr(pres_data_scaled, "scaled:center")[colnames(pres_data_scaled)=="bio1"]))

#scale_atr = list(scale = c(unlist(scale_atr$scale)), center = unlist(scale_atr$center))

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

#data_final_new= rev(data_final_new)

#all(sort(data_final_new[[3]]$X2) == sort(data_final[[3]]$X2))
#data_final[[1]]

#sets_full <- suppressWarnings(separate.data(data_final, ratio = 0.5))
#scale_atr <- data_obj$scale

sets_full = readRDS("data/sets_full.RDS")

#saveRDS(sets_full, file = paste0(outdir,"sets_full_data_",int,".RDS") )


## -----------------------------------------------------------------------------
## Priors
## -----------------------------------------------------------------------------


bd_range <- c(0.1, 1.5)
bd<-((log(bd_range[[2]])-log(bd_range[[1]]))/4)^2

GLM_only_ml <- MLglmStartpars(species_data = data_final,tree = tree, height = NULL)
heights_glm <- lapply(1:length(GLM_only_ml$start_pars_bt), function(pred) GLM_only_ml$start_pars_bt[[pred]][,3])

#saveRDS(GLM_only_ml, file = paste0(outdir,"GLM_only_ml_",int,".RDS") )


#Prior_scale <- make_all_priors(
#  Npred,
#  length(tree$tip.label),
#  bd_range = bd_range,
#  
#  ## root priors
#  root_opt_mean      = rep(0, Npred),
#  root_opt_sd        = rep(0.5, Npred),
#  root_brdth_meanlog = rep(log(0.3), Npred),
#  root_brdth_sdlog   = rep(0.2, Npred),
#  
#  ## log-normal hyper-priors
#  sigsq_opt_meanlog   = rep(log(.2), Npred),
#  sigsq_opt_sdlog     = rep(0.1, Npred),
#  sigsq_brdth_meanlog = rep(log(bd), Npred),
#  sigsq_brdth_sdlog   = rep(0.1, Npred),
#  
#  ## heights
#  heights_by_sp = heights_glm ,
#  #sample(0.95, size = Ntips, replace = TRUE),
#  
#  ## constants forwarded to makePrior_ENE
#  r = Npred,
#  p = 1,
#  plot = TRUE
#)


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
  species_data = data_final,
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
  data = data_final,
  reps_before_POE = 1000
)

startPars$sim_dat$sim_dat_bt[[1]]
startPars$sim_dat$sim_dat_bt[[2]][,3] = startPars$sim_dat$sim_dat_bt[[1]][,3]
startPars$sim_dat$sim_td_bt[[2]][,3] = startPars$sim_dat$sim_td_bt[[1]][,3]

saveRDS(startPars, file = paste0(outdir,"startPars_",int,".RDS") )


## -----------------------------------------------------------------------------
## Tuning parameters
## -----------------------------------------------------------------------------

move_details <- make_tuning(
  tree = tree,
  pred = Npred
)


## -----------------------------------------------------------------------------
## Run MCMC
## -----------------------------------------------------------------------------


#sparse_sp <- FALSE
iterations=5000000
trim_freq= iterations/5000
burnin <- 0
chain_end <- (iterations - burnin) / trim_freq




filename <- paste0(outdir,"pars_", int)


results <- BePhyNE_MCMC(
  tree,
  pa_data = sets_full$training,
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

logdf = logdf[-(1:2500),]

log_summary <- summarize_logdf(logdf)


saveRDS(log_summary, file = paste0(outdir,"log_summary_",int,".RDS") )

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
