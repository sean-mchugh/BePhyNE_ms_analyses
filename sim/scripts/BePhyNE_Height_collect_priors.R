{greeks <- 
  structure(list(grsym = c("α", "β", "γ", "δ", "ε", "ζ", 
                           "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", 
                           "ς", "σ", "τ", "υ", "φ", "χ", "ψ", "ω", "Α", "Β", "Γ", 
                           "Δ", "Ε", "Ζ", "Η", "Θ", "Ι", "Κ", "Λ", "Μ", "Ν", "Ξ", 
                           "Ο", "Π", "Ρ", "Σ", "Τ", "Υ", "Φ", "Χ", "Ψ", "Ω"), 
                 decUTF = c(945, 946, 947, 948, 949, 950, 951, 952, 953, 954, 
                            955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 
                            967, 968, 969, 913, 914, 915, 916, 917, 918, 919, 920, 921, 
                            922, 923, 924, 925, 926, 927, 928, 929, 931, 932, 933, 934, 
                            935, 936, 937), hexUTF = structure(c(945L, 946L, 947L, 948L, 
                                                                 949L, 950L, 951L, 952L, 953L, 954L, 955L, 956L, 957L, 958L, 
                                                                 959L, 960L, 961L, 962L, 963L, 964L, 965L, 966L, 967L, 968L, 
                                                                 969L, 913L, 914L, 915L, 916L, 917L, 918L, 919L, 920L, 921L, 
                                                                 922L, 923L, 924L, 925L, 926L, 927L, 928L, 929L, 931L, 932L, 
                                                                 933L, 934L, 935L, 936L, 937L), class = "hexmode"), htmlSym = c("&alpha;", 
                                                                                                                                "&beta;", "&gamma;", "&delta;", "&epsilon;", "&zeta;", "&eta;", 
                                                                                                                                "&theta;", "&iota;", "&kappa;", "&lambda;", "&mu;", "&nu;", 
                                                                                                                                "&xi;", "&omicron;", "&pi;", "&rho;", "&sigmaf;", "&sigma;", 
                                                                                                                                "&tau;", "&upsilon;", "&phi;", "&chi;", "&psi;", "&omega;", 
                                                                                                                                "&Alpha;", "&Beta;", "&Gamma;", "&Delta;", "&Epsilon;", "&Zeta;", 
                                                                                                                                "&Eta;", "&Theta;", "&Iota;", "&Kappa;", "&Lambda;", "&Mu;", 
                                                                                                                                "&Nu;", "&Xi;", "&Omicron;", "&Pi;", "&Rho;", "&Sigma;", 
                                                                                                                                "&Tau;", "&Upsilon;", "&Phi;", "&Chi;", "&Psi;", "&Omega;"
                                                                 ), Description = c("GREEK SMALL LETTER ALPHA", "GREEK SMALL LETTER BETA", 
                                                                                    "GREEK SMALL LETTER GAMMA", "GREEK SMALL LETTER DELTA", "GREEK SMALL LETTER EPSILON", 
                                                                                    "GREEK SMALL LETTER ZETA", "GREEK SMALL LETTER ETA", "GREEK SMALL LETTER THETA", 
                                                                                    "GREEK SMALL LETTER IOTA", "GREEK SMALL LETTER KAPPA", "GREEK SMALL LETTER LAMBDA", 
                                                                                    "GREEK SMALL LETTER MU", "GREEK SMALL LETTER NU", "GREEK SMALL LETTER XI", 
                                                                                    "GREEK SMALL LETTER OMICRON", "GREEK SMALL LETTER PI", "GREEK SMALL LETTER RHO", 
                                                                                    "GREEK SMALL LETTER FINAL SIGMA", "GREEK SMALL LETTER SIGMA", 
                                                                                    "GREEK SMALL LETTER TAU", "GREEK SMALL LETTER UPSILON", "GREEK SMALL LETTER PHI", 
                                                                                    "GREEK SMALL LETTER CHI", "GREEK SMALL LETTER PSI", "GREEK SMALL LETTER OMEGA", 
                                                                                    "GREEK CAPITAL LETTER ALPHA", "GREEK CAPITAL LETTER BETA", 
                                                                                    "GREEK CAPITAL LETTER GAMMA", "GREEK CAPITAL LETTER DELTA", 
                                                                                    "GREEK CAPITAL LETTER EPSILON", "GREEK CAPITAL LETTER ZETA", 
                                                                                    "GREEK CAPITAL LETTER ETA", "GREEK CAPITAL LETTER THETA", 
                                                                                    "GREEK CAPITAL LETTER IOTA", "GREEK CAPITAL LETTER KAPPA", 
                                                                                    "GREEK CAPITAL LETTER LAMBDA", "GREEK CAPITAL LETTER MU", 
                                                                                    "GREEK CAPITAL LETTER NU", "GREEK CAPITAL LETTER XI", "GREEK CAPITAL LETTER OMICRON", 
                                                                                    "GREEK CAPITAL LETTER PI", "GREEK CAPITAL LETTER RHO", "GREEK CAPITAL LETTER SIGMA", 
                                                                                    "GREEK CAPITAL LETTER TAU", "GREEK CAPITAL LETTER UPSILON", 
                                                                                    "GREEK CAPITAL LETTER PHI", "GREEK CAPITAL LETTER CHI", "GREEK CAPITAL LETTER PSI", 
                                                                                    "GREEK CAPITAL LETTER OMEGA")), .Names = c("grsym", "decUTF", 
                                                                                                                               "hexUTF", "htmlSym", "Description"), row.names = c(NA, -49L), class = "data.frame")

}


height_mean = 0.5

height_sd   = 1.0

args <- commandArgs(trailingOnly = TRUE)

height_mean <- as.numeric(args[1])
height_sd    <- as.numeric(args[2])

cat("height_mean:", height_mean, "\n")
cat("height_sd:", height_sd, "\n")

height_mean_label = gsub(pattern =".", replacement = "p", height_mean, fixed = T )
height_sd_label = gsub(pattern =".", replacement = "p", height_sd, fixed = T )

outdir = paste0("Neww_New_New_Height_runs/Height_runs_mean_",height_mean_label, "_sd_",height_sd_label)


library(stringr)

library(coda)
library(BePhyNE,lib="packages")#lambda_par_list=c(0, 1)
#lambda_list=c(0)
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
   
   
   mcmc_files_list= list()
   
   true_files_list = list()
   #lambda=lambda_list[[k]]
   MCMC_run=NA
   
   
   mcmc_files= list.files(paste0(outdir,"/"))
   #mcmc_files = mcmc_files[grep(paste("lambda_", lambda_list[[k]], sep=""),   mcmc_files)]
   mcmc_files_list =mcmc_files[grep("log.pars.log",   mcmc_files)]
   
   true_files_list =mcmc_files[grep("True",   mcmc_files)]
   
   numbers_mcmc = unlist(lapply(strsplit(x =  mcmc_files_list, split = ".", fixed = T), function(i) strsplit(i[[1]], "_")[[1]][[2]]))
   
   numbers_true = unlist(lapply(strsplit(x =  true_files_list, split = ".", fixed = T), function(i) strsplit(i[[1]], "_")[[1]][[2]]))
   
   all(numbers_mcmc ==numbers_true)
   
   mcmc_files_list = mcmc_files_list[numbers_mcmc %in% intersect(numbers_mcmc, numbers_true)]
   mcmc_files_list
   true_files_list  = true_files_list[numbers_true %in% intersect(numbers_mcmc, numbers_true)]
   true_files_list
   
   cbind( mcmc_files_list, true_files_list)
   
   traits_h_1=list()
   traits_h_2=list()
   
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
   
   cat("starting processing")
   print(outdir)
   k=1
   
   for (k in 1:length(mcmc_files_list)){
     
     
     print(k)
     # mcmc_files= list.files("Lambda_runs/")
     log_filename = mcmc_files_list[[k]]
     
     logdf = read_BePhyNE_log(paste0(outdir,"/", log_filename ))
     
     burnin = ceiling(nrow(logdf)/2)
     
     logdf = logdf[-(1:burnin),]
     
     #if(nrow(logdf)<8000){
     #  
     #  next
     #}
     logdf_summary = summarize_logdf(logdf)
     
     TruePars_scale = readRDS(paste0(outdir,"/", true_files_list[[k]] ))
     
     
     
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
   
   
   MCMC_summary_list[[length(MCMC_summary_list)+1]]=list( 
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
   
   
   print(paste("New_Height_runs/Height_runs_mean_",height_mean_label, "_sd_",height_sd_label, "_MCMC_summary_list_height.rds",  sep=""))
   saveRDS(MCMC_summary_list, paste("Neww_New_New_Height_runs/Height_runs_mean_",height_mean_label, "_sd_",height_sd_label, "_MCMC_summary_list_height.rds",  sep=""))
   
   
   MCMC_summary_list = readRDS(paste("New_Height_runs/Height_runs_mean_",height_mean_label, "_sd_",height_sd_label, "_MCMC_summary_list_height.rds",  sep=""))
   
   #xxxxxxxxxxxxxxxxxxxxx
   
   #MCMC_summary_list = readRDS(paste("MCMC_summary_list_height_p95.rds",  sep=""))
#MCMC_summary_list = readRDS(paste("Height_runs_mean_",height_mean_label, "_sd_",height_sd_label, "_MCMC_summary_list_height.rds",  sep=""))
#
###saveRDS(MCMC_summary_list, paste("Height_runs_mean_",height_mean_label, "_sd_",height_sd_label, "_MCMC_summary_list_height.rds",  sep=""))
#
#percent_error_w_1     =(MCMC_summary_list[[1]]$full_traits_w_1[,1]    -  MCMC_summary_list[[1]]$full_traits_w_1[,2]     ) #/MCMC_summary_list[[1]]$full_traits_w_1[,1]       
#percent_error_c_1     =(MCMC_summary_list[[1]]$full_traits_c_1[,1]    -  MCMC_summary_list[[1]]$full_traits_c_1[,2]     ) #/MCMC_summary_list[[1]]$full_traits_c_1[,1]       
#percent_error_w_2     =(MCMC_summary_list[[1]]$full_traits_w_2[,1]    -  MCMC_summary_list[[1]]$full_traits_w_2[,2]     ) #/MCMC_summary_list[[1]]$full_traits_w_2[,1]       
#percent_error_c_2     =(MCMC_summary_list[[1]]$full_traits_c_2[,1]    -  MCMC_summary_list[[1]]$full_traits_c_2[,2]     ) #/MCMC_summary_list[[1]]$full_traits_c_2[,1]       
#percent_error_A_C_1   =(MCMC_summary_list[[1]]$full_traits_A_C_1[,1]  -  MCMC_summary_list[[1]]$full_traits_A_C_1[,2]   ) #/MCMC_summary_list[[1]]$full_traits_A_C_1[,1]   
#percent_error_A_C_2   =(MCMC_summary_list[[1]]$full_traits_A_C_2[,1]  -  MCMC_summary_list[[1]]$full_traits_A_C_2[,2]   ) #/MCMC_summary_list[[1]]$full_traits_A_C_2[,1]   
#percent_error_Rsd_C_1 =(MCMC_summary_list[[1]]$full_traits_Rsd_C_1[,1]-  MCMC_summary_list[[1]]$full_traits_Rsd_C_1[,2] ) #/MCMC_summary_list[[1]]$full_traits_Rsd_C_1[,1] 
#percent_error_Rsd_C_2 =(MCMC_summary_list[[1]]$full_traits_Rsd_C_2[,1]-  MCMC_summary_list[[1]]$full_traits_Rsd_C_2[,2] ) #/MCMC_summary_list[[1]]$full_traits_Rsd_C_2[,1]
#percent_error_A_W_1   =(MCMC_summary_list[[1]]$full_traits_A_W_1[,1]  -  MCMC_summary_list[[1]]$full_traits_A_W_1[,2]   ) #/MCMC_summary_list[[1]]$full_traits_A_W_1[,1]   
#percent_error_A_W_2   =(MCMC_summary_list[[1]]$full_traits_A_W_2[,1]  -  MCMC_summary_list[[1]]$full_traits_A_W_2[,2]   ) #/MCMC_summary_list[[1]]$full_traits_A_W_2[,1]   
#percent_error_Rsd_W_1 =(MCMC_summary_list[[1]]$full_traits_Rsd_W_1[,1]-  MCMC_summary_list[[1]]$full_traits_Rsd_W_1[,2] ) #/MCMC_summary_list[[1]]$full_traits_Rsd_W_1[,1]
#percent_error_Rsd_W_2 =(MCMC_summary_list[[1]]$full_traits_Rsd_W_2[,1]-  MCMC_summary_list[[1]]$full_traits_Rsd_W_2[,2] ) #/MCMC_summary_list[[1]]$full_traits_Rsd_W_2[,1]
#percent_error_Rcor_1  =(MCMC_summary_list[[1]]$full_traits_Rcor_1[,1] -  MCMC_summary_list[[1]]$full_traits_Rcor_1[,2]  ) #/MCMC_summary_list[[1]]$full_traits_Rcor_1[,1] 
#percent_error_Rcor_2  =(MCMC_summary_list[[1]]$full_traits_Rcor_2[,1] -  MCMC_summary_list[[1]]$full_traits_Rcor_2[,2]  ) #/MCMC_summary_list[[1]]$full_traits_Rcor_2[,1] 
#
#
#
#error_list= list(
#     "Posterior Median - True w_1"     = percent_error_w_1    
#    ,"Posterior Median - True c_1"     = percent_error_c_1    
#    ,"Posterior Median - True w_2"     = percent_error_w_2    
#    ,"Posterior Median - True c_2"     = percent_error_c_2    
#    ,"Posterior Median - True A_C_1"   = percent_error_A_C_1  
#    ,"Posterior Median - True A_C_2"  = percent_error_A_C_2  
#    ,"Posterior Median - True Rsd_C_1" = percent_error_Rsd_C_1
#    ,"Posterior Median - True Rsd_C_2" = percent_error_Rsd_C_2
#    ,"Posterior Median - True A_W_1"   = percent_error_A_W_1  
#    ,"Posterior Median - True A_W_2"   = percent_error_A_W_2  
#    ,"Posterior Median - True Rsd_W_1" = percent_error_Rsd_W_1
#    ,"Posterior Median - True Rsd_W_2" = percent_error_Rsd_W_2
#    ,"Posterior Median - True Rcor_1"  = percent_error_Rcor_1 
#    ,"Posterior Median - True Rcor_2"  = percent_error_Rcor_2 
#  )
#
#
#dir="plots"
#
#
#pdf(file= paste(dir, "/", "Height_runs_mean_",height_mean_label, "_sd_",height_sd_label, "_par_errors.pdf", sep=""))
#
#for(i in 1:length(error_list)){
#
#  if(i<5){
#   n_per=200 #tips
#  }else{
#    n_per=1 #tips
#    
#  }
#  x= error_list[[i]]
#  group <- rep(seq_len(length(x) / n_per), each = n_per)
#  
#  # compute mean per dataset
#  means <- tapply(x, group, mean)
#  
#  
#  barplot(means, main = names(error_list)[[i]])
#
#}
#
#dev.off()
##kappa=0
#
##xxxxxxxxxxxxxxxxx
#
#
##post_out
##post_out=MCMC_summary_list[[k]][[1]]
#
#
##MCMC_summary_list[[k]]$full_traits_h_1==
##MCMC_summary_list[[k]]$full_traits_h_2
#
#
##plot(density(MCMC_summary_list[[k]]$full_traits_h_1[,1]))
#
##names(MCMC_summary_list[[k]])
#
##expression(italic(theta)[1])
#index=c(2,4,1,3,5,6,9,10, 7,8,11:14)
#
#names(MCMC_summary_list[[1]])[  index]
#
#
##label_name=list(expression(italic(ω)[italic(1)]),
##                expression(italic(θ)[italic(1)]),
##                expression(italic(ω)[italic(2)]),
##                expression(italic(θ)[italic(2)]),
##                expression(italic(A)[italic(θ1)]),
##                expression(italic(A)[italic(θ2)]),
##                expression(italic(σ)[italic(θ1)]),
##                expression(italic(σ)[italic(θ2)]),  
##                expression(italic(A)[italic(ω1)]),
##                expression(italic(A)[italic(ω2)]),
##                expression(italic(σ)[italic(ω1)]),
##                expression(italic(σ)[italic(ω2)]),
##                expression(italic(R)[italic(COR1)]),
##                expression(italic(R)[italic(COR2)]) )
##
##plot_name=label_name[[1]]
#
#
#
#index=c(3,6,2,5,7,8,11,12,9,10,13:16, 1,4)
#
##names(MCMC_summary_list[[1]])[  index]
#
#
##label_name=list(expression(italic(h)[italic(1)]),
##                expression(italic(ω)[italic(1)]),
##                expression(italic(θ)[italic(1)]),
##                expression(italic(h)[italic(2)]),
##                expression(italic(ω)[italic(2)]),
##                expression(italic(θ)[italic(2)]),
##                expression(italic(A)[italic(θ1)]),
##                expression(italic(A)[italic(θ2)]),
##                expression(italic(σ)[italic(θ1)]),
##                expression(italic(σ)[italic(θ2)]),  
##                expression(italic(A)[italic(ω1)]),
##                expression(italic(A)[italic(ω2)]),
##                expression(italic(σ)[italic(ω1)]),
##                expression(italic(σ)[italic(ω2)]),
##                expression(italic(R)[italic(COR1)]),
##                expression(italic(R)[italic(COR2)]) )
#
#
#label_name <- list(
#  expression(italic(h)[1]),
#  expression(italic(omega)[1]),
#  expression(italic(theta)[1]),
#  
#  expression(italic(h)[2]),
#  expression(italic(omega)[2]),
#  expression(italic(theta)[2]),
#  
#  expression(italic(A)[theta[1]]),
#  expression(italic(A)[theta[2]]),
#  expression(italic(sigma)[theta[1]]),
#  expression(italic(sigma)[theta[2]]),
#  
#  expression(italic(A)[omega[1]]),
#  expression(italic(A)[omega[2]]),
#  expression(italic(sigma)[omega[1]]),
#  expression(italic(sigma)[omega[2]]),
#  
#  expression(italic(R)[COR1]),
#  expression(italic(R)[COR2])
#)
#
##plot_name=label_name[[1]]
#
#
#
#plot_post_line=function(post_out, plot_name, add=F, trend_line_col="red", point_size=1){
#  min <- min(c(min(post_out[,1]),min(post_out[,2])))-.1
#  max <- max(c(max(post_out[,1]),max(post_out[,2])))+.1
#  true_lim=c(min, max)
#  final_lim=c(min, max)
#  
#  if(!add){
#    plot(NULL, NULL, ylab= "true", xlab="posterior median", main=plot_name, xlim=final_lim, ylim=true_lim, cex.main=2)
#  }
#  points(post_out[,1:2], pch=16, cex=point_size, col=trend_line_col)
#  
#  abline(lm(post_out[,2]~post_out[,1]), col=trend_line_col)
#  
#  abline(0,1, col="grey")
#  
#}
#
#
#{
#  pdf(file= paste(dir, "/", "Height_runs_mean_",height_mean_label, "_sd_",height_sd_label, "_niche_line_plot.pdf", sep=""))
#  
#  {
#    #par(mfrow = c(2, 2))
#    par(mfrow = c(4, 4))
#    
#    for(i in index){
#      
#      for(k in 1:length(MCMC_summary_list)){
#        
#        post_out=MCMC_summary_list[[k]][[i]]
#        plot_name=label_name[[i]]
#        add=(k>1)
#        trend_line_col=k+2
#        point_size=0.5
#        
#        min <- min(c(min(post_out[,1]),min(post_out[,2])))-.1
#        max <- max(c(max(post_out[,1]),max(post_out[,2])))+.1
#        true_lim=c(min, max)
#        final_lim=c(min, max)
#        
#        if(!add){
#          plot(NULL, NULL, ylab= "true", xlab="posterior median", main=plot_name, xlim=final_lim, ylim=true_lim, cex.main=2)
#        }
#        points(post_out[,1:2], pch=16, cex=point_size, col=trend_line_col)
#        
#        abline(lm(post_out[,2]~post_out[,1]), col=trend_line_col)
#        
#        abline(0,1, col="grey")
#        
#      }
#      
#      
#    }
#    ###### ##
#    
#    #plot.new()
#    #par(xpd=NA)
#    #legend("left",col=c(3:4, "grey"), legend = c( 0,1.0, "1:1 trend"), lwd=5, cex=2.4, horiz = F, title= (paste( "λ", "Transformation")), ncol = 2)
#  }
#  
#  dev.off()  
#} 
#
##bad=(200*81):(200*82)
#
##{
##  # pdf(file= paste(dir, "/",  "kappa_line_plot.pdf", sep=""))
##  
##  par(mfrow = c(4, 4))           
##  
##  for(i in index){
##    
##    for(k in 1:length(MCMC_summary_list)){
##      
##      post_out=MCMC_summary_list[[k]][[i]]
##      plot_name=label_name[[i]]
##      add=(k>1)
##      trend_line_col=k+1
##      point_size=0.5
##      
##      
##      min <- min(c(min(post_out[,1]),min(post_out[,2])))-.1
##      max <- max(c(max(post_out[,1]),max(post_out[,2])))+.1
##      true_lim=c(min, max)
##      final_lim=c(min, max)
##      
##      if(!add){
##        plot(NULL, NULL, ylab= "true", xlab="posterior median", main=plot_name, xlim=final_lim, ylim=true_lim, cex.main=2.5, cex.lab=1.8,cex.axis=1.4, cex.names=1.4)
##      }
##      points(post_out[,1:2], pch=16, cex=point_size, col=trend_line_col)
##      
##      abline(lm(post_out[,2]~post_out[,1]), col=trend_line_col)
##      
##      abline(0,1, col="grey")
##      
##    }
##    
##    
##  }
##  ###### ##
##  
##  #legend("bottomleft", col=c(2:5), legend = c(0, 0.25, 0.5, 0.75, 1))
##  
##  #plot.new()
##  #par(xpd=NA)
##  #legend("left",col=c(2:5, "grey"), legend = c(0, 0.25, 0.5, 0.75, 1.0, "1:1 trend"), lwd=5, cex=2.4, horiz = F, title= (paste("κ", "Transformation")), ncol = 2)
##  #"Κ"
##  #dev.off()  
##} 
##
##dev.off()  
##
###
###}
###post_diff_w=abs(MCMC_summary_list[[1]]$full_traits_w_2[,1]-MCMC_summary_list[[1]]$full_traits_w_2[,2])>.1
##
###post_diff_c=abs(MCMC_summary_list[[1]]$full_traits_c_2[,1]-MCMC_summary_list[[1]]$full_traits_c_2[,2])>.1
##
##
##post_diff=((post_diff_w+post_diff_c)>0)
##
##cbind((MCMC_summary_list[[1]]$full_traits_c_2[post_diff,1:2]),
##      (MCMC_summary_list[[1]]$full_traits_w_2[post_diff,1:2]),
##      ceiling((1:nrow(MCMC_summary_list[[1]]$full_traits_w_2))[post_diff]/200), 
##      (1:nrow(MCMC_summary_list[[1]]$full_traits_w_2))[post_diff]%%200)
##
##
##
##MCMC_summary_list[[1]]$full_traits_w_2[MCMC_summary_list[[1]]$full_traits_w_2[,1]>1,]
##
##
##full_percent_error_w_1    = do.call(rbind, lapply(1:length(percent_error_w_1), function(x)     percent_error_w_1[[x]]     ))
##full_percent_error_c_1    = do.call(rbind, lapply(1:length(percent_error_c_1), function(x)     percent_error_c_1[[x]]     ))
##full_percent_error_w_2    = do.call(rbind, lapply(1:length(percent_error_w_2), function(x)     percent_error_w_2[[x]]     ))
##full_percent_error_c_2    = do.call(rbind, lapply(1:length(percent_error_c_2), function(x)     percent_error_c_2[[x]]     ))
##full_percent_error_A_C_1  = do.call(rbind, lapply(1:length(percent_error_A_C_1), function(x)   percent_error_A_C_1[[x]]   ))
##full_percent_error_A_C_2  = do.call(rbind, lapply(1:length(percent_error_A_C_2), function(x)   percent_error_A_C_2[[x]]   ))
##full_percent_error_Rsd_C_1= do.call(rbind, lapply(1:length(percent_error_Rsd_C_1), function(x) percent_error_Rsd_C_1[[x]] ))
##full_percent_error_Rsd_C_2= do.call(rbind, lapply(1:length(percent_error_Rsd_C_2), function(x) percent_error_Rsd_C_2[[x]] ))
##full_percent_error_A_W_1  = do.call(rbind, lapply(1:length(percent_error_A_W_1), function(x)   percent_error_A_W_1[[x]]   ))
##full_percent_error_A_W_2  = do.call(rbind, lapply(1:length(percent_error_A_W_2), function(x)   percent_error_A_W_2[[x]]   ))
##full_percent_error_Rsd_W_1= do.call(rbind, lapply(1:length(percent_error_Rsd_W_1), function(x) percent_error_Rsd_W_1[[x]] ))
##full_percent_error_Rsd_W_2= do.call(rbind, lapply(1:length(percent_error_Rsd_W_2), function(x) percent_error_Rsd_W_2[[x]] ))
##full_percent_error_Rcor_1 = do.call(rbind, lapply(1:length(percent_error_Rcor_1), function(x)  percent_error_Rcor_1[[x]]  ))
##full_percent_error_Rcor_2 = do.call(rbind, lapply(1:length(percent_error_Rcor_2), function(x)  percent_error_Rcor_2[[x]]  ))
##
##
