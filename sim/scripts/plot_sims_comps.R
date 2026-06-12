library(vioplot)


label_name=list(expression(italic(θ)),
                expression(italic(ω)),
                expression(italic(A)[italic(θ)]),
                expression(italic(A)[italic(ω)]),
                expression(italic(σ)[italic(θ)]),
                expression(italic(σ)[italic(ω)]),
                expression(italic(R)[italic(COR)])
                )


par_names = c("full_traits_c", "full_traits_w", "full_traits_A_C", "full_traits_A_W", "full_traits_Rsd_C" , "full_traits_Rsd_W", "full_traits_Rcor")

get_sim_true_vardifs = function(simvstrue, treatments, par_names){
  
  var_diff_df_list = list()

  for (i in 1:length(par_names )){
    
    var_diffs   =  lapply(1:length(treatments), function(treat) simvstrue[[treat]][[par_names[[i]]]][[1]][,1] - simvstrue[[treat]][[par_names[[i]]]][[1]][,2] )
    
    min_length  = min(unlist(lapply(var_diffs, length)))
    
    var_diff_df = do.call(cbind, lapply(var_diffs, function(x) x[sample(min_length, replace = F)]))
    
    colnames(var_diff_df) =  treatments
    
    var_diff_df_list[[i]] = var_diff_df 
    
  }
  
  names(var_diff_df_list)=par_names
  
  return(var_diff_df_list)

}

plot_vardiff_boxplots = function(var_diff_df_list, file = "~/Downloads/boxplots.pdf", label_name, cols =c(2:length(label_name))){
  

  #pdf(file)
  
  par(mfrow=(c(4,2)))
  for (par in 1:length(label_name)){
    
    #if(par<3){
    #  ylim =c(-.1,.1)
    #}else{
      ylim=c(-1*max(abs(var_diff_df_list[[par]])),max(abs(var_diff_df_list[[par]])))
    #}
    
    boxplot( var_diff_df_list[[par]], cols =cols, ylim=ylim, col=cols, main =label_name[[par]], cex.main=2)
  }
  
  #dev.off()
  
}


label_name=list(expression(italic(θ)),
                expression(italic(ω)),
                expression(italic(A)[italic(θ)]),
                expression(italic(A)[italic(ω)]),
                expression(italic(σ)[italic(θ)]),
                expression(italic(σ)[italic(ω)]),
                expression(italic(R)[italic(COR)])
)

Kappa_simvstrue = readRDS("Kappa_simmedian_vs_true_pars.RDS")

names(Kappa_simvstrue)

#Kappa_simvstrue = Kappa_simvstrue[c("full_traits_c", "full_traits_w", "full_traits_A_C", "full_traits_A_W", "full_traits_Rsd_C" , "full_traits_Rsd_W", "full_traits_Rcor")]
kappa_pars = paste0( c("0.0", "0.25", "0.5", "0.75", "1.0"))
par_names = c("full_traits_c", "full_traits_w", "full_traits_A_C", "full_traits_A_W", "full_traits_Rsd_C" , "full_traits_Rsd_W", "full_traits_Rcor")

kappa_var_diff_df_list = get_sim_true_vardifs(Kappa_simvstrue, kappa_pars , par_names)

plot_vardiff_boxplots(kappa_var_diff_df_list,"~/Manuscripts/BePhyNE/Sys_bio_sub/resub/Kappa_boxplots.pdf" , label_name)

simvstrue = Kappa_simvstrue
{
  #pdf("~/Downloads/lambda_lines_1.pdf")
  
  par(mfrow = c(4, 2),           # Plot layout
      mar = c(2.5, 2.5, 3, 1),   # Increase top margin (3 instead of 1)
      oma = c(0, 0, 0, 0),       # No outer margin
      mgp = c(1.5, 0.5, 0),      # Axis label positions
      cex.main = 1.5)            # Increase main title font size
  
  for(par in 1:length(simvstrue[[1]])){
    for(i in 1:length(simvstrue)){
      par_simsvstrue = simvstrue[[i]][[par]][[1]]
      if(i ==1){
        plot(par_simsvstrue[,2],par_simsvstrue[,1], ylab = "True", xlab= "Posterior Median",main=label_name[[par]], col = i+1, cex=0.3, pch=19)
      }else{
        points(par_simsvstrue[,2],par_simsvstrue[,1],  col = i+1, cex=0.2)
      }
      abline(lm(par_simsvstrue[,1]~par_simsvstrue[,2]), col=i+1)
      abline(a=0,b=1, cex=0.3, col= "grey")
      
    }
  }
  
  plot.new()
  par(xpd=NA)
  legend("left",col=c(2:6, "grey"), legend =  c(0, 0.25, 0.5, 0.75, 1.0, "1:1 trend"), lwd=5, cex=1.5, horiz = F, title=(paste("κ", "Transformation")), ncol = 2)
  
  #plot.new()
  #par(xpd=NA)
  #legend("left",col=c(2:5, "grey"), legend = c(0, 0.25, 0.5, 0.75, 1.0, "1:1 trend"), lwd=5, cex=2.4, horiz = F, title= (paste("κ", "Transformation")), ncol = 2)
  
  #dev.off()
  
}



#########################
lambda_simvstrue = readRDS("Lambda_simmedian_vs_true_pars.RDS")

names(Kappa_simvstrue)

#Kappa_simvstrue = Kappa_simvstrue[c("full_traits_c", "full_traits_w", "full_traits_A_C", "full_traits_A_W", "full_traits_Rsd_C" , "full_traits_Rsd_W", "full_traits_Rcor")]
lambda_pars = paste0( c("1.0", "0.0"))
par_names = c("full_traits_c", "full_traits_w", "full_traits_A_C", "full_traits_A_W", "full_traits_Rsd_C" , "full_traits_Rsd_W", "full_traits_Rcor")

lambda_var_diff_df_list = get_sim_true_vardifs(lambda_simvstrue, lambda_pars, par_names)

plot_vardiff_boxplots(lambda_var_diff_df_list,"~/Downloads/Kappa__boxplots.pdf" , label_name, c(6:7))


simvstrue = lambda_simvstrue
{
  #pdf("~/Downloads/lambda_lines_1.pdf")

  par(mfrow = c(4, 2),           # Plot layout
      mar = c(2.5, 2.5, 3, 1),   # Increase top margin (3 instead of 1)
      oma = c(0, 0, 0, 0),       # No outer margin
      mgp = c(1.5, 0.5, 0),      # Axis label positions
      cex.main = 1.5)            # Increase main title font size
  
  for(par in 1:length(simvstrue[[1]])){
    for(i in 1:length(simvstrue)){
      par_simsvstrue = simvstrue[[i]][[par]][[1]]
      if(i ==1){
        plot(par_simsvstrue[,2],par_simsvstrue[,1], ylab = "True", xlab= "Posterior Median",main=label_name[[par]], col = i+5, cex=0.3, pch=19)
      }else{
        points(par_simsvstrue[,2],par_simsvstrue[,1],  col = i+5, cex=0.2)
      }
      abline(lm(par_simsvstrue[,1]~par_simsvstrue[,2]), col=i+5)
      abline(a=0,b=1, cex=0.3, col= "grey")
      
    }
  }
  
  plot.new()
  par(xpd=NA)
  legend("left",col=c(6:7, "grey"), legend = c(1.0, 0.0, "1:1 trend"), lwd=5, cex=1.5, horiz = F, title=(paste( "λ", "Transformation")), ncol = 2)
  
  #dev.off()
  
}




############################
back_simvstrue = readRDS("grid_Background_simmedian_vs_true_pars.RDS")

background_pars = paste0( c("True Absence", "Background"))
par_names = c("full_traits_c", "full_traits_w", "full_traits_A_C", "full_traits_A_W", "full_traits_Rsd_C" , "full_traits_Rsd_W", "full_traits_Rcor")

back_var_diff_df_list = get_sim_true_vardifs(back_simvstrue, background_pars, par_names)

par(mfrow = c(4, 2),           # Plot layout
    mar = c(2.5, 2.5, 3, 1),   # Increase top margin (3 instead of 1)
    oma = c(0, 0, 0, 0),       # No outer margin
    mgp = c(1.5, 0.5, 0),      # Axis label positions
    cex.main = 1.5)            # Increase main title font size
plot_vardiff_boxplots(back_var_diff_df_list,"~/Downloads/grid_back_boxplots.pdf" , label_name,cols = (6:7))


{
  #pdf("~/Downloads/lambda_lines_1.pdf")
  simvstrue = back_simvstrue
  
  par(mfrow = c(4, 2),           # Plot layout
      mar = c(2.5, 2.5, 3, 1),   # Increase top margin (3 instead of 1)
      oma = c(0, 0, 0, 0),       # No outer margin
      mgp = c(1.5, 0.5, 0),      # Axis label positions
      cex.main = 1.5)            # Increase main title font size
  
  for(par in 1:length(simvstrue[[1]])){
    for(i in 1:length(simvstrue)){
      par_simsvstrue = simvstrue[[i]][[par]][[1]]
      if(i ==1){
        plot(par_simsvstrue[,2],par_simsvstrue[,1], ylab = "True", xlab= "Posterior Median",main=label_name[[par]], col = i+5, cex=0.3, pch=19)
      }else{
        points(par_simsvstrue[,2],par_simsvstrue[,1],  col = i+5, cex=0.2)
      }
      abline(lm(par_simsvstrue[,1]~par_simsvstrue[,2]), col=i+5)
      abline(a=0,b=1, lwd=0.5, col="grey")
      
    }
  }
  
  plot.new()
  par(xpd=NA)
  legend("left",col=c(6:7, "grey"), legend = c("True Absence", "Background", "1:1 trend"), lwd=5, cex=1.4, horiz = F, title=(paste( "Absence Data Type")), ncol = 2)
  
  #dev.off()
  
}
