library(remotes)
remotes::install_github("grosed/ButterflyLS/R-package")

library(ButterflyLS)

data(dark_green_fritillary_weekly)
dgf_week <- dark_green_fritillary_weekly


dgf_week <- dgf_week %>% filter(SITENO%in%sample(unique(dgf_week$SITENO),30)) %>% filter(YEAR > 2008)

output <- analysis_multiyear(dgf_week,"weekly","slope")



###############################################
### output transform to lifespan

function(output,setup_level,phi_type){
  trans_phi <- output[[1]][["phi.out"]]
  if (setup_level=="weekly"){
    trans_phi_day <- trans_phi^(1/7)
    ls_phi <- (1/(1-trans_phi_day))
  } else {
    ls_phi <- (1/(1-trans_phi))
  }
  
  Hess = output[[1]][["modelfit"]][["Hessian"]]
  inv.Hess = solve(Hess)
  res = sqrt(diag(inv.Hess))
  
  if (phi_type=="slope"){
    

    phi.cov = scale(1:length(output[[1]][["mu1.est"]]))
    se_phi_int <- exp(output[[1]][["phi.slope"]]*phi.cov)*exp(output[[1]][["phi.int"]])
    se_phi_slope <- phi.cov*exp(output[[1]][["phi.int"]])*exp(output[[1]][["phi.slope"]]*phi.cov)
    
    
    
    se_lifespan = vector(mode="numeric",length=length(output[[1]][["mu1.est"]]))
    for (i in 1:length(output[[1]][["mu1.est"]])){
      se_phi_mat <- as.matrix(c(se_phi_int[i],se_phi_slope[i],se_phi_slope2[i]))
      
      se_phi_row_mat <- t(se_phi_mat)
      
      var_cor_phi <- inv.Hess[(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 1):(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 2),(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 1):(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 2)]
    
      
      se_lifespan[i] <- se_phi_row_mat %*% var_cor_phi %*% se_phi_mat
    }
    
    
    CI_lifespan = 1.96*(sqrt(se_lifespan))
   
    
    plot_df <- data.frame(lifespan=ls_phi,
                          ci_low=ls_phi-CI_lifespan,
                          ci_up=ls_phi+CI_lifespan,
                          year=1:length(output[[1]][["mu1.est"]]))
    
    
  } else if (phi_type=="variable"){
    res=res[(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 1):(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + length(output[[1]][["mu1.est"]]))]
    
    if (setup_level=="daily"){
      SE_ls_d <- (exp(output[[1]][["modelfit"]][["allval"]]))*res
    } else {
      SE_ls_d <- (1/(1-(1/(1+exp(-(output[[1]][["modelfit"]][["allval"]]))))^(1/7))^2)*((1/7)*(1/(1+exp(-(output[[1]][["modelfit"]][["allval"]])))^(-(6/7))))*((exp(-(output[[1]][["modelfit"]][["allval"]])))/(1+exp(-(output[[1]][["modelfit"]][["allval"]])))^2)*(res)
    }
    
    
    plot_df <- data.frame(lifespan=ls_phi,
                          ci_low=ls_phi-(SE_ls_d*1.96),
                          ci_up=ls_phi+(SE_ls_d*1.96),
                          year=1:length(output[[1]][["mu1.est"]]))
    
    
    
  } else if (ohi_type=="doubleslope"){
      phi.cov = scale(1:length(output[[1]][["mu1.est"]]))
      
      se_phi_int <- exp(output[[1]][["phi.slope"]]*phi.cov)*exp(output[[1]][["phi.int"]])*exp(output[[1]][["phi.slope2"]]*phi.cov^2)
      se_phi_slope <- phi.cov*exp(output[[1]][["phi.slope"]]*phi.cov)*exp(output[[1]][["phi.int"]])*exp(output[[1]][["phi.slope2"]]*phi.cov^2)
      se_phi_slope2 <- (2*phi.cov)*phi.cov*exp(output[[1]][["phi.slope"]]*phi.cov)*exp(output[[1]][["phi.int"]])*exp(output[[1]][["phi.slope2"]]*phi.cov^2)
      # se_phi_mat <- cbind(se_phi_int,se_phi_slope)
      
      
      se_lifespan = vector(mode="numeric",length=length(output[[1]][["mu1.est"]]))
      for (i in 1:length(output[[1]][["mu1.est"]])){
        se_phi_mat <- as.matrix(c(se_phi_int[i],se_phi_slope[i],se_phi_slope2[i]))
        
        se_phi_row_mat <- t(se_phi_mat)
        
        var_cor_phi <- inv.Hess[(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 1):(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 3),(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 1):(length(output[[1]][["mu1.est"]]) + length(output[[1]][["sigma.out"]]) + 3)]
        
        
        se_lifespan[i] <- se_phi_row_mat %*% var_cor_phi %*% se_phi_mat
      }
      
      
      CI_lifespan = 1.96*(sqrt(se_lifespan))
      
      
      plot_df <- data.frame(lifespan=ls_phi,
                            ci_low=ls_phi-CI_lifespan,
                            ci_up=ls_phi+CI_lifespan,
                            year=1:length(output[[1]][["mu1.est"]]))

  }
  
  return(plot_df)
}





