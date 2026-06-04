

#' Analysis code for fitting a GAI stopover model to a single year of data.
#'
#' @name analysis
#'
#' @description Code used for fitting a single year GAI stopover model. The phi parameter is used to derive lifespan estimates. Data can be added either in a "daily" or "weekly" format.
#'
#' @param data  dataframe containing the data to analysis
#' @param year  character value corresponding to the the year to analyse from the data
#' @param setup_level  character value "daily" or "weekly" dependent on how the data was sampled. 
#'
#' @return TODO  describe the ouput structure here 
#'
#' @examples
#'
#' \donttest{
#' library(ButterflyLS)
#'
#' # lazy load data
#' data(dark_green_fritillary_weekly)
#' # year to analyse
#' year <- "2018"
#' # weekly sampling
#' setup_level <- "weekly"
#'
#' # run analysis
#' results <- analysis(dark_green_fritillary_weekly,year,setup_level)
#' }
#'
#' @export
analysis <- function(data,year,setup_level)
{

        ukbms_sp <- data
  
        ukbms_sp <- ukbms_sp %>% filter(YEAR %in% year)
        
        
        
        # Set up the values needed by the GAI
        nS <- length(unique(ukbms_sp$SITENO))
      
        
        if (setup_level== "daily") {
          nT <- length(DATE)
        } else {
          nT <-  26
        }

        
        # create matrix of values of nS and nT (which will be my value of y)
        y <- counts <- matrix(ukbms_sp$COUNT, nrow=nS, ncol=nT, byrow=TRUE)
        
        
 


    
        ##########################################
        # Model fitting code-  derived from the Biometrics paper
        ##########################################
        # optimisation method either "SANN", "Nelder-Mead", "L-BFGS-B"
        # L-BFGS-B seems to be best for the multiyear stuff. "Nelder-Mead" for single year
        meth <- "Nelder-Mead"
        # dist.type can be "P", "NB" or "ZIP"
        dist.type <- "NB"    
        # a.type can be "N", "SO" or "S"
        a.type <- "SO"  
        # If a.type is "N" or "SO"
        B <- 1  ## no. of broods
        # do I want to convert from the daily to the weekly format. Either "yes" or "no"
        convert_weekly <- "no"
        
        mu.type <- "common"
        mu.diff.type <- "common"
        sigma.type <- "hom"
        w.type <- "common"
       
        
        # Specify number of random starts 
        nstart <- 3
        
        ## Grab the functions for fitting the model


######## duplicate the model file inplace !!!




if(length(year)==1){
  # Likelihood function
  ll_func <- function(parm,irep=1,Nguess=NULL){  
    if(a.type == "N" | a.type == "SO"){
      par.index <- 0	
      # mu can be constant or varying with a spatial covariate
      mu.int <- parm[par.index+1]; par.index <- par.index + 1  
      if(mu.type == "cov"){mu.slope <- parm[par.index+1]; par.index <- par.index + 1} 
      switch(mu.type,
             common = {mu.est1 <- rep(exp(mu.int),nS)},
             cov = {mu.est1 <- exp(mu.int + mu.slope*mu.cov)})
      
      ## This bit is about broods now 
      switch(B,
             "1" = {
               sigma.est <- rep(exp(parm[par.index+1]),2); par.index <- par.index + 1},
             "2" = {
               # mu_d can be constant or varying with a spatial covariate
               mu.diff.int <- parm[par.index+1]; par.index <- par.index + 1 
               if(mu.diff.type == "cov"){mu.diff.slope <- parm[par.index+1]; par.index <- par.index + 1} 
               switch(mu.diff.type,
                      common = {mu.diff.est <- rep(exp(mu.diff.int),nS)},
                      cov = {mu.diff.est <- exp(mu.diff.int + mu.diff.slope*mu.diff.cov)})	
               # if B = 2, sigma can be the same or different for each brood
               switch(sigma.type,
                      hom = {sigma.est <- rep(exp(parm[par.index+1]),2); par.index <- par.index + 1},
                      het = {sigma.est <- exp(parm[(par.index+1):(par.index+2)]); par.index <- par.index + 2})
               # w can be constant or varying with a spatial covariate
               # This is the weight
               w.int <- parm[par.index+1]; par.index <- par.index + 1
               if(w.type == "cov"){w.slope <- parm[par.index+1]; par.index <- par.index + 1}
               switch(w.type,
                      common = {w.est <- rep(expit(w.int),nS)},
                      cov = {w.est <- expit(w.int + w.slope*w.cov)})})
      
      switch(a.type,
             N = {         ## MIXTURE MODEL METHOD
               switch(B,
                      "1" = {
                        afunc <- matrix(dnorm(rep(1:nT,each=nS),mu.est1,sigma.est),nrow=nS,ncol=nT)},     ##FOR THE SEASONAL EFFECT
                      "2" = {
                        afunc <- matrix(rep(w.est,nT)*dnorm(rep(1:nT,each=nS),mu.est1,sigma.est[1]) + (1-rep(w.est,nT))*dnorm(rep(1:nT,each=nS),mu.est1+mu.diff.est,sigma.est[2]),nrow=nS,ncol=nT)})}, ## ED: here we have mu2 = mu.est1+mu.diff.est
             SO = {        ## STOPOVER METHOD
               phi.est <-  expit(parm[par.index+1]); par.index <- par.index + 1
               switch(B,
                      "1" = {	
                        betta.est <- matrix(c(pnorm(1,mean=mu.est1,sd=sigma.est[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1,length(2:(nT-1))),sd=sigma.est[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1,length(1:(nT-2))),sd=sigma.est[1]),1-pnorm(nT-1,mean=mu.est1,sd=sigma.est[1])),nrow=nS,ncol=nT)},
                      "2" = {
                        betta.est <- matrix(rep(w.est,each=nT)*c(pnorm(1,mean=mu.est1,sd=sigma.est[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1,length(2:(nT-1))),sd=sigma.est[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1,length(1:(nT-2))),sd=sigma.est[1]),1-pnorm(rep(nT-1,each=nS),mean=mu.est1,sd=sigma.est[1])) + (1-rep(w.est,each=nT))*c(pnorm(1,mean=mu.est1 + mu.diff.est,sd=sigma.est[2]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1+ mu.diff.est,length(2:(nT-1))),sd=sigma.est[2])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1+ mu.diff.est,length(1:(nT-2))),sd=sigma.est[2]),1-pnorm(nT-1,mean=mu.est1+ mu.diff.est,sd=sigma.est[2])) ,nrow=nS,ncol=nT)})
               
               
               afunc <- betta.est				
               for(j in 2:nT){
                 for(b in 1:(j-1)){
                   afunc[,j] <- afunc[,j] + betta.est[,b]*phi.est^length(b:(j-1))
                 }
               }})
      afunc[is.na(y)] <- NA
      
      
    } else if (a.type == "S"){         ## SPLINE METHOD
      if(dist.type == "P"){alpha.est <- parm} else {alpha.est <- parm[1:(length(parm)-1)]; par.index <- length(parm)-1}      ## POISSON DIST.
      bsbasis <- bs(1:nT,df=degf,degree=deg,intercept=TRUE)
      afunc <- exp(matrix(bsbasis%*%alpha.est,nrow=nS,ncol=nT,byrow=TRUE))
      afunc <- afunc/rowSums(afunc)
      afunc[is.na(y)] <- NA
    }
    
    if(dist.type == "NB"){                   ## NEGATIVE BINOMIAL DIST.
      r.est <- exp(parm[par.index+1])}
    if(dist.type == "ZIP"){                  ## ZERO-INFLATED POISSON DIST.
      psi.est <- expit(parm[par.index+1])}
    
    # Concentrated likelihood formulation
    switch(dist.type, 
           P = { 
             llik <- dpois(y,lambda=afunc*rep(apply(y,1,sum,na.rm=TRUE)/(apply(afunc,1,sum,na.rm=TRUE)),nT),log=TRUE)},## ED: focussing on the Poisson for now, here we essentially have lambda = a*N, where N is formulated using the concentrated likelihood. We want to bring in covariates such as weather by having lambda = a*N*p
           NB = { 
             if(irep > 1){
               llik <- dnbinom(y,mu=afunc*matrix(Nguess,nrow=nS,ncol=nT),size=r.est,log=TRUE)
             } else {
               llik <- dnbinom(y,mu=afunc*rep(apply(y,1,sum,na.rm=TRUE)/(apply(afunc,1,sum,na.rm=TRUE)),nT),size=r.est,log=TRUE)}}, 
           ZIP = {
             if(irep > 1){
               llik <- dpois(y,lambda=afunc*matrix(Nguess,nrow=nS,ncol=nT),log=FALSE)
               llik[y==0 & !is.na(y)] <- log((1-psi.est) + psi.est*llik[y==0 & !is.na(y)])
               llik[y!=0 & !is.na(y)] <- log(psi.est*llik[y!=0 & !is.na(y)])
             } else {
               llik <- dpois(y,lambda=afunc*rep(apply(y,1,sum,na.rm=TRUE)/apply(afunc,1,sum,na.rm=TRUE),nT),log=FALSE)
               llik[y==0 & !is.na(y)] <- log((1-psi.est) + psi.est*llik[y==0 & !is.na(y)])
               llik[y!=0 & !is.na(y)] <- log(psi.est*llik[y!=0 & !is.na(y)])}})
    -1*sum(llik,na.rm=TRUE) 
  }
  
  
  
  #Starting values   
  start_val_func <- function(){
    psi.st <- NULL
    if(dist.type == "ZIP"){
      psi.st <- logit(sample(seq(0.5,0.9,0.1),1))}	
    r.st <- NULL
    if(dist.type == "NB"){
      r.st <- log(sample(1:5,1))}
    
    if(a.type == "S"){
      parm <- c(psi.st,sample(seq(-2,2,.1),degf),r.st)
    } else {
      if (setup_level=="weekly" |convert_weekly=="yes"){
        samp <- sort(sample(5:15,2))
        mu1.int.st <- log(samp[1])
      } else {
        samp <- sort(sample(35:105,2))
        mu1.int.st <- log(samp[1])
      }
      mu1.slope.st <- NULL
      if(mu.type == "cov") mu1.slope.st <- 0
      
      mu.diff.int.st <- mu.diff.slope.st <- NULL       
      if(B == 2){
        mu.diff.int.st <- log(samp[2]-samp[1])
        if(exp(mu.diff.int.st)<7)mu.diff.int.st <- log(7)       
        if(mu.diff.type == "cov") mu.diff.slope.st <- 0
      }
      
      if(B == 1){
        if (setup_level=="weekly"|convert_weekly=="yes"){
          sigma.st <- log(sample(2:3,1))
        } else {
          sigma.st <-log(11)
        }
        w.int.st <- w.slope.st <- NULL
      } else {
        switch(sigma.type,
               hom = {sigma.st <- log(sample(2:3,1))},
               het = {sigma.st <- rep(sample(2:3,1),2)})
        w.int.st <- logit(sample(seq(.2,.8,.1),1))
        w.slope.st <- NULL
        if(w.type == "cov") w.slope.st <- 0
      }
      switch(a.type,
             N = {
               parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st)},
             SO = {
               phi.st <- logit(sample(seq(.3,.9,.1),1))
               parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,phi.st,r.st,psi.st)
             })
    }		
    return(parm)
  }		
  
  
  
  # Equation 4 in the paper
  der_func <- function(x,model,isite,bmat){
    switch(dist.type,
           NB = {
             sum(model$modelfit$y[isite,]/x - (model$r.out+model$modelfit$y[isite,])*model$afunc.out[isite,]/(model$r.out+x*model$afunc.out[isite,]),na.rm=TRUE)},
           ZIP = {
             sum((1-bmat[isite,])*(-model$afunc.out[isite,]*model$psi.est[1]*exp(-model$afunc.out[isite,]*x))/(1-model$psi.est[1]+model$psi.est[1]*exp(-model$afunc.out[isite,]*x)) - bmat[isite,]*model$afunc.out[isite,] + bmat[isite,]*model$modelfit$y[isite,]/x,na.rm=TRUE)})
  }	
  
  # Wrapper for fitting GAI the model for multiple starts, including the iterative approach for NB and ZIP	
  fit_it_model <- function(){
    if(!(dist.type %in% c("P","NB","ZIP")))
      stop("Distribution must be P, NB or ZIP")
    if(!(a.type %in% c("N","SO","S")))
      stop("Function for {a} must be N, SO or S")
    
    fit_k <- list(); fit_k.ll <- rep(NA,nstart)
    
    for(k in 1:nstart){	
      st <- proc.time()
      irep <- 1
      fit1 <- try(fit_model(irep = irep),silent=FALSE)
      
      # If dist.type is "NB" or "ZIP" the iterative procedure is required
      if(dist.type %in% c("NB","ZIP")){
        if(dist.type == "ZIP"){
          # A matrix b indicating where y_{i,j} > 0
          bmat <- matrix(1,nrow=nS,ncol=nT)	
          bmat[is.na(y)] <- NA
          bmat[!is.na(y) & y==0] <- 0
        }
        lld <- 1
        fit <- list()
        fit[[1]] <- fit1
        ll <- NA
        ll[1] <- fit1$ll.val
        uppvals <- lowvals <- NULL
        # Iterate until convergence (here defined when the difference in likelihoods is sufficiently small)
        while(lld > 0.01){
          irep <- irep + 1
          Nest <- rep(NA,nS)
          for(isite in 1:nS){
            if(max(y[isite,],na.rm=TRUE) == 0){
              low <- 0
            } else {
              low <- 0.1
            }
            lowvals <- c(lowvals,low)
            uppvals <- c(uppvals,2500)
            # Find each N_i numerically
            temp <- try(uniroot(der_func,lower=low,upper=2500,model=fit1,isite=isite,bmat=bmat)$root,silent=TRUE)
            dtemp <- 1
            # Different options for the interval in the 1-d root finding
            lowval <- c(0,0,0,0,.1,.1,.1,.1,.01,.01,.01,.01,.01)
            uppval <- c(1000,2500,5000,10000,1000,2500,5000,10000,1000,2500,5000,10000)
            while(class(temp) == "try-error" & dtemp < 9){
              temp <- try(uniroot(der_func,lower=lowval[dtemp],upper=uppval[dtemp],bmat=bmat,model=fit1,isite=isite)$root,silent=TRUE)
              dtemp <- dtemp + 1
              lowvals <- c(lowvals,lowval)
              uppvals <- c(uppvals,uppval)
            }
            Nest[isite] <- unlist(temp)
          }
          Nest <- as.numeric(Nest)
          vals <- fit1$modelfit$allval
          fit[[irep]] <- try(fit_model(irep=irep,Nguess=Nest,vals=vals),silent=FALSE)
          ll[irep] <- fit[[irep]]$ll.val
          fit1 <- fit[[irep]]
          lld <- abs(ll[irep]-ll[irep-1])
        }
        fit1 <- fit[[irep]]
        fit1$iterations <- list(fit)
      }
      et <- proc.time()
      fit1$time <- (et-st)[3]
      if (is.na(fit1[1])){
        fit_k[[k]] <- NA
        fit_k.ll[k] <- NA
      } else{
        fit_k[[k]] <- fit1
        fit_k.ll[k] <- fit_k[[k]]$ll.val
      }
      
    }
    output <- list(fit_k[[min(c(1:nstart)[fit_k.ll==max(fit_k.ll,na.rm=T)],na.rm=T)]],fit_k,fit_k.ll)
    
    return(output)		
  }
  
  
  # Fit the GAI model
  fit_model <- function(irep=1,Nguess=NULL,vals=NULL){
    
    if(irep==1){parm <- start_val_func()} else {parm <- vals}	
    
    if(a.type == "S") { meth <- "BFGS" } else {meth <- meth}
    
    this.fit <- optim(par=parm,
                      fn=ll_func,Nguess=Nguess,irep=irep,hessian=TRUE,method= meth,
                      control=list(trace=TRUE,maxit=45000, REPORT=10, temp=30, tmax=25,pgtol=1e-7,factr=1e-9))
    
    if(is.list(this.fit) & class(try(solve(this.fit$hessian),silent=TRUE))[1] != "try-error"){	
      # Model output
      N.out <- psi.out <- r.out <- mu1.out <- mu1.int.out <- mu1.slope.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out <- sigma.out <- alpha.out <- phi.out <- betta.out <- bsbasis  <- NULL
      
      out.index <- 0
      
      if(a.type == "N" | a.type == "SO"){ 
        alpha.out <- NULL
        switch(mu.type,
               common = {mu1.out <- rep(exp(this.fit$par[out.index+1]),nS); mu1.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1},
               cov = {mu1.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu.cov); mu1.int.out <- this.fit$par[out.index+1]; mu1.slope.out <- this.fit$par[out.index+2]; out.index <- out.index + 2})
        
        switch(B,
               "1" = {
                 sigma.out <- rep(exp(this.fit$par[out.index + 1]),2); out.index <- out.index + 1}, 
               "2" = {
                 switch(mu.diff.type,
                        common = {mu.diff.out <- rep(exp(this.fit$par[out.index + 1]),nS); mu.diff.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1},
                        cov = {mu.diff.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu.diff.cov); mu.diff.int.out <- this.fit$par[out.index+1]; mu.diff.slope.out <- this.fit$par[out.index+2]; out.index <- out.index + 2})
                 
                 switch(sigma.type,
                        hom = {sigma.out <- rep(exp(this.fit$par[out.index + 1]),2); out.index <- out.index + 1},
                        het = {sigma.out <- exp(this.fit$par[(out.index + 1):(out.index + 2)]); out.index <- out.index + 2})
                 
                 switch(w.type,
                        common = {w.out <- rep(expit(this.fit$par[out.index + 1]),nS); w.int.out <- this.fit$par[out.index+1];out.index <- out.index + 1},
                        cov = {w.out <- expit(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*w.cov); w.int.out <- this.fit$par[out.index+1]; w.slope.out <- this.fit$par[out.index+2];out.index <- out.index + 2})})
        
        if(a.type == "SO"){
          phi.out <-  expit(this.fit$par[out.index+1]); out.index <- out.index + 1
        }			
      } else if(a.type == "S"){
        if(dist.type == "P"){alpha.out <- this.fit$par} else {alpha.out <- this.fit$par[1:(length(this.fit$par)-1)]; out.index <- length(this.fit$par)-1}
        bsbasis <- bs(1:nT,df=degf,degree=deg,intercept=TRUE)
        afunc.out <- exp(matrix(bsbasis%*%alpha.out,nrow=nS,ncol=nT,byrow=TRUE))
        afunc.out <- afunc.out/rowSums(afunc.out)	
      }
      
      if(dist.type == "NB"){
        r.out <- exp(this.fit$par[out.index+1])}
      if(dist.type == "ZIP"){
        psi.out <- expit(this.fit$par[out.index+1])}
      
      switch(a.type,
             N = {
               switch(B,
                      "1" = {
                        afunc.out <- matrix(dnorm(rep(1:nT,nS),mu1.out,sigma.out),nrow=nS,ncol=nT,byrow=TRUE)},
                      "2" = {
                        afunc.out <- matrix(rep(w.out,nT)*dnorm(rep(1:nT,nS),mu1.out,sigma.out[1]) + (1-rep(w.out,nT))*dnorm(rep(1:nT,nS),mu1.out + mu.diff.out,sigma.out[2]),nrow=nS,ncol=nT,byrow=TRUE)})},
             SO = {
               switch(B,
                      "1" = {	
                        betta.out <- matrix(c(pnorm(1,mean=mu1.out,sd=sigma.out[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out,length(2:(nT-1))),sd=sigma.out[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out,length(1:(nT-2))),sd=sigma.out[1]),1-pnorm(nT-1,mean=mu1.out,sd=sigma.out[1])),nrow=nS,ncol=nT)},
                      "2" = {
                        betta.out <- matrix(rep(w.out,nT)*c(pnorm(1,mean=mu1.out,sd=sigma.out[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out,length(2:(nT-1))),sd=sigma.out[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out,length(1:(nT-2))),sd=sigma.out[1]),1-pnorm(nT-1,mean=mu1.out,sd=sigma.out[1])) + (1-rep(w.out,nT))*c(pnorm(1,mean=mu1.out + mu.diff.out,sd=sigma.out[2]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out+ mu.diff.out,length(2:(nT-1))),sd=sigma.out[2])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out+ mu.diff.out,length(1:(nT-2))),sd=sigma.out[2]),1-pnorm(nT-1,mean=mu1.out+ mu.diff.out,sd=sigma.out[2])) ,nrow=nS,ncol=nT)})
               
               afunc.out <- betta.out				
               for(j in 2:nT){
                 for(b in 1:(j-1)){
                   afunc.out[,j] <- afunc.out[,j] + betta.out[,b]*phi.out^length(b:(j-1))
                 }
               }})
      
      afunc.outNA <- afunc.out
      afunc.outNA[is.na(y)] <- NA
      N.out <- apply(y,1,sum,na.rm=TRUE)/(apply(afunc.outNA,1,sum,na.rm=TRUE))
      
      
      
      output <- list(ll.val=-this.fit$value,
                     npar=length(this.fit$par),
                     N.est=N.out,
                     w.est=w.out,
                     w.int=w.int.out,
                     w.slope=w.slope.out,
                     mu1.est=mu1.out,
                     mu1.int=mu1.int.out,
                     mu1.slope=mu1.slope.out,
                     mu.diff.est=mu.diff.out,mu.diff.int=mu.diff.int.out,
                     mu.diff.slope=mu.diff.slope.out,
                     sigma.out=sigma.out,
                     r.out=r.out,
                     psi.est=psi.out,
                     phi.out=phi.out,
                     afunc.out=afunc.out,
                     afunc.outNA=afunc.outNA,
                     betta.out=betta.out,
                     bsbasis=bsbasis,
                     alpha.out=alpha.out,
                     #                   Stopover=Stopover,
                     modelfit=list(Hessian=this.fit$hessian,
                                   starts=parm,
                                   allval=this.fit$par,
                                   nS=nS,nT=nT,
                                   convergence=this.fit$convergence,
                                   y=y))
      
      output$Fitted <- output$afunc.out*output$N.est
      if(dist.type == "ZIP") output$Fitted <- output$psi.est*output$Fitted
      output$dev <- switch(dist.type,
                           P={
                             2*(sum((output$modelfit$y*log(output$modelfit$y/output$Fitted)-(output$modelfit$y-output$Fitted))[!is.na(output$modelfit$y) & output$modelfit$y!=0])+sum(output$Fitted[output$modelfit$y== 0 & !is.na(output$modelfit$y)]))},
                           NB={
                             2*(sum((output$modelfit$y*log(output$modelfit$y/output$Fitted)-((output$r.out+output$modelfit$y)*log((output$r.out+output$modelfit$y)/(output$r.out+output$Fitted))))[output$modelfit$y>0 & !is.na(output$modelfit$y)])-sum((output$r.out*log(output$r.out/(output$r.out+output$Fitted)))[!is.na(output$modelfit$y) & output$modelfit$y==0]))},
                           ZIP={NA})
      output$D <- output$dev/(length(output$modelfit$y[!is.na(output$modelfit$y)])-output$npar)
      output
    } else {NA}
  }	
} else if (length(year>1)){
  
  # Likelihood function
  ll_func <- function(parm,irep=1,Nguess=NULL){ 
    if(a.type == "N" | a.type == "SO"){
      par.index <- 0	
      # mu can be constant or varying with a spatial covariate 
      if (mu.type=="common"&vary_mu=="yes"){
        mu.int <- vector(mode="numeric", length=length(year))
        for (i in 1:length(year)){
          mu.int[i] <- parm[par.index+1]; par.index <- par.index + 1
        }
      }else if (mu.type=="common" & vary_mu=="no"|mu.type=="cov"|mu.type=="double"){
        mu.int <- parm[par.index+1]; par.index <- par.index + 1
      } else if (mu.type=="within"){
        mu.int <- vector(mode="list", length=length(year))
        tau.int <- vector(mode="numeric", length=length(year))
        for (i in 1:length(year)){
          tau.int[i] <- parm[par.index+1]; par.index <- par.index + 1
        }
        tau.slope <- parm[par.index+1]; par.index <- par.index + 1
        for (i in 1:length(year)){
          mu.int[[i]] <- exp(tau.int[i] + tau.slope*tau.cov[[i]])
        }
      } else if (mu.type=="doublewithin"){
        mu.int <- vector(mode="list", length=length(year))
        tau.int <- vector(mode="numeric", length=length(year))
        for (i in 1:length(year)){
          tau.int[i] <- parm[par.index+1]; par.index <- par.index + 1
        }
        tau.slope <- parm[par.index+1]; par.index <- par.index + 1
        tau.slope2 <- parm[par.index+1]; par.index <- par.index + 1
        for (i in 1:length(year)){
          mu.int[[i]] <- exp(tau.int[i] + tau.slope*tau.cov[[i]]+tau.slope2*tau.cov2[[i]])
        }
      }else if (mu.type=="interaction"){
        mu.int <- vector(mode="list", length=length(year))
        tau.int <- vector(mode="numeric", length=length(year))
        for (i in 1:length(year)){
          tau.int[i] <- parm[par.index+1]; par.index <- par.index + 1
        }
        tau.slope <- parm[par.index+1]; par.index <- par.index + 1
        tau.slope2 <- parm[par.index+1]; par.index <- par.index + 1
        tau.slope3 <- parm[par.index+1]; par.index <- par.index + 1
        for (i in 1:length(year)){
          mu.int[[i]] <- exp(tau.int[i] + tau.slope*tau.cov[[i]]+tau.slope2*tau.cov2[[i]]+tau.slope3*tau.cov[[i]]*tau.cov2[[i]])
        }
      }
      
      
      
      if(mu.type == "cov"){mu.slope <- parm[par.index+1]; par.index <- par.index + 1;mu.est1 <-  vector(mode = "numeric", length = length(year))} 
      if(mu.type == "double"){mu.slope <- parm[par.index+1]; par.index <- par.index + 1;mu.slope2 <- parm[par.index+1]; par.index <- par.index + 1;mu.est1 <-  vector(mode = "numeric", length = length(year))}
      if(mu.type == "within"|mu.type=="doublewithin"|mu.type=="interaction"){mu.slope <- parm[par.index+1]; par.index <- par.index + 1;mu.est1 <-  vector(mode = "list", length = length(year))}
      if(vary_mu=="yes") {mu.est1 <- vector(mode="list",length=length(year))}  
      for (i in 1:length(year)){
        switch(mu.type,
               common = {if(vary_mu=="yes"){
                 mu.est1[[i]] <- rep(exp(mu.int[i]),nS[i])
               }else{
                 mu.est1 <- rep(exp(mu.int),nS[i])
               }
               }
               ,
               
               cov = {mu.est1[i] <- exp(mu.int + mu.slope*mu.cov[i])},
               double={mu.est1[i] <- exp(mu.int + mu.slope*mu.cov[i] + mu.slope2*mu.cov2[i])},
               within={mu.est1[[i]] <- mu.int[[i]] + mu.slope*mu.cov[i]},
               doublewithin={mu.est1[[i]] <- mu.int[[i]] + mu.slope*mu.cov[i]},
               interaction={mu.est1[[i]] <- mu.int[[i]] + mu.slope*mu.cov[i]}
        )  
      }
      
      
      ## This bit is about broods now 
      switch(B,
             "1" = {
               if (sigma_type=="variable"){
                 sigma.est <-  vector(mode = "numeric", length = length(year))
                 for (i in 1:length(year)){
                   sigma.est[i] <- exp(parm[par.index+1]); par.index <- par.index + 1
                 }
               } else{
                 sigma.est <- rep(exp(parm[par.index+1]),2); par.index <- par.index + 1
               }
             },
             "2" = {
               # mu_d can be constant or varying with a spatial covariate
               mu.diff.int <- parm[par.index+1]; par.index <- par.index + 1    
               if(mu.diff.type == "cov"){mu.diff.slope <- parm[par.index+1]; par.index <- par.index + 1} 
               switch(mu.diff.type,
                      common = {mu.diff.est <- rep(exp(mu.diff.int),nS)},
                      cov = {mu.diff.est <- exp(mu.diff.int + mu.diff.slope*mu.diff.cov)})	
               # if B = 2, sigma can be the same or different for each brood
               switch(sigma.type,
                      hom = {sigma.est <- rep(exp(parm[par.index+1]),2); par.index <- par.index + 1},
                      het = {sigma.est <- exp(parm[(par.index+1):(par.index+2)]); par.index <- par.index + 2})
               # w can be constant or varying with a spatial covariate
               # This is the weight
               w.int <- parm[par.index+1]; par.index <- par.index + 1
               if(w.type == "cov"){w.slope <- parm[par.index+1]; par.index <- par.index + 1}
               switch(w.type,
                      common = {w.est <- rep(expit(w.int),nS)},
                      cov = {w.est <- expit(w.int + w.slope*w.cov)})}
      )
      
      switch(a.type,
             N = {     ## MIXTURE MODEL METHOD
               
               afunc <-  vector(mode = "list", length = length(year))
               
               
               switch(B,
                      "1" = {
                        for (i in 1:length(year)){
                          if (sigma_type=="variable" & mu.type=="cov"){
                            afunc[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu.est1[i],sigma.est[[i]][1]),nrow=nS[i],ncol=nT)
                          } else if (sigma_type=="variable" & vary_mu=="yes" & mu.type=="common") {
                            afunc[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu.est1[[i]],sigma.est[[i]][1]),nrow=nS[i],ncol=nT)
                          }else if (sigma_type=="fixed" & mu.type=="cov") {
                            afunc[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu.est1[i],sigma.est[1]),nrow=nS[i],ncol=nT)
                          }else if (sigma_type=="fixed" & vary_mu=="yes" & mu.type=="common") {
                            afunc[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu.est1[[i]],sigma.est[1]),nrow=nS[i],ncol=nT)
                          }else if(sigma_type=="fixed" & mu.type=="within"|sigma_type=="fixed" & mu.type=="doublewithin"|sigma_type=="fixed" & mu.type=="interaction") {
                            afunc[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu.est1[[i]],sigma.est[1]),nrow=nS[i],ncol=nT)
                          }else if(sigma_type=="variable" & mu.type=="within"|sigma_type=="variable" & mu.type=="doublewithin"|sigma_type=="variable" & mu.type=="interaction") {
                            afunc[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu.est1[[i]],sigma.est[[i]][1]),nrow=nS[i],ncol=nT)
                          }
                        }
                      },     
                      "2" = {
                        afunc <- matrix(rep(w.est,nT)*dnorm(rep(1:nT,each=nS),mu.est1,sigma.est[1]) + (1-rep(w.est,nT))*dnorm(rep(1:nT,each=nS),mu.est1+mu.diff.est,sigma.est[2]),nrow=nS,ncol=nT)})
             }, 
             SO = {        ## STOPOVER METHOD
               afunc <-  vector(mode = "list", length = length(year))
               betta.est <- vector(mode="list", length = length(year))
               
               if (phi_type=="variable"){
                 phi.est <-  vector(mode = "numeric", length = length(year))
                 for (i in 1:length(year)){
                   phi.est[i] <-  expit(parm[par.index+1]); par.index <- par.index + 1
                 }
               } else if (phi_type=="fixed"){
                 phi.est <- expit(parm[par.index+1]); par.index <- par.index + 1
               } else if (phi_type=="slope"){
                 phi.est <-  vector(mode = "numeric", length = length(year))
                 phi.int <- parm[par.index+1]; par.index <- par.index + 1
                 phi.slope <- parm[par.index+1]; par.index <- par.index + 1
                 for (i in 1:length(year)){
                   phi.est[i] <- expit(phi.int + phi.slope*phi.cov[i])
                 }
               } else if (phi_type=="doubleslope"){
                 phi.est <-  vector(mode = "numeric", length = length(year))
                 phi.int <- parm[par.index+1]; par.index <- par.index + 1
                 phi.slope <- parm[par.index+1]; par.index <- par.index + 1
                 phi.slope2 <- parm[par.index+1]; par.index <- par.index + 1
                 for (i in 1:length(year)){
                   phi.est[i] <- expit(phi.int + phi.slope*phi.cov[i] + phi.slope2*phi.cov2[i])
                 }
               }
               
               
               for (i in 1:length(year)){
                 switch(B,
                        "1" = {	
                          if (sigma_type=="variable" & mu.type=="cov"|sigma_type=="variable" & mu.type=="double"){
                            betta.est[[i]] <- matrix(c(pnorm(1,mean=rep(mu.est1[i],nS[i]),sd=sigma.est[i]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu.est1[i],nS[i]*length(2:(nT-1))),
                                                             sd=sigma.est[i])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu.est1[i],nS[i]*length(1:(nT-2))),
                                                               sd=sigma.est[i]),
                                                       1-pnorm(nT-1,
                                                               mean=rep(mu.est1[i],nS[i]),
                                                               sd=sigma.est[i])),nrow=nS[i],ncol=nT)
                          }else if(sigma_type=="variable" & vary_mu=="yes" & mu.type=="common"){
                            betta.est[[i]] <- matrix(c(pnorm(1,mean=mu.est1[[i]],sd=sigma.est[[i]][1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu.est1[[i]],length(2:(nT-1))),
                                                             sd=sigma.est[[i]][1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu.est1[[i]],length(1:(nT-2))),
                                                               sd=sigma.est[[i]][1]),
                                                       1-pnorm(nT-1,
                                                               mean=mu.est1[[i]],
                                                               sd=sigma.est[[i]][1])),nrow=nS[i],ncol=nT)                         
                          }else if (sigma_type=="fixed" & mu.type=="cov"|sigma_type=="fixed" & mu.type=="double"){
                            betta.est[[i]] <- matrix(c(pnorm(1,mean=rep(mu.est1[i],nS[i]),sd=sigma.est[1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu.est1[i],nS[i]*length(2:(nT-1))),
                                                             sd=sigma.est[1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu.est1[i],nS[i]*length(1:(nT-2))),
                                                               sd=sigma.est[1]),
                                                       1-pnorm(nT-1,
                                                               mean=rep(mu.est1[i],nS[i]),
                                                               sd=sigma.est[1])),nrow=nS[i],ncol=nT)
                          } else if(sigma_type=="fixed"&vary_mu=="yes"&mu.type=="common"){
                            betta.est[[i]] <- matrix(c(pnorm(1,mean=mu.est1[[i]],sd=sigma.est[1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu.est1[[i]],length(2:(nT-1))),
                                                             sd=sigma.est[1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu.est1[[i]],length(1:(nT-2))),
                                                               sd=sigma.est[1]),
                                                       1-pnorm(nT-1,
                                                               mean=mu.est1[[i]],
                                                               sd=sigma.est[1])),nrow=nS[i],ncol=nT)
                          }else if(sigma_type=="fixed"&vary_mu=="no"&mu.type=="common"){
                            betta.est[[i]] <- matrix(c(pnorm(1,mean=mu.est1,sd=sigma.est[1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu.est1,length(2:(nT-1))),
                                                             sd=sigma.est[1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu.est1,length(1:(nT-2))),
                                                               sd=sigma.est[1]),
                                                       1-pnorm(nT-1,
                                                               mean=mu.est1,
                                                               sd=sigma.est[1])),nrow=nS[i],ncol=nT)
                          }
                          
                        },
                        "2" = {
                          betta.est <- matrix(rep(w.est,each=nT)*c(pnorm(1,mean=mu.est1,sd=sigma.est[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1,length(2:(nT-1))),sd=sigma.est[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1,length(1:(nT-2))),sd=sigma.est[1]),1-pnorm(rep(nT-1,each=nS),mean=mu.est1,sd=sigma.est[1])) + (1-rep(w.est,each=nT))*c(pnorm(1,mean=mu.est1 + mu.diff.est,sd=sigma.est[2]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1+ mu.diff.est,length(2:(nT-1))),sd=sigma.est[2])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1+ mu.diff.est,length(1:(nT-2))),sd=sigma.est[2]),1-pnorm(nT-1,mean=mu.est1+ mu.diff.est,sd=sigma.est[2])) ,nrow=nS,ncol=nT)})
                 
                 
                 afunc[[i]] <- betta.est[[i]]	
                 
                 if (phi_type=="variable"|phi_type=="slope"|phi_type=="doubleslope"){
                   for(j in 2:nT){
                     for(b in 1:(j-1)){
                       afunc[[i]][,j] <- afunc[[i]][,j] + betta.est[[i]][,b]*phi.est[i]^length(b:(j-1))
                     }
                   }
                 }else {
                   for(j in 2:nT){
                     for(b in 1:(j-1)){
                       afunc[[i]][,j] <- afunc[[i]][,j] + betta.est[[i]][,b]*phi.est^length(b:(j-1))
                     }
                   }
                 }
                 
                 afunc[[i]][is.na(y[[i]])] <- NA
               }
               
             }
      )
      
      
      
    } else if (a.type == "S"){         ## SPLINE METHOD
      if(dist.type == "P"){alpha.est <- parm} else {alpha.est <- parm[1:(length(parm)-1)]; par.index <- length(parm)-1}      ## POISSON DIST.
      bsbasis <- bs(1:nT,df=degf,degree=deg,intercept=TRUE)
      afunc <- exp(matrix(bsbasis%*%alpha.est,nrow=nS,ncol=nT,byrow=TRUE))
      afunc <- afunc/rowSums(afunc)
      afunc[is.na(y)] <- NA
    }
    
    if(dist.type == "NB"){ ## NEGATIVE BINOMIAL DIST.
      if (vary_r=="no"){
        r.est <- exp(parm[par.index+1])
      }else{
        r.est <-  vector(mode = "numeric", length = length(year))
        for (i in 1:length(year)){
          r.est[i] <-  exp(parm[par.index+1]); par.index <- par.index + 1
        }
      }
    }
    if(dist.type == "ZIP"){                  ## ZERO-INFLATED POISSON DIST.
      psi.est <- expit(parm[par.index+1])}
    
    # Concentrated likelihood formulation
    llik <- vector(mode="list", length=length(year))
    for (i in 1:length(year)){
      switch(dist.type, 
             P = { 
               llik[[i]] <- dpois(y[[i]],lambda=afunc[[i]]*rep(apply(y[[i]],1,sum,na.rm=TRUE)/(apply(afunc[[i]],1,sum,na.rm=TRUE)),nT),log=TRUE)},
             
             NB = { 
               if (vary_r=="no"){
                 llik[[i]] <- dnbinom(y[[i]],mu=afunc[[i]]*rep(apply(y[[i]],1,sum,na.rm=TRUE)/apply(afunc[[i]],1,sum,na.rm=TRUE),nT),size=r.est,log=TRUE)#},
               } else {
                 llik[[i]] <- dnbinom(y[[i]],mu=afunc[[i]]*rep(apply(y[[i]],1,sum,na.rm=TRUE)/apply(afunc[[i]],1,sum,na.rm=TRUE),nT),size=r.est[i],log=TRUE)
               }
             },
             ZIP = {
               if(irep > 1){
                 llik <- dpois(y,lambda=afunc*matrix(Nguess,nrow=nS,ncol=nT),log=FALSE)
                 llik[y==0 & !is.na(y)] <- log((1-psi.est) + psi.est*llik[y==0 & !is.na(y)])
                 llik[y!=0 & !is.na(y)] <- log(psi.est*llik[y!=0 & !is.na(y)])
               } else {
                 llik <- dpois(y,lambda=afunc*rep(apply(y,1,sum,na.rm=TRUE)/apply(afunc,1,sum,na.rm=TRUE),nT),log=FALSE)
                 llik[y==0 & !is.na(y)] <- log((1-psi.est) + psi.est*llik[y==0 & !is.na(y)])
                 llik[y!=0 & !is.na(y)] <- log(psi.est*llik[y!=0 & !is.na(y)])}})
    }
    -1*sum(sapply(llik, sum, na.rm=TRUE))
  }
  
  
  
  #Starting values   
  start_val_func <- function(){
    psi.st <- NULL
    r.st <- NULL
    if(dist.type == "ZIP"){
      psi.st <- logit(sample(seq(0.5,0.9,0.1),1))	
      r.st <- NULL
    }
    if(dist.type == "NB"){
      if (vary_r=="no"){
        r.st <- log(sample(1:3,1))
      } else {
        r.st <-  vector(mode = "numeric", length = length(year))
        for (i in 1:length(year)){
          r.st[i] <- log(sample(1:5,1))
        }
      }
    }
    
    if(a.type == "S"){
      parm <- c(psi.st,sample(seq(-2,2,.1),degf),r.st)
    } else {
      if (setup_level=="weekly"|convert_weekly=="yes"){
        if (mu.type=="common"&vary_mu=="no"|mu.type=="cov"|mu.type=="double"){
          samp <- sort(sample(5:15,2))
          mu1.int.st <- log(samp[1])
        }else if (mu.type=="common"& vary_mu=="yes") {
          mu1.int.st <-  vector(mode = "numeric", length = length(year))
          for (i in 1:length(year)){
            samp <- sort(sample(5:15,2))
            mu1.int.st[i] <- log(samp[1]) 
          }
        } else if (mu.type=="within"){
          samp <- sort(sample(5:15,2))
          tau.int.st <- rep(log(samp[1]), times=length(year))
          
          
          tau.slope.st <- 0.01
        }else if (mu.type=="doublewithin"){
          
          samp <- sort(sample(5:10,2))
          tau.int.st <- rep(log(samp[1]), times=length(year)) 
          
          
          tau.slope.st <- 0.01
          tau.slope.st2 <- 0.01
        }else if (mu.type=="interaction"){
          samp <- sort(sample(5:15,2))
          tau.int.st <- rep(log(samp[1]),times=length(year))
          
          tau.slope.st <- 0.001
          tau.slope.st2 <- 0.001
          tau.slope.st3 <- 0.001
        }
      } else {
        if (vary_mu=="no"){
          samp <- sort(sample(30:100,2))
          mu1.int.st <- log(samp[1])
        }else{
          mu1.int.st <-  vector(mode = "numeric", length = length(year))
          for (i in 1:length(year)){
            samp <- sort(sample(30:100,2))
            mu1.int.st[i] <- log(samp[1]) 
          }
        }
      }
      mu1.slope.st <- NULL
      mu1.slope2.st <- NULL
      if(mu.type == "cov"|mu.type=="within"|mu.type=="doublewithin"|mu.type=="interaction") {mu1.slope.st <- -0.01} 
      if(mu.type=="double"){
        mu1.slope.st <- sample(seq(-0.1,0.1,length.out=10),1)
        mu1.slope2.st <- sample(seq(-0.1,0.1,length.out=10),1)
      }
      mu.diff.int.st <- mu.diff.slope.st <- NULL       
      if(B == 2){
        mu.diff.int.st <- log(samp[2]-samp[1])
        if(exp(mu.diff.int.st)<7)mu.diff.int.st <- log(7)      
        if(mu.diff.type == "cov") mu.diff.slope.st <- 0
      }
      
      if(B == 1){
        if (setup_level=="weekly"|convert_weekly=="yes"){
          if (sigma_type=="variable"){
            sigma.st <- rep(log(sample(2:4,1)), times=length(year))
          }else{
            sigma.st <- log(sample(2:4,1))
          }
        } else {
          if (sigma_type=="variable"){
            sigma.st <-  vector(mode = "numeric", length = length(year))
            for (i in 1:length(year)){
              sigma.st[i] <- log(sample(5:20,1))
            }
          }else{
            sigma.st <- log(sample(5:20,1))
          }
        }
        w.int.st <- w.slope.st <- NULL
      } else {
        switch(sigma.type,
               hom = {sigma.st <- log(sample(2:3,1))},
               het = {sigma.st <- rep(sample(2:3,1),2)})
        w.int.st <- logit(sample(seq(.2,.8,.1),1))
        w.slope.st <- NULL
        if(w.type == "cov") w.slope.st <- 0
      }
      switch(a.type,
             N = {
               if (sigma_type=="variable" & mu.type=="cov"){
                 parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st)
               } else if(sigma_type=="variable" & vary_mu=="yes"&mu.type=="common"){
                 parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st)
               } else if(sigma_type=="fixed" & mu.type=="cov"){
                 parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st)
               }else if(sigma_type=="fixed" & vary_mu=="yes" & mu.type=="common"){
                 parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st)
               }else if(sigma_type=="fixed" & mu.type=="within"){
                 parm <- c(tau.int.st[1:length(year)], tau.slope.st, mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st)
               } else if (sigma_type=="variable" & mu.type=="within"){
                 parm <- c(tau.int.st[1:length(year)],tau.slope.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st)
               }else if(sigma_type=="fixed" & mu.type=="doublewithin"){
                 parm <- c(tau.int.st[1:length(year)], tau.slope.st, tau.slope.st2,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st)
               } else if (sigma_type=="variable" & mu.type=="doublewithin"){
                 parm <- c(tau.int.st[1:length(year)],tau.slope.st, tau.slope.st2,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st)
               }else if (sigma_type=="variable" & mu.type=="interaction"){
                 parm <- c(tau.int.st[1:length(year)],tau.slope.st, tau.slope.st2,tau.slope.st3,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st)
               }else if(sigma_type=="fixed" & mu.type=="interaction"){
                 parm <- c(tau.int.st[1:length(year)], tau.slope.st, tau.slope.st2, tau.slope.st3,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st)
               }
               
             },
             SO = {
               if (phi_type=="variable"){
                 phi.st <-  vector(mode = "numeric", length = length(year))
                 for (i in 1:length(year)){
                   phi.st[i] <- sample(seq(.3,.9,.1),1)#logit(sample(seq(.3,.9,.1),1))
                 }
                 if (sigma_type=="variable"&mu.type=="cov"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,phi.st[1:length(year)],r.st,psi.st)
                 }else if (sigma_type=="variable"& vary_mu=="yes"& vary_r=="no" & mu.type=="common"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,phi.st[1:length(year)],r.st,psi.st)
                 }else if (sigma_type=="variable"& vary_mu=="yes" & vary_r=="yes"&mu.type=="common"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,phi.st[1:length(year)],r.st[1:length(year)],psi.st)
                 }else if(sigma_type=="fixed"& mu.type=="cov"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,phi.st[1:length(year)],r.st,psi.st)
                 }else if (sigma_type=="fixed" & vary_mu=="yes" & mu.type=="common"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,phi.st[1:length(year)],r.st,psi.st)
                 }else if (sigma_type=="fixed" & mu.type=="double"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu1.slope2.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,phi.st[1:length(year)],r.st,psi.st)
                 }else if (sigma_type=="variable" & mu.type=="double"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu1.slope2.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,phi.st[1:length(year)],r.st,psi.st)
                 }else{
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,phi.st[1:length(year)],r.st,psi.st)
                 }
                 
               } else if (phi_type=="fixed"){
                 phi.st <- logit(sample(seq(.3,.9,.1),1))
                 if (sigma_type=="variable"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,phi.st,r.st,psi.st)
                 } else{
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,phi.st,r.st,psi.st)
                 }
               } else if (phi_type=="slope"){
                 if(setup_level=="weekly"){
                   phi.int.st <- sample(seq(.3,.9,.1),1)
                   phi.slope.st <- sample(seq(-0.3,0.2,length.out=10),1)
                 }else {
                   phi.int.st <- logit(sample(seq(.8,.9,.85),1))
                   phi.slope.st <- sample(seq(-0.1,0.1,length.out=10),1)
                 }
                 if (sigma_type=="variable"&mu.type=="cov"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st)
                 }else if (sigma_type=="variable"&mu.type=="double"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu1.slope2.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st)
                 }else if (sigma_type=="variable"& vary_mu=="yes"& vary_r=="no" & mu.type=="common"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st)
                 }else if (sigma_type=="variable"& vary_mu=="yes" & vary_r=="yes" & mu.type=="common"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st[1:length(year)],psi.st,phi.int.st,phi.slope.st)
                 }else if(sigma_type=="fixed"& mu.type=="cov"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st)
                 }else if (sigma_type=="fixed" & vary_mu=="yes" & mu.type=="common"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st)
                 }else if (sigma_type=="fixed" & mu.type=="double"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu1.slope2.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st)
                 }else if (sigma_type=="fixed" & vary_mu=="no" & mu.type=="common"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st)
                 }
               } else if (phi_type=="doubleslope"){
                 if(setup_level=="weekly"){
                   phi.int.st <- logit(sample(seq(.3,.9,.1),1))
                   phi.slope.st <- sample(seq(-0.3,0.2,length.out=10),1)
                   phi.slope.st2 <- sample(seq(-0.3,0.2,length.out=10),1)
                 }else {
                   phi.int.st <- logit(sample(seq(.8,.9,.85),1))
                   phi.slope.st <- sample(seq(-0.1,0.1,length.out=10),1)
                   phi.slope.st2 <- sample(seq(-0.1,0.1,length.out=10),1)
                 }
                 if (sigma_type=="variable"&mu.type=="cov"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st, phi.slope.st2)
                 }else if (sigma_type=="variable"&vary_mu=="yes"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st[1:length(year)],w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st, phi.slope.st2)
                 }else if(sigma_type=="fixed"&mu.type=="cov"){
                   parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st,phi.slope.st2)
                 }else if (sigma_type=="fixed" & vary_mu=="yes"){
                   parm <- c(mu1.int.st[1:length(year)],mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st,phi.int.st,phi.slope.st,phi.slope.st2)
                 }
               }
             })
    }		
    return(parm)
  }		
  
  # Equation 4 in the paper
  der_func <- function(N, # parameter we are trying to solve for
                       model, # gai output from previous iteration
                       isite,
                       iyear, 
                       bmat # only relevant for ZIP
  ){
    switch(dist.type,
           NB = {
             sum(model$modelfit$y[[iyear]][isite,]/N - (model$r.out+model$modelfit$y[[i]][isite,])*model$afunc.out[[i]][isite,]/(model$r.out+N*model$afunc.out[[i]][isite,]),na.rm=TRUE)},
           ZIP = {
             sum((1-bmat[isite,])*(-model$afunc.out[isite,]*model$psi.est[1]*exp(-model$afunc.out[isite,]*x))/(1-model$psi.est[1]+model$psi.est[1]*exp(-model$afunc.out[isite,]*x)) - bmat[isite,]*model$afunc.out[isite,] + bmat[isite,]*model$modelfit$y[isite,]/x,na.rm=TRUE)})
  }	
  
  # Wrapper for fitting GAI the model for multiple starts, including the iterative approach for NB and ZIP	
  fit_it_model <- function(){
    if(!(dist.type %in% c("P","NB","ZIP")))
      stop("Distribution must be P, NB or ZIP")
    if(!(a.type %in% c("N","SO","S")))
      stop("Function for {a} must be N, SO or S")
    
    fit_k <- list(); fit_k.ll <- rep(NA,nstart)
    
    for(k in 1:nstart){	
      st <- proc.time()
      irep <- 1
      fit1 <- try(fit_model(irep = irep),silent=FALSE)
      while(is.list(fit1)==FALSE){fit1 <- try(fit_model(irep = irep),silent=FALSE)}
      # If dist.type is "NB" or "ZIP" the iterative procedure is required
      if(dist.type %in% c("NB","ZIP")){
        if(dist.type == "ZIP"){
          # A matrix b indicating where y_{i,j} > 0
          bmat <- matrix(1,nrow=nS,ncol=nT)	
          bmat[is.na(y)] <- NA
          bmat[!is.na(y) & y==0] <- 0
        }
        lld <- 1 
        fit <- list() 
        fit[[1]] <- fit1 
        ll <- NA
        ll[1] <- fit1$ll.val 
        uppvals <- lowvals <- NULL 
        # Iterate until convergence (here defined when the difference in likelihoods is sufficiently small)
        
        while(lld > 0.001){
          Nest <- vector(mode="list", length=length(year)) 
          irep <- irep + 1 # Keep track of which iteration this is, partly so ll_func knows we are inputting N
          for (iyear in 1:length(year)){
            Nest[[i]] <- rep(NA,nS[iyear])
            for(isite in 1:nS[iyear]){
              lowvals <- c(lowvals, 0)
              uppvals <- c(uppvals, 2500)
              # Find each N_i numerically
              temp <- try(uniroot(der_func,lower=0,upper=2500,
                                  model=fit1,
                                  isite=isite,iyear=iyear,
                                  bmat=bmat)$root,silent=TRUE)  
              dtemp <- 1 
              lowval <- c(0,0,0,0,.1,.1,.1,.1,.01,.01,.01,.01,.01)
              uppval <- c(1000,2500,5000,10000,1000,2500,5000,10000,1000,2500,5000,10000)
              while(class(temp) == "try-error" & dtemp < 9){       
                temp <- try(uniroot(der_func,lower=lowval[dtemp],upper=uppval[dtemp],
                                    model=fit1,
                                    isite=isite, iyear=iyear,
                                    bmat=bmat)$root,silent=TRUE)
                dtemp <- dtemp + 1 
                lowvals <- c(lowvals,lowval)
                uppvals <- c(uppvals,uppval)
              }
              Nest[[iyear]][isite] <- as.numeric(unlist(temp))
            }
            
            vals <- fit1$modelfit$allval
          }

          fit[[irep]] <- try(fit_model(irep=irep, Nguess=Nest, vals=vals),silent=FALSE)
          ll[irep] <- fit[[irep]]$ll.val
          fit1 <- fit[[irep]]
          lld <- abs(ll[irep]-ll[irep-1]) 
        } 
        fit1 <- fit[[irep]]
        fit1$iterations <- list(fit)
      }
      et <- proc.time()
      fit1$time <- (et-st)[3]
      if (is.na(fit1[1])){
        fit_k[[k]] <- NA
        fit_k.ll[k] <- NA
      } else{
        fit_k[[k]] <- fit1
        fit_k.ll[k] <- fit_k[[k]]$ll.val
      }
      
    }
    output <- list(fit_k[[min(c(1:nstart)[fit_k.ll==max(fit_k.ll,na.rm=T)],na.rm=T)]],fit_k,fit_k.ll)
    
    return(output)
  }
  
  
  # Fit the GAI model
  fit_model <- function(irep=1,Nguess=NULL,vals=NULL){
    
    if(irep==1){parm <- start_val_func()} else {parm <- vals}	
    
    if(a.type == "S") { meth <- "BFGS" } else {meth <- meth}
    
    ##pgtol is to improve precision of L-BFGS-B method (default is 1e-5)
    ## temp and tmax are for the SANN method
    this.fit <- try(optim(par=parm,
                          fn=ll_func,Nguess=Nguess,irep=irep,hessian=TRUE,method= meth,
                          lower = rep(-Inf, length(parm)), upper = rep(Inf,length(parm)),
                          control=list(trace=TRUE,maxit=50000, REPORT=10, temp=30, tmax=25,pgtol=1e-7,factr=1e-9)),
                    silent = FALSE)
    
    
    if(is.list(this.fit) & class(try(solve(this.fit$hessian),silent=TRUE))[1] != "try-error"){	
      # Model output
      betta.out <- vector(mode="list", length=length(year))
      afunc.out <- vector(mode="list", length=length(year))
      afunc.outNA <- vector(mode="list", length=length(year))
      N.out <- vector(mode="list", length=length(year))
      
      if (mu.type=="common" & vary_mu=="yes"){
        mu1.out <- vector(mode = "list", length=length(year))
      } else if (mu.type=="cov"|mu.type=="double")  {
        mu1.out <- vector(mode = "numeric", length=length(year))
      }else if (mu.type=="within"|mu.type=="doublewithin"|mu.type=="interaction"){
        tau.int.out <- vector(mode = "numeric", length=length(year))
        mu1.out <- vector(mode = "list", length=length(year))
        mu1.int.out <- vector(mode = "list", length=length(year))
      } else if (vary_mu=="no" & mu.type=="common"){
        mu1.out <- NULL
      }
      
      if(a.type=="N"){
        if(mu.type=="common"|mu.type=="cov"|mu.type=="double"){
          sigma.out <- vector(mode = "list", length = length(year))
          tau.int.out <- tau.slope.out <- tau.slope.out2 <- tau.slope.out3 <- NULL
          phi.int.out <- phi.out <- phi.slope.out <- phi.slope.out2 <-   NULL
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out <- mu1.slope2.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis  <- NULL
        } else if (mu.type=="within"){
          sigma.out <- vector(mode = "list", length = length(year))
          phi.int.out <- phi.out <- phi.slope.out <- phi.slope.out2 <- tau.slope.out2 <- tau.slope.out3 <-   NULL
          psi.out <- r.out  <- tau.slope.out <- mu1.slope.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis  <- NULL
        }else if (mu.type=="doublewithin"){
          sigma.out <- vector(mode = "list", length = length(year))
          phi.int.out <- phi.out <- phi.slope.out <- phi.slope.out2 <- tau.slope.out3 <-  NULL
          psi.out <- r.out  <- tau.slope.out <- tau.slope.out2 <-  mu1.slope.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis  <- NULL
        }else if (mu.type=="interaction"){
          sigma.out <- vector(mode = "list", length = length(year))
          phi.int.out <- phi.out <- phi.slope.out <- phi.slope.out2 <-   NULL
          psi.out <- r.out  <- tau.slope.out <- tau.slope.out2 <- tau.slope.out3 <-  mu1.slope.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis  <- NULL
        }
      }  
      
      if(a.type=="SO"){  
        tau.int.out <- tau.slope.out <- tau.slope.out2 <- tau.slope.out3 <- NULL
        
        if (phi_type=="variable" & sigma_type=="variable" & vary_r=="no"){
          phi.out <- vector(mode="numeric", length=length(year))
          sigma.out <- vector(mode = "list", length = length(year))
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out <- mu1.slope2.out <-  mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis  <- NULL
        }else if (phi_type=="variable" & sigma_type=="variable" & vary_r=="yes"){
          phi.out <- vector(mode="numeric", length=length(year))
          sigma.out <- vector(mode = "list", length = length(year))
          r.out <- vector(mode = "numeric", length = length(year))
          psi.out   <- mu1.int.out <-  mu1.slope.out <- mu1.slope2.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis  <- NULL
        }else if (phi_type=="variable" & sigma_type == "fixed") {
          phi.out <- vector(mode="numeric", length=length(year))
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out<- mu1.slope2.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <-sigma.out <-  w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- NULL
        } else if (phi_type=="fixed" & sigma_type == "variable") {
          sigma.out <- vector(mode = "list", length = length(year))
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out<- mu1.slope2.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- phi.out <- NULL
        } else if (phi_type=="fixed" & sigma_type == "fixed"){
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out<- mu1.slope2.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- sigma.out <-  w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- phi.out <- NULL
        }else if (phi_type=="slope" & sigma_type=="variable"& vary_r=="no"){
          phi.out <- vector(mode="numeric", length=length(year))
          sigma.out <- vector(mode = "list", length = length(year))
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out<- mu1.slope2.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- phi.int.out <- phi.slope.out <- NULL
        } else if (phi_type=="slope" & sigma_type=="variable"& vary_r=="yes"){
          phi.out <- vector(mode="numeric", length=length(year))
          sigma.out <- vector(mode = "list", length = length(year))
          r.out <- vector(mode = "numeric", length = length(year))
          psi.out  <- mu1.int.out <-  mu1.slope.out <- mu1.slope2.out<- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- phi.int.out <- phi.slope.out <- phi.slope.out2 <-  NULL
        }else if (phi_type=="slope" & sigma_type=="fixed"){
          phi.out <- vector(mode="numeric", length=length(year))
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out <- mu1.slope2.out<- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- sigma.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- phi.int.out <- phi.slope.out  <- NULL
        }else if (phi_type=="doubleslope" & sigma_type=="variable"){
          phi.out <- vector(mode="numeric", length=length(year))
          sigma.out <- vector(mode = "list", length = length(year))
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out<- mu1.slope2.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- phi.int.out <- phi.slope.out <- phi.slope.out2 <- NULL
        } else if (phi_type=="doubleslope" & sigma_type=="fixed"){
          phi.out <- vector(mode="numeric", length=length(year))
          psi.out <- r.out  <- mu1.int.out <-  mu1.slope.out <- mu1.slope2.out<- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- sigma.out <- w.out <- w.int.out <- w.slope.out  <- alpha.out  <- bsbasis <- phi.int.out <- phi.slope.out <- phi.slope.out2 <- NULL
        }
      }
      
      out.index <- 0
      
      if(a.type == "N" | a.type == "SO"){ 
        alpha.out <- NULL
        
        switch(mu.type,
               common = {
                 if (vary_mu=="yes"){
                   mu1.int.out <- vector(mode="numeric", length=length(year))
                   for (i in 1:length(year)){
                     mu1.out[[i]] <- rep(exp(this.fit$par[out.index+1]),nS[i]); mu1.int.out[i] <- this.fit$par[out.index+1]; out.index <- out.index + 1
                   }
                 }else {
                   mu1.out <- exp(this.fit$par[out.index+1]); mu1.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1
                 }
               }
               ,
               cov = {
                 for (i in 1:length(year)){
                   mu1.out[i] <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index+2]*mu.cov[i])}; 
                 mu1.int.out <- this.fit$par[out.index+1]; mu1.slope.out <- this.fit$par[out.index + 2]; 
                 out.index <- out.index + 2}
               ,
               double=
                 {
                   for (i in 1:length(year)){
                     mu1.out[i] <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index+2]*mu.cov[i] +this.fit$par[out.index+3]*mu.cov2[i])}; 
                   mu1.int.out <- this.fit$par[out.index+1]; mu1.slope.out <- this.fit$par[out.index + 2]; 
                   mu1.slope2.out <- this.fit$par[out.index + 3];
                   out.index <- out.index + 3}
               ,
               within={
                 for (i in 1:length(year)){
                   tau.int.out[i] <- this.fit$par[out.index+1];out.index <- out.index + 1
                 }
                 tau.slope.out <- this.fit$par[out.index + 1]; out.index <- out.index + 1
                 
                 for (i in 1:length(year)){
                   mu1.int.out[[i]] <- exp(tau.int.out[i] + tau.slope.out*tau.cov[[i]])
                 }
                 
                 for (i in 1:length(year)){
                   mu1.out[[i]] <- mu1.int.out[[i]] + this.fit$par[out.index+1]*mu.cov[i]
                 }
                 
                 mu1.slope.out <- this.fit$par[out.index+1];out.index <- out.index + 1
                 
               },
               doublewithin={
                 for (i in 1:length(year)){
                   tau.int.out[i] <- this.fit$par[out.index+1];out.index <- out.index + 1
                 }
                 tau.slope.out <- this.fit$par[out.index + 1]; out.index <- out.index + 1
                 tau.slope.out2 <- this.fit$par[out.index + 1]; out.index <- out.index + 1
                 
                 for (i in 1:length(year)){
                   mu1.int.out[[i]] <- exp(tau.int.out[i] + tau.slope.out*tau.cov[[i]] + tau.slope.out2*tau.cov2[[i]])
                 }
                 
                 for (i in 1:length(year)){
                   mu1.out[[i]] <- mu1.int.out[[i]] + this.fit$par[out.index+1]*mu.cov[i]
                 }
                 
                 mu1.slope.out <- this.fit$par[out.index+1];out.index <- out.index + 1
                 
               },
               interaction={
                 for (i in 1:length(year)){
                   tau.int.out[i] <- this.fit$par[out.index+1];out.index <- out.index + 1
                 }
                 tau.slope.out <- this.fit$par[out.index + 1]; out.index <- out.index + 1
                 tau.slope.out2 <- this.fit$par[out.index + 1]; out.index <- out.index + 1
                 tau.slope.out3 <- this.fit$par[out.index + 1]; out.index <- out.index + 1
                 
                 for (i in 1:length(year)){
                   mu1.int.out[[i]] <- exp(tau.int.out[i] + tau.slope.out*tau.cov[[i]] + tau.slope.out2*tau.cov2[[i]] + tau.slope.out3*tau.cov[[i]]*tau.cov2[[i]])
                 }
                 
                 for (i in 1:length(year)){
                   mu1.out[[i]] <- mu1.int.out[[i]] + this.fit$par[out.index+1]*mu.cov[i]
                 }
                 
                 mu1.slope.out <- this.fit$par[out.index+1];out.index <- out.index + 1
                 
               }
        )
        switch(B,
               "1" = {
                 if (sigma_type=="variable"){
                   for (i in 1:length(year)){
                     sigma.out[[i]] <- rep(exp(this.fit$par[out.index + 1]),2); out.index <- out.index + 1
                   }
                 }else{
                   sigma.out <- rep(exp(this.fit$par[out.index + 1]),2); out.index <- out.index + 1
                 }
               },
               "2" = {
                 switch(mu.diff.type,
                        common = {mu.diff.out <- rep(exp(this.fit$par[out.index + 1]),nS); mu.diff.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1},
                        cov = {mu.diff.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu.diff.cov); mu.diff.int.out <- this.fit$par[out.index+1]; mu.diff.slope.out <- this.fit$par[out.index+2]; out.index <- out.index + 2})
                 
                 switch(sigma.type,
                        hom = {sigma.out <- rep(exp(this.fit$par[out.index + 1]),2); out.index <- out.index + 1},
                        het = {sigma.out <- exp(this.fit$par[(out.index + 1):(out.index + 2)]); out.index <- out.index + 2})
                 
                 switch(w.type,
                        common = {w.out <- rep(expit(this.fit$par[out.index + 1]),nS); w.int.out <- this.fit$par[out.index+1];out.index <- out.index + 1},
                        cov = {w.out <- expit(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*w.cov); w.int.out <- this.fit$par[out.index+1]; w.slope.out <- this.fit$par[out.index+2];out.index <- out.index + 2})})
        
        if(a.type == "SO"){
          if (phi_type=="variable"){
            for (i in 1:length(year)){
              phi.out[i] <-  expit(this.fit$par[out.index+1]); out.index <- out.index + 1
              phi.int.out <- phi.slope.out <- phi.slope.out2 <-  NULL
            }
          }else if (phi_type=="fixed"){
            phi.out <-  expit(this.fit$par[out.index+1]); out.index <- out.index + 1
            phi.int.out <- phi.slope.out <- phi.slope.out2 <-   NULL
          }else if (phi_type=="slope"){
            phi.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1
            phi.slope.out <- this.fit$par[out.index+1]; out.index <- out.index + 1
            phi.slope.out2 <- NULL
            for (i in 1:length(year)){
              phi.out[i] <- expit(phi.int.out + phi.slope.out*phi.cov[i])
            }
          }else if (phi_type=="doubleslope"){
            phi.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1
            phi.slope.out <- this.fit$par[out.index+1]; out.index <- out.index + 1
            phi.slope.out2 <- this.fit$par[out.index+1]; out.index <- out.index + 1
            for (i in 1:length(year)){
              phi.out[i] <- expit(phi.int.out + phi.slope.out*phi.cov[i]+phi.slope.out2*phi.cov2[i])
            }
          }
        }
      } else if(a.type == "S"){
        if(dist.type == "P"){alpha.out <- this.fit$par} else {alpha.out <- this.fit$par[1:(length(this.fit$par)-1)]; out.index <- length(this.fit$par)-1}
        bsbasis <- bs(1:nT,df=degf,degree=deg,intercept=TRUE)
        afunc.out <- exp(matrix(bsbasis%*%alpha.out,nrow=nS,ncol=nT,byrow=TRUE))
        afunc.out <- afunc.out/rowSums(afunc.out)	
      }
      
      if(dist.type == "NB"){
        if (vary_r=="no"){
          r.out <- exp(this.fit$par[out.index+1])
        }else{
          for (i in 1:length(year)){
            r.out[i] <-  exp(this.fit$par[out.index+1]); out.index <- out.index + 1
          }
        }
      }
      if(dist.type == "ZIP"){
        psi.out <- expit(this.fit$par[out.index+1])}
      
      
      switch(a.type,
             N = {
               switch(B,
                      "1" = {
                        for (i in 1:length(year)){
                          if (sigma_type=="variable" & mu.type=="cov"){
                            afunc.out[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu1.out[i],sigma.out[[i]][1]),nrow=nS[i],ncol=nT,byrow=F)
                          } else if (sigma_type=="variable" & vary_mu=="yes" & mu.type=="common") {
                            afunc.out[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu1.out[[i]],sigma.out[[i]][1]),nrow=nS[i],ncol=nT,byrow=F)
                          }else if (sigma_type=="fixed" & mu.type=="cov") {
                            afunc.out[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu1.out[i],sigma.out[1]),nrow=nS[i],ncol=nT,byrow=F)
                          }else if (sigma_type=="fixed" & vary_mu=="yes" & mu.type=="common") {
                            afunc.out[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu1.out[[i]],sigma.out[1]),nrow=nS[i],ncol=nT,byrow=F)
                          }else if (sigma_type=="fixed" & mu.type=="within"|sigma_type=="fixed" & mu.type=="doublewithin"|sigma_type=="fixed" & mu.type=="interaction") {
                            afunc.out[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu1.out[[i]],sigma.out[1]),nrow=nS[i],ncol=nT,byrow=F)
                          }else if (sigma_type=="variable" & mu.type=="within"|sigma_type=="variable" & mu.type=="doublewithin"|sigma_type=="variable" & mu.type=="interaction") {
                            afunc.out[[i]] <- matrix(dnorm(rep(1:nT,each=nS[i]),mu1.out[[i]],sigma.out[[i]][1]),nrow=nS[i],ncol=nT,byrow=F)
                          }
                          
                          afunc.outNA[[i]] <- afunc.out[[i]]
                          afunc.outNA[[i]][is.na(y[[i]])] <- NA
                          N.out[[i]] <- apply(y[[i]],1,sum,na.rm=TRUE)/(apply(afunc.outNA[[i]],1,sum,na.rm=TRUE))
                        }
                      },
                      "2" = {
                        afunc.out <- matrix(rep(w.out,nT)*dnorm(rep(1:nT,nS),mu1.out,sigma.out[1]) + (1-rep(w.out,nT))*dnorm(rep(1:nT,nS),mu1.out + mu.diff.out,sigma.out[2]),nrow=nS,ncol=nT,byrow=TRUE)})},
             SO = {
               switch(B,
                      "1" = {	
                        if (sigma_type=="variable"&mu.type=="cov"|sigma_type=="variable"&mu.type=="double"){
                          for (i in 1:length(year)){
                            betta.out[[i]] <- matrix(c(pnorm(1,mean=rep(mu1.out[i],nS[i]),sd=sigma.out[[i]][1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu1.out[i],nS[i]*length(2:(nT-1))),
                                                             sd=sigma.out[[i]][1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu1.out[i],nS[i]*length(1:(nT-2))),
                                                               sd=sigma.out[[i]][1]),
                                                       1-pnorm(nT-1,
                                                               mean=rep(mu1.out[i],nS[i]),sd=sigma.out[[i]][1])),
                                                     nrow=nS[i],ncol=nT)}
                        }else if(sigma_type=="variable"&vary_mu=="yes" & mu.type=="common"){
                          for (i in 1:length(year)){
                            betta.out[[i]] <- matrix(c(pnorm(1,mean=mu1.out[[i]],sd=sigma.out[[i]][1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu1.out[[i]],length(2:(nT-1))),
                                                             sd=sigma.out[[i]][1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu1.out[[i]],length(1:(nT-2))),
                                                               sd=sigma.out[[i]][1]),
                                                       1-pnorm(nT-1,
                                                               mean=mu1.out[[i]],sd=sigma.out[[i]][1])),
                                                     nrow=nS[i],ncol=nT)}
                        }else if(sigma_type=="fixed"&mu.type=="cov"|sigma_type=="fixed"&mu.type=="double") {
                          for (i in 1:length(year)){
                            betta.out[[i]] <- matrix(c(pnorm(1,mean=rep(mu1.out[i],nS[i]),sd=sigma.out[1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu1.out[i],nS[i]*length(2:(nT-1))),
                                                             sd=sigma.out[1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu1.out[i],nS[i]*length(1:(nT-2))),
                                                               sd=sigma.out[1]),
                                                       1-pnorm(nT-1,
                                                               mean=rep(mu1.out[i],nS[i]),sd=sigma.out[1])),
                                                     nrow=nS[i],ncol=nT)}
                        }else if(sigma_type=="fixed" & vary_mu=="yes" & mu.type=="common"){
                          for (i in 1:length(year)){
                            betta.out[[i]] <- matrix(c(pnorm(1,mean=mu1.out[[i]],sd=sigma.out[1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu1.out[[i]],length(2:(nT-1))),
                                                             sd=sigma.out[1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu1.out[[i]],length(1:(nT-2))),
                                                               sd=sigma.out[1]),
                                                       1-pnorm(nT-1,
                                                               mean=mu1.out[[i]],sd=sigma.out[1])),
                                                     nrow=nS[i],ncol=nT)}
                        }else if(sigma_type=="fixed" & vary_mu=="no" & mu.type=="common"){
                          for (i in 1:length(year)){
                            betta.out[[i]] <- matrix(c(pnorm(1,mean=mu1.out,sd=sigma.out[1]),
                                                       pnorm(rep(2:(nT-1),each=nS[i]),
                                                             mean=rep(mu1.out,length(2:(nT-1))),
                                                             sd=sigma.out[1])-
                                                         pnorm(rep(1:(nT-2),each=nS[i]),
                                                               mean=rep(mu1.out,length(1:(nT-2))),
                                                               sd=sigma.out[1]),
                                                       1-pnorm(nT-1,
                                                               mean=mu1.out,sd=sigma.out[1])),
                                                     nrow=nS[i],ncol=nT)}
                        }
                      },
                      "2" = {
                        betta.out <- matrix(rep(w.out,nT)*c(pnorm(1,mean=mu1.out,sd=sigma.out[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out,length(2:(nT-1))),sd=sigma.out[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out,length(1:(nT-2))),sd=sigma.out[1]),1-pnorm(nT-1,mean=mu1.out,sd=sigma.out[1])) + (1-rep(w.out,nT))*c(pnorm(1,mean=mu1.out + mu.diff.out,sd=sigma.out[2]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out+ mu.diff.out,length(2:(nT-1))),sd=sigma.out[2])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out+ mu.diff.out,length(1:(nT-2))),sd=sigma.out[2]),1-pnorm(nT-1,mean=mu1.out+ mu.diff.out,sd=sigma.out[2])) ,nrow=nS,ncol=nT)
                      })
               
               
               
               for (i in 1:length(year)){
                 if (phi_type=="variable"|phi_type=="slope"|phi_type=="doubleslope"){
                   afunc.out[[i]] <- betta.out[[i]]				
                   for(j in 2:nT){
                     for(b in 1:(j-1)){
                       afunc.out[[i]][,j] <- afunc.out[[i]][,j] + betta.out[[i]][,b]*phi.out[i]^length(b:(j-1))
                     }
                   }
                   afunc.outNA[[i]] <- afunc.out[[i]]
                   afunc.outNA[[i]][is.na(y[[i]])] <- NA
                   N.out[[i]] <- apply(y[[i]],1,sum,na.rm=TRUE)/(apply(afunc.outNA[[i]],1,sum,na.rm=TRUE))
                 } else {
                   afunc.out[[i]] <- betta.out[[i]]				
                   for(j in 2:nT){
                     for(b in 1:(j-1)){
                       afunc.out[[i]][,j] <- afunc.out[[i]][,j] + betta.out[[i]][,b]*phi.out^length(b:(j-1))
                     }
                   }
                   afunc.outNA[[i]] <- afunc.out[[i]]
                   afunc.outNA[[i]][is.na(y[[i]])] <- NA
                   N.out[[i]] <- apply(y[[i]],1,sum,na.rm=TRUE)/(apply(afunc.outNA[[i]],1,sum,na.rm=TRUE))
                 }
               }
             })
      
      
      output <- list(ll.val=-this.fit$value,
                     npar=length(this.fit$par),
                     N.est=N.out,
                     w.est=w.out,
                     w.int=w.int.out,
                     w.slope=w.slope.out,
                     mu1.est=mu1.out,
                     mu1.int=mu1.int.out,
                     mu1.slope=mu1.slope.out,
                     mu1.slope2=mu1.slope2.out,
                     mu.diff.est=mu.diff.out,mu.diff.int=mu.diff.int.out,
                     mu.diff.slope=mu.diff.slope.out,
                     tau.int=tau.int.out,
                     tau.slope=tau.slope.out,
                     tau.slope2=tau.slope.out2,
                     tau.slope3=tau.slope.out3,
                     sigma.out=sigma.out,
                     r.out=r.out,
                     psi.est=psi.out,
                     phi.out=phi.out,
                     phi.int=phi.int.out,
                     phi.slope=phi.slope.out,
                     phi.slope2=phi.slope.out2,
                     afunc.out=afunc.out,
                     afunc.outNA=afunc.outNA,
                     betta.out=betta.out,
                     bsbasis=bsbasis,
                     alpha.out=alpha.out,
                     modelfit=list(Hessian=this.fit$hessian,
                                   starts=parm,
                                   allval=this.fit$par,
                                   nS=nS,nT=nT,
                                   convergence=this.fit$convergence,
                                   y=y))
      
      len_y <- vector(mode="numeric", length=length(year))
      for (i in 1:length(year)){
        output$Fitted[[i]] <- output$afunc.out[[i]]*output$N.est[[i]]
        if(dist.type == "ZIP") output$Fitted <- output$psi.est*output$Fitted
        output$dev[i] <- switch(dist.type,
                                P={
                                  2*(sum((output$modelfit$y[[i]]*log(output$modelfit$y[[i]]/output$Fitted[[i]])-(output$modelfit$y[[i]]-output$Fitted[[i]]))[!is.na(output$modelfit$y[[i]]) & output$modelfit$y[[i]]!=0])+sum(output$Fitted[[i]][output$modelfit$y[[i]]== 0 & !is.na(output$modelfit$y[[i]])]))},
                                NB={
                                  if (vary_r=="no"){
                                    2*(sum((output$modelfit$y[[i]]*log(output$modelfit$y[[i]]/output$Fitted[[i]])-((output$r.out+output$modelfit$y[[i]])*log((output$r.out+output$modelfit$y[[i]])/(output$r.out+output$Fitted[[i]]))))[output$modelfit$y[[i]]>0 & !is.na(output$modelfit$y[[i]])])-sum((output$r.out*log(output$r.out/(output$r.out+output$Fitted[[i]])))[!is.na(output$modelfit$y[[i]]) & output$modelfit$y[[i]]==0]))
                                  }else{
                                    2*(sum((output$modelfit$y[[i]]*log(output$modelfit$y[[i]]/output$Fitted[[i]])-((output$r.out[i]+output$modelfit$y[[i]])*log((output$r.out[i]+output$modelfit$y[[i]])/(output$r.out[i]+output$Fitted[[i]]))))[output$modelfit$y[[i]]>0 & !is.na(output$modelfit$y[[i]])])-sum((output$r.out[i]*log(output$r.out[i]/(output$r.out[i]+output$Fitted[[i]])))[!is.na(output$modelfit$y[[i]]) & output$modelfit$y[[i]]==0]))
                                  }
                                },
                                ZIP={NA})
        len_y[i] <- length(output$modelfit$y[[i]][!is.na(output$modelfit$y[[i]])])
      } 
      output$D <- sum(output$dev)/(sum(len_y)-output$npar)
      
      output
    } else {NA}
  }
}







######## end duplicate the model file inplace !!!

        # Output is a list of length 3: the best model of multiple starts, output from the multiple starts, and the log-likelihood values for the multiple starts
        
        output <- fit_it_model()

     }   


        
        
        
        
      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

     
