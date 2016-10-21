### TSHT.R
### Function: Implements the two-stage hard thresholing
###           procedure for estimation and inference of 
###           treatment effects in the presence of invalid IVs
###           as developed in Guo, Kang, Cai, and Small (2016)             
### Maintainer: Hyunseung Kang
### E-mail: hskang@stanford.edu


### TSHT.hdim
### Function: Point estimate and CI for treatment effect with invalid IVs ###
###           using individual-level data                                 ###
### Input: Y, continuous outcome (n by 1 vector)
###        D, continuous or discrete treatment (n by 1 vector)
###        Z, continuous or discrete instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continuous or discrete covariate matrix containing 
###           p_x covarites (n by p_x matrix)
###        intercept, should the intercept term be included? (TRUE)
###        alpha, a numeric value indicating the significance level for 
###               confidence interval
###        estimator, choose either Javanmard & Montanari or HDI package
### Output: a list of (a) beta (scalar numeric value:
###                             the estimate of the treatment effect)
###                   (b) se (scalar numeric value:
###                           the standard error of beta)
###                   (c) ci (two dimensional vector, 
###                           the 1-alpha confidence interval for beta
###                           with the lower and upper endpoints)
###                   (d) V (numeric vector denoting the set of valid and relevant IVs)
###                   (e) S (numeric vector denoting the set of relevant IVs)
TSHT.hdim <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,estimator=c("hdi","JM")){
  # Constants
  n = length(Y)
  
  # Setting Up W Matrix
  if(missing(X)) {
    W = Z
  } else{
    W = cbind(Z,X)
  }
  #if(intercept) {
  #  W = cbind(W,1)
  #}
  
  # First Part (ITT Estimation)#
  if(estimator[1] == "hdi") {
    model_Y = lasso.proj(W,Y)
    ITT_Y = model_Y$bhat[1:ncol(Z)]
    residY = (Y - W %*% model_Y$betahat)
    SigmaSqY = sum(residY^2)/n
    
    model_D = lasso.proj(W,D)
    ITT_D = model_D$bhat[1:ncol(Z)]
    residD = (D - W %*% model_D$betahat)
    SigmaSqD = sum(residD^2)/n
    
    SigmaYD = sum(residD * residY)/n
    
    # ZIJIAN: How do you estimate this from the hdi comment
    # How do you estimate covariance matrix of $W$?
    covW = t(W) %*% W / n
  }
  if(estimator[1] == "JM") {
    pz=ncol(Z)
    model_Y <- SSLasso(W,Y,intercept,verbose=FALSE)
    model_D = SSLasso(W,D,intercept,verbose=FALSE) 
    # ZIJIAN: How do you estimate this from Javanmard's code? Should we extract the $M$ from InverseLinfinity?
    #covW = t(W) %*% W / n
    ITT_Y = model_Y$unb.coef[1:(pz)]
    ITT_D = model_D$unb.coef[1:(pz)]
    if(intercept==TRUE){
    covWinv =(model_D$proj %*% t(model_D$proj)/n)[2:(pz+1),2:(pz+1)]
    W.whole<-cbind(rep(1,n),W)
    covW<- t(W.whole)%*%(W.whole)/n
    resid_Y<-Y-W.whole%*%model_Y$coef
    resid_D<-D-W.whole%*%model_D$coef
    Resid<-cbind(resid_Y,resid_D)
    Theta.hat<-t(Resid)%*%Resid/n
    SigmaSqY=Theta.hat[1,1]
    SigmaSqD=Theta.hat[2,2]
    SigmaYD =Theta.hat[1,2]
  } else{
    covWinv =(model_D$proj %*% t(model_D$proj)/n)[1:(pz),1:(pz)]
    W.whole<-W
    covW<- t(W.whole)%*%(W.whole)/n
    resid_Y<-Y-W.whole%*%model_Y$coef
    resid_D<-D-W.whole%*%model_D$coef
    Resid<-cbind(resid_Y,resid_D)
    Theta.hat<-t(Resid)%*%Resid/n
    SigmaSqY=Theta.hat[1,1]
    SigmaSqD=Theta.hat[2,2]
    SigmaYD =Theta.hat[1,2]
  }
  return(TSHT.summ(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD = SigmaSqD,SigmaSqY = SigmaSqY,SigmaYD=SigmaYD,covW=covW,covWinv=covWinv,n=n,alpha=alpha))
}

### TSHT.ldim
### Function: Point estimate and CI for treatment effect with invalid IVs ###
###           using individual-level data                                 ###
### Input: Y, continuous outcome (n by 1 vector)
###        D, continuous or discrete treatment (n by 1 vector)
###        Z, continuous or discrete instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continous or discrete covariate matrix containing 
###           p_x covarites (n by p_x matrix)
###        intercept, should the intercept term be included? (TRUE)
###        alpha, a numeric value indicating the significance level for 
###               confidence interval
### Output: a list of (a) beta (scalar numeric value:
###                             the estimate of the treatment effect)
###                   (b) se (scalar numeric value:
###                           the standard error of beta)
###                   (c) ci (two dimensional vector, 
###                           the 1-alpha confidence interval for beta
###                           with the lower and upper endpoints)
###                   (d) V (numeric vector denoting the set of valid and relevant IVs)
###                   (e) S (numeric vector denoting the set of relevant IVs)
TSHT.ldim <- function(Y,D,Z,X,intercept=FALSE,alpha=0.05) {
  # Constants
  n = length(Y)
  
  # Setting Up W Matrix
  if(missing(X)) {
    W = Z
  } else{
    W = cbind(Z,X)
  }
  if(intercept) {
    W = cbind(W,1)
  }
  
  # First Part (OLS Estimation)
  qrW = qr(W)
  ITT_Y = qr.coef(qrW,Y)[1:ncol(Z)]
  ITT_D = qr.coef(qrW,D)[1:ncol(Z)]
  SigmaSqY = sum(qr.resid(qrW,Y)^2)/n
  SigmaSqD = sum(qr.resid(qrW,D)^2)/n
  SigmaYD = sum(qr.resid(qrW,Y) * qr.resid(qrW,D)) / n
  covW = t(W) %*% W / n
  covWinv = solve(covW)
  pz=ncol(Z)
  # Second Part (Invalid IV Estimation)
  return(TSHT.summ(ITT_Y = ITT_Y,ITT_D = ITT_D,SigmaSqD = SigmaSqD,SigmaSqY = SigmaSqY,SigmaYD=SigmaYD,covW=covW,covWinv = covWinv,n=n,alpha=alpha))
}

### TSHT.summ
### Function: Point estimate and CI for treatment effect with invalid IVs ###
###           using summary statistics 
### Input: ITT_Y, a p_z by 1 vector of estimated values of each IV's 
###               effect on the outcome
###        ITT_D, a p_z by 1 vector of estimated values of each IV's 
###               effect on the treatment
###        SigmaSqD, a positive scalar value of the estimated 
###               variance of the treatment
###        SigmaSqY, a positive scalar value of the estimated 
###               variance of the outcome
###        SigmaYD, a positive scalar value of the estimated 
###                 covariance between the outcome and treatment
###        covW, the estimated p by p covariance matrix of the 
###              concatenated matrix consisting of p_z IVs and 
###              p_x covariates (i.e. p by p matrix where
###              p = p_z + px)
###        covWinv, the estimated p by p inverse covariance 
###                 of the concatenated matrix consisting of 
###                 p_z IVs and p_x covariats (i.e. p by p
###                 matrix where p = p_z + p_x).
###                 Either covW or covWinv must be provided.
###        n, sample size used in the estimates above
###        weightedBeta, should we weigh the beta? Only applies to case where n > p and covW is provided
###        alpha, a numeric value indicating the significance level for 
###               confidence interval
###        tuning, a constant value greater than 2
### Output: a list of (a) beta (scalar numeric value:
###                             the estimate of the treatment effect)
###                   (b) se (scalar numeric value:
###                           the standard error of beta)
###                   (c) ci (two dimensional vector, 
###                           the 1-alpha confidence interval for beta
###                           with the lower and upper endpoints)
###                   (d) V (numeric vector denoting the set of valid and relevant IVs)
###                   (e) S (numeric vector denoting the set of relevant IVs)
TSHT.summ <- function(ITT_Y,ITT_D,SigmaSqD,SigmaSqY,SigmaYD,covW,covWinv,n,weightedBeta = FALSE,alpha=0.05,tuning=2.01) {
  pz = length(ITT_Y)
  # ZIJIAN: Why do you distinguish the intercept here? Can't we residualize the intercept by de-meaning the estimates?
  #if(missing(covWinv)&intercept==TRUE) {
  #  covWinv = solve(covW)[2:(pz+1),2:(pz+1)]
  #}
  #if(missing(covWinv)&intercept==FALSE) {
  #  covWinv = solve(covW)[1:(pz),1:(pz)]
  #}
  
  # First Stage
  Stilde = (1:length(ITT_D))[(abs(ITT_D) >= sqrt(SigmaSqD) * sqrt(diag(covWinv)) * sqrt(tuning*log(max(pz,n))/n))]
  
  # Second Stage
  invalidIV.candidate = matrix(FALSE,length(ITT_Y),length(Stilde))
  pi.candidate = matrix(0,length(ITT_Y),length(Stilde))
  for(i in 1:length(Stilde)) {
    j = Stilde[i]
    betaj = ITT_Y[j] / ITT_D[j]
    pi.candidate.noHT = ITT_Y - ITT_D * betaj
    sigmaSq.candidate = SigmaSqY + betaj^2 * SigmaSqD - 2*betaj*SigmaYD #Is it possible to use reduced-form sigma directly?
	pi.candidate.HT.index = (1:length(ITT_Y))[(abs(pi.candidate.noHT) >= tuning * sqrt(sigmaSq.candidate) * 
                                                 sqrt(diag(covWinv)- 2 *ITT_D / ITT_D[j] * covWinv[j,] + ITT_D^2 / ITT_D[j]^2 * covWinv[j,j]) * 
                                                 sqrt(log(max(pz,n))/n))]
    invalidIV.candidate[intersect(Stilde,pi.candidate.HT.index),i] = TRUE
    pi.candidate[,i] = 0; 
    pi.candidate[invalidIV.candidate[,i],i] = pi.candidate.noHT[invalidIV.candidate[,i]]
  }
  
  # L0 Penalization
  invalidIV.candidate.l0 = colSums(invalidIV.candidate)
  jstar = which(min(invalidIV.candidate.l0) == invalidIV.candidate.l0)
  if(length(jstar) > 1) {
    # L1 Penalization if non-unique form exists
    invalidIV.candidate.l1 = colSums(abs(pi.candidate[,jstar]))
    jstar = jstar[which(min(invalidIV.candidate.l1) == invalidIV.candidate.l1)]
    if(length(jstar) > 1) jstar = jstar[1]
  }

  # Vtilde and Aweight
  Vtilde = (Stilde)[!(Stilde %in% which(invalidIV.candidate[,jstar] == TRUE))]
  if(ncol(covW) < (n/10) & weightedBeta & !missing(covW)) {
    Aweight = covW[Vtilde,Vtilde] - 
              covW[Vtilde,-Vtilde] %*% (solve(covW))[-Vtilde,-Vtilde] %*% covW[-Vtilde,Vtilde]
    
    ### betaE, standard error, and conf int.
    commonDenom = (t(ITT_D[Vtilde]) %*% Aweight %*% ITT_D[Vtilde])
    betaE = (t(ITT_D[Vtilde]) %*% Aweight %*% ITT_Y[Vtilde]) / commonDenom
    VE = (SigmaSqY + betaE^2 * SigmaSqD - 2*betaE*SigmaYD)/ commonDenom
    confInt = c(betaE - qnorm(1-alpha/2) * sqrt(VE/n),betaE + qnorm(1-alpha/2) * sqrt(VE/n))
    return(list(beta=betaE,se = sqrt(VE/n),ci = confInt,V=Vtilde,S = Stilde))  
  } else {
    Aweight = diag(1,length(Vtilde),length(Vtilde))
    commonDenom = (t(ITT_D[Vtilde]) %*% Aweight %*% ITT_D[Vtilde])
    
    ### General beta estimator betaG, standard error and 
    betaG = (t(ITT_D[Vtilde]) %*% Aweight %*% ITT_Y[Vtilde]) / commonDenom
    VG = (SigmaSqY + betaG^2 * SigmaSqD - 2*betaG*SigmaYD)*(t(ITT_D[Vtilde]) %*% covWinv[Vtilde,Vtilde]%*%ITT_D[Vtilde])/ (commonDenom)^2
    confInt = c(betaG - qnorm(1-alpha/2) * sqrt(VG/n),betaG + qnorm(1-alpha/2) * sqrt(VG/n))
    return(list(beta=betaG,se = sqrt(VG/n),ci = confInt,V=Vtilde,S = Stilde))
  }
}
  
