### TSHT.R
### Function: Implements the two-stage hard thresholing
###           procedure for estimation and inference of 
###           treatment effects in the presence of invalid IVs
###           as developed in Guo, Kang, Cai, and Small (2016)             
### Maintainer: Hyunseung Kang
### E-mail: hskang@stanford.edu

### InputCheck
### FUNCTION: checks the input for TSHT() to make sure they are properly
###           formatted. This is an alternative to the model.matrix()
###           command.
### INPUT: Y, continuous outcome vector (u by 1 vector)
###        D, continuous or discrete treatment vector (n by 1 vector)
###        Z, continuous or discrete instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continuous or discrete matrix containing p_x 
###           covariates (n by p_x matrix)
### OUTPUT: a list containing
###           Y, a numeric outcome vector (n by 1 vector)
###           D, a numeric treatment vector (n by 1 vector)
###           Z, a numeric instrument matrix containing p_z instruments
###              (n by p_z matrix)
###           X, a numeric covariate matrix containing p_x
###              covariates (n by p_x matrix). 
###              If X is missing, a NULL value is returned.
InputCheck <- function(Y,D,Z,X) {
  # Check Y
  stopifnot(missing(Y),(is.numeric(Y) || is.logical(Y)))
  if((is.matrix(Y) || is.data.frame(Y)) & ncol(Y) > 1) stop("Y must be a vector!")
  Y = as.numeric(Y)
  Y.NAorNANindex = (is.na(Y) | is.nan(Y))
  
  # Check D
  stopifnot(missing(D),(is.numeric(D) || is.logical(D) || is.character(D)))
  if((is.matrix(D) || is.data.frame(D)) & ncol(D) > 1) stop("D must be a vector!")
  if(is.character(D) || is.factor(D)) {
    print("Treatment vector is factor or character-type. 
          Checking to see if there are only two unique character values.")
    tempD = factor(D)
    if(nlevels(tempD) != 2) {
      stop("There aren't two unique values!")
    } else{
      print("There are two unique values. Converting them to a binary vector...")
      print(paste("--->",levels(tempD)[1],"is assigned 0"))
      print(paste("--->",levels(tempD)[2],"is assigned 1"))
      D = as.numeric(tempD) - 1
    }
  } 
  D = as.numeric(D)
  D.NAorNANindex = (is.na(D) | is.nan(D))
  
  if(missing(X)) {
    return(list(Y =Y,D=D,Z=Z,X=NULL))
  } else{
    return(list(Y =Y,D=D,Z=Z,X=X))
  }
}


### TSHT
### FUNCTION: Point estimate, SE, and CI for treatment effect with invalid IVs 
###           with individual-level data using 
###           two-stage hard thresholding (TSHT)
### INPUT: Y, continuous outcome vector (u by 1 vector)
###        D, continuous or discrete treatment vector (n by 1 vector)
###        Z, continuous or discrete instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continuous or discrete matrix containing p_x 
###           covariates (n by p_x matrix)
###        intercept, should the intercept term be included? 
###                   (TRUE/FALSE)
###        alpha, a numeric value indicating the significance level for 
###               the confidence interval
###        tuning, a numeric value tuning parameter for TSHT greater 
###                than 2
###        ..., additional arguments to be passed to SSLasso and USolver
###             in JMCode.r
### OUTPUT: a list (a) beta (scalar numeric value:
###                             the estimate of the treatment effect)
###                (b) se (scalar numeric value:
###                           the standard error of beta)
###                (c) p-value (scalar numeric value:
###                           p-value of no treatment effect for beta)
###                (d) ci (two dimensional vector, 
###                           the 1-alpha confidence interval for beta
###                           with the lower and upper endpoints)
###                (e) V (numeric vector denoting the set of valid and relevant IVs)
###                (f) S (numeric vector denoting the set of relevant IVs)
TSHT <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,tuning=2.01,...) {
  # Check and Clean Input Type #
  # Check Y
  stopifnot(missing(Y),(is.numeric(Y) || is.logical(Y)))
  if((is.matrix(Y) || is.data.frame(Y)) & ncol(Y) > 1) stop("Y must be a vector!")
  Y = as.numeric(Y)
  
  # Check D
  stopifnot(missing(D),(is.numeric(D) || is.logical(D) || is.character(D)))
  if((is.matrix(D) || is.data.frame(D)) & ncol(D) > 1) stop("D must be a vector!")
  if(is.character(D) || is.factor(D)) {
    print("Treatment vector is factor or character-type. 
          Checking to see if there are only two unique character values.")
    tempD = factor(D)
    if(nlevels(tempD) != 2) {
      stop("There aren't two unique values!")
    } else{
      print("There are two unique values. Converting them to a binary vector...")
      print(paste("--->",levels(tempD)[1],"is assigned 0"))
      print(paste("--->",levels(tempD)[2],"is assigned 1"))
      D = as.numeric(tempD) - 1
    }
  } 
  D = as.numeric(D)
  
  # Check Z and X. 
  # Currently, we leverage the model.matrix function. 
  # While it's slow, it's easier to read and robust towards many use cases
  if(!missing(X)) {
    model.matrix(~.,data=data.frame(Y,D,Z,X))[,-1]
  } else {
    model.matrix(~.,data=data.frame(Y,D,Z))[,-1]
  }
  
  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) != 1,alpha <= 1,alpha >= 0)
  stopifnot(is.numeric(tuning),length(tuning) != 1, tuning >=2)
  
  # Constants
  n = length(Y); pz=ncol(Z)
  
  # Setting Up W Matrix
  if(missing(X)) {
    W = Z
  } else{
    W = cbind(Z,X)
  }
  p = ncol(W)
}

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
TSHT.hdim <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,tuning=2.01){
  # Constants
  n = length(Y); pz=ncol(Z)
  
  # Setting Up W Matrix
  if(missing(X)) {
    W = Z
  } else{
    W = cbind(Z,X)
  }
  p = ncol(W)
  
  # Fit Reduced-Form Model for Y and D
  model_Y <- SSLasso(X=W,y=Y,intercept=intercept,verbose=FALSE)
  model_D = SSLasso(X=W,y=D,intercept=intercept,verbose=FALSE) 
  ITT_Y = model_Y$unb.coef[1:(pz)]
  ITT_D = model_D$unb.coef[1:(pz)]
  resid_Y = model_Y$resid.lasso; resid_D = model_D$resid.lasso
  SigmaSqY=sum(resid_Y^2)/n
  SigmaSqD=sum(resid_D^2)/n
  SigmaYD =sum(resid_Y * resid_D)/n
  
  WUMat = model_D$WUMat
  if(intercept) {
    W = cbind(W,1)
    covW = t(W) %*% W/n
  } else{
    covW = t(W) %*% W/n
  }
  
  return(TSHT.helper(ITT_Y = ITT_Y,ITT_D = ITT_D,
                     SigmaSqD = SigmaSqD,SigmaSqY = SigmaSqY,SigmaYD=SigmaYD,
                     covW = covW,WUMat = WUMat,
                     weightedBeta=FALSE,alpha=alpha,tuning=tuning))
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
TSHT.ldim <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,tuning=2.01) {
  # Constants
  n = length(Y); pz=ncol(Z)
  
  # Setting Up W Matrix
  if(missing(X)) {
    W = Z
  } else{
    W = cbind(Z,X)
  }
  p = ncol(W)
  
  # Include intercept
  if(intercept) {
    W = cbind(W,1)
  } 
  
  # Compute covariance of W and W %*% U
  covW = t(W) %*% W /n
  WUMat = W %*% (solve(covW))[,1:p]
  
  # First Part (OLS Estimation)
  qrW = qr(W)
  ITT_Y = qr.coef(qrW,Y)[1:pz]
  ITT_D = qr.coef(qrW,D)[1:pz]
  SigmaSqY = sum(qr.resid(qrW,Y)^2)/n
  SigmaSqD = sum(qr.resid(qrW,D)^2)/n
  SigmaYD = sum(qr.resid(qrW,Y) * qr.resid(qrW,D)) / n

  # Second Part (Invalid IV Estimation)
  TSHT.helper(ITT_Y = ITT_Y,ITT_D = ITT_D,
              SigmaSqD = SigmaSqD,SigmaSqY = SigmaSqY,SigmaYD=SigmaYD,
              covW = covW,WUMat = WUMat,weightedBeta = TRUE,alpha=alpha,tuning=tuning)
}

### TSHT.helper
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
###        covW, estimmate of E(W^T W), used for weighting beta 
###        WUMat, the matrix W %*% U in the thresholding procedure.
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
TSHT.helper <- function(ITT_Y,ITT_D,
                        SigmaSqD,SigmaSqY,SigmaYD,
                        covW,WUMat,
                        weightedBeta = FALSE,alpha=0.05,tuning=2.01) {
  n = nrow(WUMat)
  pz = length(ITT_Y)
  
  # First Stage
  Stilde = (1:length(ITT_D))[abs(ITT_D) >= (sqrt(SigmaSqD) * sqrt(colSums(WUMat^2)/n) * sqrt(tuning*log(max(pz,n))/n))]
  if(length(Stilde) == 0) {
    print("The IVs are individually weak and TSHT using these IVs will give misleading CIs, SEs, and p-values.")
    response = readline("Proceed anyway treating all the IVs as strong (yes/no)?: ")
    if(response == "yes") {
      Stilde = 1:length(ITT_D)
      print("Again, please be cautious interpreting the CIs, SEs, and p-values...")
    }
  }
  suppgamma = rep(FALSE,length(ITT_D))
  suppgamma[Stilde] = TRUE
  
  # Second Stage
  # pi.candidate is the estimated value of pi across different candidates
  pi.candidate = matrix(0,length(ITT_Y),length(Stilde)) 
  for(i in 1:length(Stilde)) {
    j = Stilde[i]
    betaj = ITT_Y[j] / ITT_D[j]
    pi.candidate.noHT = ITT_Y - ITT_D * betaj #HT stands for hard thresholding
    sigmaSq.candidate = SigmaSqY + betaj^2 * SigmaSqD - 2*betaj*SigmaYD
	  pi.candidate.HT.index = (1:length(ITT_Y))[suppgamma & #this line here essentially does k \in Stilde in our paper
	                                            (abs(pi.candidate.noHT) >= tuning * sqrt(sigmaSq.candidate) * 
	                                              sqrt(colSums( (WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) *
                                                sqrt(log(max(pz,n))/n))]
    pi.candidate[pi.candidate.HT.index,i] = pi.candidate.noHT[pi.candidate.HT.index]
  }
  
  # L0 Penalization
  pi.candidate.l0 = colSums(abs(pi.candidate) > 0) #The inner operation doesn't turn abs(pi.candidate) > 0 into a vector.
  jstar = which(min(pi.candidate.l0) == pi.candidate.l0)
  if(length(jstar) > 1) {
    # L1 Penalization if non-unique form exists
    pi.candidate.l1 = colSums(abs(pi.candidate[,jstar]))
    minL1amongNonUniqueL0 = which(min(pi.candidate.l1) == pi.candidate.l1) 
    jstar = jstar[minL1amongNonUniqueL0]
    if(length(jstar) > 1) jstar = jstar[1] #If it's still non-unique, we just choose the first one.
  }
  
  # Estimated support of pi 
  supppi = which(abs(pi.candidate[,jstar]) > 0) 
  
  # Vtilde and Aweight
  Vtilde = setdiff(Stilde,supppi) #this is a setminus operation with Stilde \setminus supp(Pi)
  if((ncol(covW) < (n/10) || weightedBeta) && !missing(covW)) {
    VtildeComplement = setdiff(1:ncol(covW),Vtilde) #Remember, there can be an intercept in covW and the intercept term
                                                    #is always the last term of W
    if(length(VtildeComplement) == 0) {
      Aweight = covW[Vtilde,Vtilde]
    } else{
      Aweight = covW[Vtilde,Vtilde] - 
                (covW[Vtilde,VtildeComplement,drop=FALSE] %*% 
                (solve(covW[VtildeComplement,VtildeComplement])) %*% 
                covW[VtildeComplement,Vtilde,drop=FALSE])
    }
    commonDenom = (t(ITT_D[Vtilde]) %*% Aweight %*% ITT_D[Vtilde])
    
    ### betaE, standard error, and conf int.
    betaE = (t(ITT_D[Vtilde]) %*% Aweight %*% ITT_Y[Vtilde]) / commonDenom
    sigmaSq = SigmaSqY + betaE^2 * SigmaSqD - 2*betaE*SigmaYD
    V = sigmaSq / commonDenom
    confInt = c(betaE - qnorm(1-alpha/2) * sqrt(V/n),betaE + qnorm(1-alpha/2) * sqrt(V/n))
    return(list(beta=betaE,se = sqrt(V/n),ci = confInt,V=Vtilde,S = Stilde))  
  } else {
    commonDenom = sum(ITT_D[Vtilde]^2)
    
    ### General beta estimator betaG
    betaG = sum(ITT_D[Vtilde] * ITT_Y[Vtilde]) / commonDenom
    sigmaSq = SigmaSqY + betaG^2 * SigmaSqD - 2*betaG*SigmaYD
    V = sigmaSq *
         sum((WUMat[,Vtilde] %*% ITT_D[Vtilde] / sqrt(n))^2) / (commonDenom)^2
    confInt = c(betaG - qnorm(1-alpha/2) * sqrt(V/n),betaG + qnorm(1-alpha/2) * sqrt(V/n))
    return(list(beta=betaG,se = sqrt(V/n),ci = confInt,V=Vtilde,S = Stilde))
  }
}
  
