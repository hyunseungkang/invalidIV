### TSHT.R
### Function: Implements the two-stage hard thresholing
###           procedure for estimation and inference of 
###           treatment effects in the presence of invalid IVs
###           as developed in Guo, Kang, Cai, and Small (2017)             
### Maintainer: Hyunseung Kang
### E-mail: hskang@stanford.edu

### TSHT
### FUNCTION: Point estimate, SE, and CI for treatment effect with invalid IVs 
###           using a single-sample, individual-level data via TSHT. 
### INPUT: Y, continuous, non-missing, numeric outcome vector (u by 1 vector)
###        D, continuous or discrete, non-missing, numeric treatment vector (n by 1 vector)
###        Z, continuous or discrete, non-missing, numeric instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continuous or discrete, non-missing, numeric matrix containing p_x 
###           covariates (n by p_x matrix)
###        intercept, a boolean scalar asking "should the intercept term be included?" 
###                   with TRUE/FALSE (default = TRUE)
###        alpha, a numeric scalar value between 0 and 1 indicating the significance level for 
###               the confidence interval (default = 0.05)
###        tuning, a numeric scalar value tuning parameter for TSHT greater 
###                than 2 (default = 2.01)
###        method, a character scalar declaring the method used to estimate the inputs in TSHT
###                (default = "OLS")
### OUTPUT: a list (a) VHat (numeric vector denoting the set of valid and relevant IVs)
###                (b) SHat (numeric vector denoting the set of relevant IVs)
###                (c) betaHat (scalar numeric value:
###                             the estimate of the treatment effect)
###                (d) varHat (scalar numeric value:
###                            estimated variance of betaHat)
###                (e) ci (two dimensional vector, 
###                           the 1-alpha confidence interval for beta
###                           with the lower and upper endpoints)
TSHT <- function(Y,D,Z,X,intercept=TRUE,alpha=0.05,tuning=2.01,method=c("OLS","DeLasso")) {
  # Check and Clean Input Type #
  # Check Y
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  Y = as.numeric(Y)
  
  # Check D
  stopifnot(!missing(D),(is.numeric(D) || is.logical(D)),(is.matrix(D) || is.data.frame(D)) && ncol(D) == 1)
  stopifnot(all(!is.na(D)))
  D = as.numeric(D)

  # Check Z
  stopifnot(!missing(Z),(is.numeric(Z) || is.logical(Z)),is.matrix(Z))
  stopifnot(all(!is.na(Z)))

  # Check dimesions
  stopifnot(length(Y) == length(D), length(Y) == nrow(Z))

  # Check X, if present
  if(!missing(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),is.matrix(X) && nrow(X) == nrow(Z))
    stopifnot(all(!is.na(X)))

    W = cbind(Z,X)
  } else {
    W = Z
  }

  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.numeric(alpha),length(alpha) == 1,alpha <= 1,alpha >= 0)
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  stopifnot(is.character(method))

  # Derive Inputs for TSHT
  n = length(Y); pz=ncol(Z)
  if(method[1] == "OLS") {
    inputs = TSHT.OLS(Y,D,W,pz,intercept)
    A = t(inputs$WUMat) %*% inputs$WUMat / n
  } else{
    inputs = TSHT.DeLasso(Y,D,W,pz,intercept)
    A = diag(pz)
  }

  # Estimate Valid IVs
  SetHats = TSHT.VHat(ITT_Y = inputs$ITT_Y,ITT_D = inputs$ITT_D,WUMat = inputs$WUMat,
                      SigmaSqD = inputs$SigmaSqD,SigmaSqY = inputs$SigmaSqY,SigmaYD=inputs$SigmaYD,tuning=tuning)
  VHat = SetHats$VHat; SHat = SetHats$SHat

  # Obtain point est, se, and ci
  AVHat = solve(A[VHat,VHat])
  betaHat = (t(inputs$ITT_Y[VHat]) %*% AVHat %*% inputs$ITT_D[VHat]) / (t(inputs$ITT_D[VHat]) %*% AVHat %*% inputs$ITT_D[VHat])
  SigmaSq = inputs$SigmaSqY + betaHat^2 * inputs$SigmaSqD - 2*betaHat * inputs$SigmaYD
  betaVarHat = SigmaSq * (t(inputs$ITT_D[VHat]) %*% AVHat %*% (t(inputs$WUMat) %*% inputs$WUMat/ n)[VHat,VHat] %*% AVHat %*% inputs$ITT_D[VHat]) / (t(inputs$ITT_D[VHat]) %*% AVHat %*% inputs$ITT_D[VHat])^2
  ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat / n),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat/n))

  return(list(VHat = VHat,SHat=SHat,betaHat=betaHat,betaVarHat = betaVarHat,ci=ci))
}

TSHT.OLS <- function(Y,D,W,pz,intercept=TRUE) {
  # Include intercept
  if(intercept) {
    W = cbind(W,1)
  }
  p = ncol(W) 
  
  # Compute covariance of W and W %*% U
  covW = t(W) %*% W /n #this should automatically turn covW into a matrix
  WUMat = W %*% (solve(covW))[,1:pz]
  
  # First Part (OLS Estimation)
  qrW = qr(W)
  ITT_Y = qr.coef(qrW,Y)[1:pz]
  ITT_D = qr.coef(qrW,D)[1:pz]
  SigmaSqY = sum(qr.resid(qrW,Y)^2)/(n -p)
  SigmaSqD = sum(qr.resid(qrW,D)^2)/(n -p)
  SigmaYD = sum(qr.resid(qrW,Y) * qr.resid(qrW,D)) / (n - p)

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
}

TSHT.DeLasso <- function(Y,D,W,pz,intercept=TRUE) {
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

  return(list(ITT_Y = ITT_Y,ITT_D = ITT_D,WUMat = WUMat,SigmaSqY = SigmaSqY,SigmaSqD = SigmaSqD,SigmaYD = SigmaYD))
}


### TSHT.VHat
### FUNCTION: Estimates the set of valid and relevant IVs
### INPUT: ITT_Y: a numeric, non-missing vector estimating each instrument's effect on the outcome
###        ITT_D: a numeric, non-missing vector estimating each instrument's effect on the treatment
###        WUMat: an n by pz numeric, non-missing matrix that computes the precision of W 
###               (ex1 t(WUMat) %*% WUMat / n = (W^T W/n)^(-1))
###               (ex2 U = (W^T W/n)^(-1))
###        SigmaSqY: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of Y in the reduced-form model of Y
###        SigmaSqD: a numeric, non-missing, positive scalar value estimating the 
###                  error variance of D in the reduced-form model of D
###        SigmaYD: a numeric, non-missing, scalar value estimating the covariance
###                 of the error terms in the reduced-form models of Y and D
###        tuning, a numeric scalar value tuning parameter for TSHT greater 
###                than 2 (default = 2.01)
### OUTPUT: a list (a) VHat (numeric vector denoting the set of valid and relevant IVs)
###                (b) SHat (numeric vector denoting the set of relevant IVs)
TSHT.VHat <- function(ITT_Y,ITT_D,WUMat,SigmaSqY,SigmaSqD,SigmaYD,tuning = 2.01) {
  # Check ITT_Y and ITT_D
  stopifnot(!missing(ITT_Y),!missing(ITT_D),length(ITT_Y) == length(ITT_D))
  stopifnot(all(!is.na(ITT_Y)),all(!is.na(ITT_D)))
  ITT_Y = as.numeric(ITT_Y); ITT_D = as.numeric(ITT_D)
  
  # Check WUMat
  stopifnot(!missing(WUMat),is.matrix(WUMat), nrow(WUMat) > 1, ncol(WUMat) == length(ITT_Y)) 
  stopifnot(all(!is.na(WUMat)))

  # Check Sigmas 
  stopifnot(!missing(SigmaSqY), is.numeric(SigmaSqY), length(SigmaSqY) == 1, !is.na(SigmaSqY), SigmaSqY > 0)
  stopifnot(!missing(SigmaSqD), is.numeric(SigmaSqD), length(SigmaSqD) == 1, !is.na(SigmaSqD),SigmaSqD > 0)
  stopifnot(!missing(SigmaYD),is.numeric(SigmaYD), length(SigmaYD) == 1,!is.na(SigmaYD))

  # Other Input check
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)

  # Constants
  n = nrow(WUMat); pz = length(ITT_Y)
  
  # First Stage
  SHat = (1:pz)[(abs(ITT_D) >= (sqrt(SigmaSqD * colSums(WUMat^2) /n) * sqrt(tuning*log(pz)/n)))]
  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE

  # Second Stage
  # pi.candidate is the estimated value of pi across different candidates
  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat
  for(j in SHat) {
    beta.j = ITT_Y[j] / ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    sigmasq.j = SigmaSqY + beta.j^2 * SigmaSqD - 2* beta.j * SigmaYD
    PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j * colSums( (WUMat - outer(WUMat[,j]/ITT_D[j], ITT_D))^2)/n) * 
                          sqrt(tuning^2 * log(pz)/n)
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }

  # Voting
  VM = rowSums(VHats.bool)
  VM.m = rownames(VHats.bool)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.bool)[max(VM) == VM] #Plurality winners
  VHat = as.numeric(union(VM.m,VM.p))

  # Error check
  if(length(VHat) == 0){
    warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
    warning("Defaulting to all IVs being valid")
    VHat = 1:pz
  }

  return(list(VHat = VHat,SHat=SHat))
}  
