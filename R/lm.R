create.stimulus <- function(scans=1 ,onsets=c(1) ,length=1, rt=3, mean=TRUE) {
  numberofonsets <- length(onsets)
  stimulus <- rep(0, scans)
  
  for (i in 1:numberofonsets) {
    onset <- onsets[i]
    for (j in onset:(onset+length-1)) {
      stimulus[j] <- 1
    }
  }

  mygamma <- function(x, a1, a2, b1, b2, c) {
    d1 <- a1 * b1
    d2 <- a2 * b2
    c1 <- ( x/d1 )^a1
    c2 <- c * ( x/d2 )^a2
    res <- c1 * exp(-(x-d1)/b1) - c2 * exp(-(x-d2)/b2)
    res
  }

  hrf <- convolve(stimulus,mygamma(scans:1, 6, 12, 0.9/rt, 0.9/rt, 0.35))
  dim(hrf) <- c(scans,1)
  
  if (mean) {
    hrf - mean(hrf)
  } else {
    hrf
  }
}



create.designmatrix <- function(hrf, order=0) {
  stimuli <- dim(hrf)[2]
  scans <- dim(hrf)[1]

  z <- matrix(0, scans, stimuli+order+1)

  for (i in 1:stimuli) {
    z[,i] <- hrf[,i]
  }

  ortho <- matrix(0, stimuli, stimuli)
  for (i in 1:stimuli) {
    for (j in 1:stimuli) {
      ortho[i,j] <- z[,i]%*%z[,j]
    }
  }

  z[,stimuli+1] <- 1

  if (order != 0) {
    for (i in (stimuli+2):(stimuli+order+1)) {
      z[,i] <- (1:scans)^(i-stimuli-1)
      hz <- numeric(stimuli)
      for (j in 1:stimuli) {
        hz[j] <- z[,j]%*%z[,i]
      }
      tmp <- lm(-hz~ortho-1)
      z[,i] <- z[,i] + as.vector(hrf %*% as.vector(tmp$coeff))
    }
  }
  
  z
}



create.arcorrection <- function(scans, rho=0) {
  rho0 <- 1/sqrt(1-rho^2)
  a <- numeric(scans*scans)
  
  a[1] <- 1
  ind <- (2:scans) *(scans+1) - 2*scans
  a[ind] <- -rho*rho0
  a[ind+scans] <- rho0
  dim(a) <- c(scans,scans)
  
  a
}



calculate.lm <- function(ttt,z,actype="smooth",hmax=3.52,vtype="var",step=0.01,contrast=c(1)) {
  require(aws)
  cat("calculate.lm: entering function with:",actype, hmax, vtype, "\n")

  # first consider contrast vector! NO test whether it is real contrast!!
  if (length(contrast) <= dim(z)[2]) contrast <- c(contrast,rep(0,dim(z)[2]-length(contrast)))
  length(contrast) <- dim(z)[2]

  # first get the SVD for the design matrix
  svdresult <- svd(z)
  u <- svdresult$u
  v <- svdresult$v
  vt <- t(v)
  lambda1 <- diag(1/svdresult$d)
  lambda2 <- diag(1/svdresult$d^2)
  xtx <- v %*% lambda2 %*% vt
  # now we have z = u lambda1^(-1) vt

  # define some variables and make ttt a matrix
  dy <- dim(ttt)
  voxelcount <- prod(dy[1:3])
  dim(ttt) <- c(prod(dy[1:3]),dy[4])
  arfactor <- rep(0,length=prod(dy[1:3]))
  variance <- rep(0,length=prod(hrf,dy[1:3]))

  # calculate matrix R for bias correction in correlation coefficient
  # estimate
  R <- diag(1,dy[4]) - u %*% t(u)
  m00 <- dy[4] - dim(z)[2]
  m01 <- 0
  for (k in 1:dy[4]) {m01 <- m01 + sum(R[k,-1]*R[k,-dy[4]])}
  m10 <- m01
  m11 <- 0
  for (k in 1:(dy[4]-1)) {m11 <- m11 + sum(R[k+1,-dy[4]]*R[k,-1] + R[k+1,-1]*R[k,-dy[4]]) }
  Minv <- matrix(c(m11,-m01,-m10,m00),2,2)/(m00*m11-m01*m10) # inverse
  
  # calculate the paramters and residuals for all voxels
  beta <- ttt %*% u %*% lambda1 %*% vt
  residuals <- ttt - beta %*% t(z)

  # actype == "smooth" ... calc AC, smooth AC, calc prewhitened model
  # actype == "accalc" ... calc AC, calc prewhitened model
  # actype == "ac"     ... calc AC only
  # "accalc" is actually a special case of "smooth" (hmax=1), but we
  # leave it here for clearer function interface, so parameter set is
  # only needed for "smooth"
  if ((actype == "smooth") || (actype == "accalc") || (actype == "ac")) { 
    progress = 0
    cat("calculate.lm: calculating AR(1) model\n")
    for (i in 1:voxelcount) {
      if (i > progress/100*voxelcount) {
        cat(progress,"% . ",sep="")
        progress = progress + 10
      }
      # calculate the Koeff of ACR(1) time series model
      a0 <- residuals[i,] %*% residuals[i,]
      a1 <- residuals[i,-1] %*% residuals[i,-dim(z)[1]]
      an <- Minv %*% c(a0,a1)
      if (an[1] != 0) {
        arfactor[i] <- an[2]/an[1]
      } else {
        arfactor[i] <- 0
      }

      ### this method does have a bias!
#      if (sum(abs(residuals[i,]))) {
#        arfactor[i] <- residuals[i,2:dim(z)[1]] %*% residuals[i,1:(dim(z)[1]-1)] / residuals[i,] %*% residuals[i,]
#      } else {
#        arfactor[i] <- 0 # avoid NaN
#      }
      ### leave it for historical reasons
      
    }
    cat("\n")
    cat("calculate.lm: finished\n")
    
    if (actype == "smooth") {
      cat("calculate.lm: smoothing with (hmax):",hmax,"\n")
      dim(arfactor) <- dy[1:3]
      # now smooth (if actype is such) with AWS
      hinit <- 1
      arfactor <- gkernsm(arfactor,rep(hmax,3)*0.42445)$gkernsm
      dim(arfactor) <- voxelcount
      cat("calculate.lm: finished\n")
    }

    if ((actype == "smooth") || (actype == "accalc")) {
      progress = 0
      cat("calculate.lm: re-calculating linear model with prewithened data\n")
      # re- calculated the linear model with prewithened data
      # NOTE: sort arfactors and bin them! this combines calculation in
      # NOTE: voxels with similar arfactor, see Worsley
      variancepart <- rep(0,length=prod(dy[1:3]))
      arlist <- seq(range(arfactor)[1]-step/2,range(arfactor)[2]+step/2,step)
      for (i in 1:(length(arlist)-1)) {
        if (i > progress/100*length(arlist)) {
          cat(progress,"% . ",sep="")
          progress = progress + 10
        }
        indar <- as.logical((arfactor > arlist[i]) * (arfactor <=
                                                      arlist[i+1]))
        if(sum(indar)>0){
          a <- create.arcorrection(dy[4],mean(arlist[i:(i+1)])) # create prewhitening matrix
          zprime <- a %*% z
          svdresult <- svd(zprime) # calc SVD of prewhitened design
          v <- svdresult$v
          vt <- t(v)
          xtx <- v %*% diag(1/svdresult$d^2) %*% vt # xtx * <estimate of varince> of prewhitened noise is variance of parameter estimate
          tttprime <- ttt[indar,] %*% t(a)
          beta[indar,] <- tttprime %*% svdresult$u %*% diag(1/svdresult$d) %*% vt # estimate parameter
          residuals[indar,] <- tttprime - beta[indar,] %*% t(zprime) # calculate residuals
          variancepart[indar] <- t(contrast) %*% xtx %*% contrast # variance estimate
        }
      }
      b <- rep(1/dy[4],length=dy[4])
      variance <- ((residuals^2 %*% b) * dim(z)[1] / (dim(z)[1]-dim(z)[2])) * variancepart
      cat("\n")
      cat("calculate.lm: finished\n")
    } else {
      b <- rep(1/dy[4],length=dy[4])
      cxtx <- t(contrast) %*% xtx %*% contrast
      variance <- (residuals^2 %*% b) * dy[4] / (dy[4]-dim(z)[2]) * cxtx
    }
  } else { # actype == "noac"
    # estimate variance, add more paramters if needed
    b <- rep(1/dy[4],length=dy[4])
    cxtx <- t(contrast) %*% xtx %*% contrast
    variance <- (residuals^2 %*% b) * dy[4] / (dy[4]-dim(z)[2]) * cxtx
  }

  cbeta <- beta %*% contrast
  # re-arrange dimensions
  dim(beta) <- c(dy[1:3],dim(z)[2])
  dim(cbeta) <- dy[1:3]
  dim(variance) <- dy[1:3]
  dim(arfactor) <- dy[1:3]
  dim(residuals) <- dy

  cat("calculate.lm: exiting function\n")
  
  list(beta = beta, cbeta = cbeta, var = variance, res = residuals, arfactor = arfactor)
}

degrees <- function(z, fwhmcorr, bw) {
  x <- (z %*% (t(z) %*% z))[,1]
  tau1 <- sum(x[-1] %*% x[-length(x)])/sum(x^2)
  f <- (1+2*(fwhmcorr^2)/(bw^2))^(-1.5)
  (dim(z)[1]-4)/(1+2*f*tau1^2)
}
