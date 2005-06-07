library(awsnew)


create.stimulus <- function(scans=1 ,onsets=c(1) ,length=1, mean=TRUE) {
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

  mystimulus <- function(t) {
    if (t <= 0) {
      0
    } else {
      if (t > scans) {
        0
      } else {
        stimulus[t]
      }
    }
  }

  myhrf <- function(t) {
    mysum <- 0
    for (k in 1:scans) {
      mysum <- mysum + mygamma(k, 6, 12, 0.3, 0.3, 0.35) * mystimulus(t-k)
    }
    mysum
  }

  hrf <- rep(0, scans)
  for (i in 1:scans) {
    hrf[i] <- myhrf(i)
  }

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

  z[,stimuli+1] <- 1

  if (order != 0) {
    for (i in (stimuli+2):(stimuli+order+1)) {
      summe <- rep(0, scans)
      for (k in 1:(i-1)) {
        summe <- summe + z[,k]*(z[,k] %*% (1:scans)^(i-stimuli-1))/(z[,k] %*% z[,k])
      }
      z[,i] <- (1:scans)^(i-stimuli-1) - summe
    }
  }
  
  z
}


create.arcorrection <- function(scans, rho=0) {
  a <- array(0,dim=c(scans,scans))
  
  a[1,1] <- 1

  rho0 <- sqrt(1-rho^2)
  rho1 <- -rho/rho0
  rho2 <- 1/rho0
  
  for (i in 2:scans) {
    a[i,i-1] <- rho1
    a[i,i] <- rho2
  }

  a
}


calculate.lm <- function(ttt,z,ar,vtype="var") {
  zprime <- a %*% z
  svdresult <- svd(zprime)
  u <- svdresult$u
  v <- svdresult$v
  vt <- t(v)
  lambda1 <- diag(1/svdresult$d)
  lambda2 <- diag(1/svdresult$d^2)

  sigma <- v %*% lambda2 %*% vt

  dy <- dim(ttt)
  varianz <- array(0,dim=c(dy[1:3]))
  dim(ttt) <- c(prod(dy[1:3]),dy[4])
  beta <- ttt %*% t(a) %*% u %*% lambda1 %*% vt
  residuals <- ttt - beta %*% t(z)
  b <- rep(1/dim(z)[1],dim(z)[1])
  if (vtype == "var") {
    varianz <- (residuals^2 %*% b)*dim(z)[1]/(dim(z)[1]-dim(z)[2])*sigma[1,1]
  }
  dim(beta)<-c(dy[1:3],dim(z)[2])
  dim(varianz) <- c(dy[1:3])
  dim(residuals) <- dy
  
  result <- list(beta=beta,varianz=varianz,residuals=residuals)
  result
}

perform.aws <- function(beta,varianz,hmax=4,hinit=1,weights=c(1,1,1),qlambda=NULL) {
  varianz[varianz < quantile(varianz,0.25)] <- quantile(varianz,0.25)

  ttthat <- chaws(beta[,,], sigma2=varianz, hmax=hmax, hinit=hinit,
  qlambda=qlambda, wghts=weights)
  tttvar <- ttthat$ni2 / ttthat$ni^2

  z <- list(hat=ttthat$theta, var=tttvar)
  z
}

plot.fmri <- function(signal, mask, anatomic, zlim=0, device="X11", file="plot.png") {
  alim <- range(anatomic)
  zlim <- max(zlim,signal)

  signal[mask==0] <- 0
  anatomic[anatomic<0] <- 0
  
  if (device == "png") {
    png(filename=file, width = 480, height = 480, pointsize=12,
    bg="transparent", res=NA)
  } else {
    X11()
  }
  
  oldpar <- par(mfrow=c(5,6), mar=c(0,0,0,.25), mgp=c(2,1,0))
  
  for (i in 1:dim(anatomic)[3]) {
    image(anatomic[,,i], xaxt="n", yaxt="n", zlim=alim, col=grey(1:255/255))
    if (any(signal[,,i]))
      image(signal[,,i], zlim=c(0,zlim) ,col=c(0,rainbow(512)[350:512]), add=TRUE)
  }
  
  par(oldpar)
  as.numeric(dev.cur())
}
