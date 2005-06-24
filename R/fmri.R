read.ANALYZE <- function(prefix = "", picstart = 0, numbpic = 1) {
  if (require(AnalyzeFMRI)) {

    ttt <- f.read.analyze.volume(paste(prefix, picstart, ".img", sep=""));
    dt <- dim(ttt)
    cat(".")

    if (numbpic > 1) { 
      for (i in (picstart+1):(picstart+numbpic-1)) {
        a <- f.read.analyze.volume(paste(prefix, i, ".img", sep=""))
        if (sum() != 0)
          cat("Error: wrong spatial dimension in picture",i)
        ttt <- c(ttt,a);
        dt[4] <- dt[4] + dim(a)[4]
        cat(".")
      }
    }

    cat("\n")
    dim(ttt) <- c(dt)
    ttt
  } else {
    cat("Error: library AnalyzeFMRI not found\n")
    NA
  }
}

write.ANALYZE <- function(ttt, name = "data", size="int", voxelsize = c(2,2,2)) {
  if (require(AnalyzeFMRI)) {
    f.write.analyze(ttt, file=name, size=size, voxelsize)
  } else {
    cat("Error: library AnalyzeFMRI not found\n")
    NA
  }
}

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
      z[,i] <- (1:scans)^(i-stimuli-1)
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
  variance <- array(0,dim=c(dy[1:3]))
  dim(ttt) <- c(prod(dy[1:3]),dy[4])
  beta <- ttt %*% t(a) %*% u %*% lambda1 %*% vt
  residuals <- ttt - beta %*% t(z)
  b <- rep(1/dim(z)[1],dim(z)[1])
  if (vtype == "var") {
    variance <- (residuals^2 %*% b)*dim(z)[1]/(dim(z)[1]-dim(z)[2])*sigma[1,1]
  }
  dim(beta)<-c(dy[1:3],dim(z)[2])
  dim(variance) <- c(dy[1:3])
  dim(residuals) <- dy
  
  result <- list(beta=beta,var=variance,residuals=residuals)
  result
}

perform.aws <- function(beta,variance,hmax=4,hinit=1,weights=c(1,1,1),vweights=NULL,qlambda=1) {
  require(aws)
  variance[variance < quantile(variance,0.25)] <- quantile(variance,0.25)

  ttthat <- vaws(beta, sigma2=variance, hmax=hmax, hinit=hinit,
                 qlambda=qlambda, qtau=1,wghts=weights,vwghts=vweights)

  z <- list(hat=ttthat$theta, var=ttthat$var)
  z
}

calculate.threshold <- function(var1,var2,alpha=0.95,gamma=0.5) {
  delta <- var1 * qchisq(gamma,1)
  calpha <- qchisq(alpha,1 + delta/var2)
  test <- (var2 + 2 * delta)/(1 + delta/var2) * calpha
  test
}


plot.fmri <- function(signal, mask, anatomic, x=-1, y=-1, z=-1, zlim=0, device="X11", file="plot.png") {
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

  partition <- as.integer(sqrt(dim(anatomic)[3])) + 1
  oldpar <- par(mfrow=c(partition,partition), mar=c(0,0,0,.25), mgp=c(2,1,0))
  
  for (i in 1:dim(anatomic)[3]) {
    image(anatomic[,,i], xaxt="n", yaxt="n", zlim=alim, col=grey(1:255/255))
    if (any(signal[,,i]))
      image(signal[,,i], zlim=c(0,zlim) ,col=c(0,rainbow(512)[350:512]), add=TRUE)
    if (i == z) {
      lines(c(0,1),c(y,y)/dim(anatomic)[2],col=2)
      lines(c(x,x)/dim(anatomic)[1],c(0,1),col=2)
    }     
  }
  
  par(oldpar)
  as.numeric(dev.cur())
}


pvalue <- function(z,i,j,k,rx,ry,rz) {
  rho0 <- 1 - pnorm(z)
  rho1 <- (4 * log(2))^.5 / (2 * pi) * exp(-z*z/2)
  rho2 <- (4 * log(2)) / (2 * pi)^1.5 * exp(-z*z/2) * z
  rho3 <- (4 * log(2))^1.5 / (2 * pi)^2 * exp(-z*z/2) * (z*z -1)
  r0 <- 1
  r1 <- (i-1)*rx + (j-1)*ry +(k-1)*rz
  r2 <- (i-1)*(j-1)*rx*ry + (j-1)*(k-1)*ry*rz +(i-1)*(k-1)*rx*rz
  r3 <- (i-1)*(j-1)*(k-1)*rx*ry*rz

  rho0 * r0 + rho1 * r1 + rho2 * r2 + rho3 * r3
}

pvalueprime <- function(z,i,j,k,rx,ry,rz) {
  rho0prime <- dnorm(z)
  rho1prime <- -(4 * log(2))^.5 / (2 * pi) * exp(-z*z/2) * z
  rho2prime <- -(4 * log(2)) / (2 * pi)^1.5  * exp(-z*z/2) * (z*z -1)
  rho3prime <- -(4 * log(2))^1.5 / (2 * pi)^2 * exp(-z*z/2) * (z*z*z - z)
  r0 <- 1
  r1 <- (i-1)*rx + (j-1)*ry +(k-1)*rz
  r2 <- (i-1)*(j-1)*rx*ry + (j-1)*(k-1)*ry*rz +(i-1)*(k-1)*rx*rz
  r3 <- (i-1)*(j-1)*(k-1)*rx*ry*rz

  rho0prime * r0 + rho1prime * r1 + rho2prime * r2 + rho3prime * r3
}

resel <- function(voxeld, hmax, hv=0.919) hv * voxeld / hmax

threshold <- function(p,i,j,k,rx,ry,rz) {
  x <- 3
  repeat {
    x <- x - (pvalue(x,i,j,k,rx,ry,rz)-p)/pvalueprime(x,i,j,k,rx,ry,rz)
    if (abs(pvalue(x,i,j,k,rx,ry,rz)-p) < 0.000001) break
  }
  x
}

gkernsm <- function(y,h=1) {
  grid <- function(d) {
    d0 <- d%/%2+1
    gd <- seq(0,1,length=d0)
    if (2*d0==d+1) gd <- c(gd,-gd[d0:2]) else gd <- c(gd,-gd[(d0-1):2])
    gd
  }
  dy <- dim(y)
  if (is.null(dy)) dy<-length(y)
  ldy <- length(dy)
  if (length(h)!=ldy) h <- rep(h[1],ldy)
  kern <- switch(ldy,dnorm(grid(dy),0,h/dy),
                 outer(dnorm(grid(dy[1]),0,h[1]/dy[1]),
                       dnorm(grid(dy[2]),0,h[2]/dy[2]),"*"),
                 outer(outer(dnorm(grid(dy[1]),0,h[1]/dy[1]),
                             dnorm(grid(dy[2]),0,h[2]/dy[2]),"*"),
                       dnorm(grid(dy[3]),0,h[3]/dy[3]),"*"))
  kern <- kern/sum(kern)
  kernsq <- sum(kern^2)
  list(gkernsm=convolve(y,kern,conj=TRUE),kernsq=kernsq)
}

correlation <- function(res) {
  dr <- dim(res)
  x <- mean(res[-1,,,]*res[-dr[1],,,]/sqrt(var(res[-1,,,]) * var(res[-dr[1],,,])))
  y <- mean(res[,-1,,]*res[,-dr[2],,]/sqrt(var(res[,-1,,]) * var(res[,-dr[2],,])))
  z <- mean(res[,,-1,]*res[,,-dr[3],]/sqrt(var(res[,,-1,]) * var(res[,,-dr[3],])))
  c(x,y,z)
}

smoothness <- function(cor, dim) {
  h <- 0
  field <- array(rnorm(prod(dim)),dim)
  for (i in 1:100) {
    h <- h + 0.1
    cat(h,"\n")
    z <- gkernsm(field, h)$gkernsm;
    corg <- mean(z[,-1,]*z[,-dim[2],]/sqrt(var(z[,-1,]) * var(z[,-dim[2],])))
    cat(corg)
    if (corg > cor) break
  }
  h
}
