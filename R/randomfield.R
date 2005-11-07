pvalue <- function(z,i,j,k,rx,ry,rz,type="norm",df=4,dim=3) {
  rho0 <- 0
  rho1 <- 0
  rho2 <- 0
  rho3 <- 0
  
  if (type == "norm") {
    rho0 <- 1 - pnorm(z)
    rho1 <- (4 * log(2))^.5 / (2 * pi) * exp(-z*z/2)
    rho2 <- (4 * log(2)) / (2 * pi)^1.5 * exp(-z*z/2) * z
    rho3 <- (4 * log(2))^1.5 / (2 * pi)^2 * exp(-z*z/2) * (z*z -1)
  }
  if (type == "t") {
    rho0 <- 1 - pt(z,df)
    rho1 <- (4 * log(2))^.5 / (2 * pi) * (1+z*z/df)^(-0.5*(df-1))
    rho2 <- (4 * log(2)) / (2 * pi)^1.5 * gamma((df+1)/2)/gamma(df/2)/(df/2)^0.5 * (1+z*z/df)^(-0.5*(df-1)) * z
    rho3 <- (4 * log(2))^1.5 / (2 * pi)^2 * (1+z*z/df)^(-0.5*(df-1)) * ((df-1)/df*z*z-1)
  }
  if (type == "chisq") {
    rho0 <- 1 - pchisq(z,df)
    rho1 <- (4 * log(2))^.5 / (2 * pi)^.5 * chihelp(z,df,1)
    rho2 <- (4 * log(2)) / (2 * pi) * chihelp(z,df,2) * (z-df+1)
    rho3 <- (4 * log(2))^1.5 / (2 * pi)^1.5 * chihelp(z,df,3) * (z*z - (2*df-1)*z + (df-1) * (df-2))    
  }
  
  r0 <- 1
  r1 <- (i-1)*rx + (j-1)*ry +(k-1)*rz
  r2 <- (i-1)*(j-1)*rx*ry + (j-1)*(k-1)*ry*rz +(i-1)*(k-1)*rx*rz
  r3 <- (i-1)*(j-1)*(k-1)*rx*ry*rz
  
  if (dim == 2) {
    rho0 * r0 + rho1 * r1 + rho2 * r2
  } else {
    rho0 * r0 + rho1 * r1 + rho2 * r2 + rho3 * r3
  }
}



chihelp <- function(z,df,d) {
  z^((df-d)/2) * exp(-z/2) / 2^((df-2)/2) / gamma(df/2)
}



resel <- function(voxeld, hmax, hv=0.95) hv * voxeld / hmax # hv=0.919 for larger bandwidths than typically fmri



#threshold <- function(p,i,j,k,rx,ry,rz,type="norm",df=4,dim=3) {
#  f <- function(x,par){
#    abs(pvalue(x,par$i,par$j,par$k,par$rx,par$ry,par$rz,par$type,par$df,par$dim)-par$p)
#  }
#  optimize(f=f,lower=2,upper=10,tol=1e-6,par=list(p=p,i=i,j=j,k=k,rx=rx,ry=ry,rz=rz,type=type,df=df,dim=dim))$minimum
#}

threshold <- function(p,i,j,k,rx,ry,rz,type="norm",df=4,dim=3,step=.01) {
  n <- length(rx)
  thr <- numeric(n)
  fixed <- logical(n)
  x <- 3
  while (any(!fixed)) {  
    pv <- pvalue(x,i,j,k,rx,ry,rz,type,df,dim)
    ind <- (1:n)[!fixed][pv[!fixed]<p]
    thr[ind] <- x
    fixed[ind] <- TRUE
    x <- x+step
  }
  thr
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
  kern <- switch(ldy,dnorm(grid(dy),0,2*h/dy),
                 outer(dnorm(grid(dy[1]),0,2*h[1]/dy[1]),
                       dnorm(grid(dy[2]),0,2*h[2]/dy[2]),"*"),
                 outer(outer(dnorm(grid(dy[1]),0,2*h[1]/dy[1]),
                             dnorm(grid(dy[2]),0,2*h[2]/dy[2]),"*"),
                       dnorm(grid(dy[3]),0,2*h[3]/dy[3]),"*"))
  kern <- kern/sum(kern)
  kernsq <- sum(kern^2)
  list(gkernsm=convolve(y,kern,conj=TRUE),kernsq=kernsq)
}

get.corr.gauss <- function(h) {
  #
  #   Calculates the correlation of 
  #
  #   colored noise that was produced by smoothing with "gaussian" kernel and bandwidth h
  #
  #   Result does not depend on d for "Gaussian" kernel !!
  #
  h <- h/2.3548
  ih <- trunc(4*h+1)
  dx <- 2*ih+1
  penl <- dnorm(((-ih):ih)/h)
  sum(penl[-1]*penl[-dx])/sum(penl^2)
}

get.bw.gauss <- function(corr, step = 1.01) {
  # get the   bandwidth for lkern corresponding to a given correlation
  h <- .1
  z <- 0
  # 
  #  keep it simple result does not depend on d
  #
  while (z<corr) {
    h <- h*step
    z <- get.corr.gauss(h)
  }
  h
}

correlation <- function(res,mask) {
  meanpos <- function(a) mean(a[a!=0])
  varpos <- function(a) var(a[a!=0])

  dr <- dim(res)
  if ( (length(dim(mask)) == length(dr))
      && (sum(dim(mask)[1:3] != dr[1:3]) == 0)
      && (dim(mask)[4] == 1)) {
    mask <- rep(mask,dr[4])
    dim(mask) <- dr
    vrm <- varpos(res*mask)
    x <- meanpos(res[-1,,,]*res[-dr[1],,,]*mask[-1,,,])/vrm
    y <- meanpos(res[,-1,,]*res[,-dr[2],,]*mask[,-1,,])/vrm
    z <- meanpos(res[,,-1,]*res[,,-dr[3],]*mask[,,-1,])/vrm
    c(x,y,z)
  } else {
    cat("Error: dimension of mask and residui matrices do not match\n")    
  }
}

correlation2 <- function(res,mask) {
  meanpos <- function(a) mean(a[a!=0])
  varpos <- function(a) var(a[a!=0])

  dr <- dim(res)
  if ( (length(dim(mask)) == length(dr))
      && (sum(dim(mask)[1:3] != dr[1:3]) == 0)
      && (dim(mask)[4] == 1)) {
    mask <- rep(mask,dr[4])
    dim(mask) <- dr
    vrm <- varpos(res*mask)
    x <- meanpos(res[-c(1,2),,,]*res[-dr[1]+c(0,1),,,]*mask[-c(1,dr[1]),,,])/vrm
    y <- meanpos(res[,-c(1,2),,]*res[,-dr[2]+c(0,1),,]*mask[,-c(1,dr[2]),,])/vrm
    z <- meanpos(res[,,-c(1,2),]*res[,,-dr[3]+c(0,1),]*mask[,,-c(1,dr[3]),])/vrm
    c(x,y,z)
  } else {
    cat("Error: dimension of mask and residui matrices do not match\n")    
  }
}


correlation3 <- function(res,mask) {
  meanpos <- function(a) mean(a[a!=0])
  varpos <- function(a) var(a[a!=0])

  dr <- dim(res)
  if ( (length(dim(mask)) == length(dr))
      && (sum(dim(mask)[1:3] != dr[1:3]) == 0)
      && (dim(mask)[4] == 1)) {
    mask <- rep(mask,dr[4])
    dim(mask) <- dr
    vrm <- varpos(res*mask)
    x <- meanpos(res[-c(1,2,3),,,]*res[-dr[1]+c(0,1,2),,,]*mask[-c(1,2,dr[1]),,,])/vrm
    y <- meanpos(res[,-c(1,2,3),,]*res[,-dr[2]+c(0,1,2),,]*mask[,-c(1,2,dr[2]),,])/vrm
    z <- meanpos(res[,,-c(1,2,3),]*res[,,-dr[3]+c(0,1,2),]*mask[,,-c(1,2,dr[3]),])/vrm
    c(x,y,z)
  } else {
    cat("Error: dimension of mask and residui matrices do not match\n")    
  }
}



smoothness <- function(cor) {
  for (i in 1:dim(cor2bw)[1]) {
    if (cor2bw[i,2]>cor) break
  }
  cor2bw[i,1]
}



bandwidth <- function(res,mask) { # second argument !!!!
  require(aws)
  cxyz <- correlation(res,mask)
  bwx <- geth.gauss(cxyz[1])
  bwy <- geth.gauss(cxyz[2])
  bwz <- geth.gauss(cxyz[3])
  # 0.42445 is factor for changing units from FWHM to grid (gkernsm
  # works with grid)
  meingauss<-gkernsm(array(rnorm(prod(dim(res)[1:3])),dim=dim(res)[1:3]),c(bwx,bwy,bwz)*0.42445)
  list(bwcorr=1/meingauss$kernsq,bw=c(bwx,bwy,bwz),corr=cxyz)
}


#
# the following function is possibly obsolete since we do not
# implement an algorithm using the derivative 
#
pvalueprime <- function(z,i,j,k,rx,ry,rz,type="norm",df=4,dim=3) {
  rho0prime <- 0
  rho1prime <- 0
  rho2prime <- 0
  rho3prime <- 0
  
  if (type == "norm") {
    rho0prime <- dnorm(z)
    rho1prime <- -(4 * log(2))^.5 / (2 * pi) * exp(-z*z/2) * z
    rho2prime <- -(4 * log(2)) / (2 * pi)^1.5  * exp(-z*z/2) * (z*z -1)
    rho3prime <- -(4 * log(2))^1.5 / (2 * pi)^2 * exp(-z*z/2) * (z*z*z - z)
  }
  if (type == "t") {
    rho0 <- 1 - pt(z,df)
    rho1 <- 0
    rho2 <- 0
    rho3 <- 0
  }
  if (type == "chisq") {
    rho0 <- 1 - pchisq(z,df)
    rho1 <- 0
    rho2 <- 0
    rho3 <- 0    
  }
    
  r0 <- 1
  r1 <- (i-1)*rx + (j-1)*ry +(k-1)*rz
  r2 <- (i-1)*(j-1)*rx*ry + (j-1)*(k-1)*ry*rz +(i-1)*(k-1)*rx*rz
  r3 <- (i-1)*(j-1)*(k-1)*rx*ry*rz

  if (dim == 2) {
    rho0prime * r0 + rho1prime * r1 + rho2prime * r2
  } else {
    rho0prime * r0 + rho1prime * r1 + rho2prime * r2 + rho3prime * r3
  }
}

# simply local thresholding not used!
calculate.threshold <- function(var1,var2,alpha=0.95,gamma=0.5) {
  delta <- var1 * qchisq(gamma,1)
  calpha <- qchisq(alpha,1 + delta/var2)
  test <- (var2 + 2 * delta)/(1 + delta/var2) * calpha
  test
}

#
# end of obsolete functions
#

