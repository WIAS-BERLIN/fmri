ngca <- function(x,L=1000,T=10,m=3,eps=1.5){
#
#  NGCA algorithm 
#
#  
#  x - data matrix  (Nxd)
#
xdim <- dim(x)
d <- xdim[2]
n <- xdim[1]
xmean <- apply(x,2,mean)
xvar <- var(x)
y <- sweep(x,2,xmean)
z <- svd(xvar)
y <- y%*%z$u%*%diag(z$d^(-.5))
y <- t(y)
#
#  thats the standardized version of x
#
s <- matrix(0,L,4)
s[,1] <- seq(.5,5,length=L) 
s[,2] <- seq(5/L,5,length=L) 
s[,3] <- seq(4/L,4,length=L) 
s[,4] <- seq(0,4,length=L) 
#
#   now fast ICA
#
omega <- matrix(rnorm(4*L*d),d,L*4)
omega <- sweep(omega,2,sqrt(apply(omega^2,2,sum)),"/")
fz <- .Fortran("fastica",
              as.double(y),
              as.double(omega),
              as.integer(d),
              as.integer(n),
              as.integer(L),
              as.integer(T),
              double(d),
              v=double(d*L*4),
              normv=double(L*4),
              as.double(s),
              DUP=FALSE,
              PACKAGE="fmri")[c("v","normv")]
dim(fz$v) <- c(d,4*L)
fz$v <- t(fz$v[,fz$normv>eps])
jhat <- prcomp(fz$v)
ihat <- z$u%*%jhat$rotation[,1:m]
list(ihat=ihat,sdev=jhat$sdev[1:m])
}
