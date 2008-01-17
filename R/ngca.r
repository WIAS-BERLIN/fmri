ngca <- function(data,L=1000,T=10,m=3,eps=1.5,npca=min(dim(x)[2],dim(x)[1])-1,method="spatial",sweepmean=NULL,keepv=FALSE){
#
#  NGCA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(all(class(data)=="fmridata")) {
   x <- extract.data(data)
   mask <- data$mask
   fmriobj <- TRUE
} else if (class(data)%in%c("matrix","array")){
   x <- data
   mask <- TRUE
   fmriobj <- FALSE
} else {
   warning("data has incompatible class argument")
   return(data)
}
#  
#  x - data matrix  (Nxd)
#
set.seed(1)
xdim <- dim(x)
lxdim <- length(xdim)
d <- xdim[lxdim]
n <- nn <- prod(xdim[1:(lxdim-1)])
dim(x) <- c(n,d)
mask <- as.vector(mask)
if(length(mask)==1) mask <- rep(mask,n)
x <- x[mask,]
n <- sum(mask)
if(is.null(sweepmean)) sweepmean <- switch(method,"spatial"="temporal","temporal"="spatial")
if(is.null(npca)||npca >= min(d,n)) npca <- min(d,n)-1
x <- switch(sweepmean,"none"=x,"global"=x - mean(x),"spatial"=sweep(x,2,apply(x,1,mean)),
            "temporal"=sweep(x,2,apply(x,2,mean)))
if(method=="spatial"){
x <- t(x)
}
svdx <- svd(x,nu=0,nv=npca)
#xvar <- var(x)
#z <- svd(xvar)
cat("Dimension reduced to:",npca,"\n")
y <- t(x%*%svdx$v%*%diag(1/svdx$d[1:npca]))
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
omega <- matrix(rnorm(4*L*npca),npca,L*4)
omega <- sweep(omega,2,sqrt(apply(omega^2,2,sum)),"/")
fz <- .Fortran("fastica",
              as.double(y),
              as.double(omega),
              as.integer(npca),
              as.integer(n),
              as.integer(L),
              as.integer(T),
              double(npca),
              v=double(npca*L*4),
              normv=double(L*4),
              as.double(s),
              DUP=FALSE,
              PACKAGE="fmri")[c("v","normv")]
v <- fz$v
normv <- fz$normv
dim(v) <- c(npca,4*L)
v <- t(v[,normv>eps])
jhat <- prcomp(v)
ihat <- svdx$v%*%diag(svdx$d[1:npca])%*%jhat$rotation[,1:m]
xhat <- x%*%ihat
if(fmriobj){
if(method=="spatial"){
z <- matrix(0,nn,m)
z[mask,]<-ihat
ihat <- array(z,c(xdim[1:3],m))
} else {
z <- matrix(0,nn,m)
z[mask,]<-xhat
xhat <- array(z,c(xdim[1:3],m))
}
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat)
if(keepv) {
z$v <- v
z$normv <- normv
}
class(z) <- "fmringca"
} else {
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat)
if(keepv) {
z$v <- v
z$normv <- normv
}
class(z) <- "ngca"
}
z
}
