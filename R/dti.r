create.designmatrix.dti <- function(bvec, bvalue=1) {
  dimension <- dim(bvec)[2] # should be 3
  if (dimension != 3) {
    warning("Error: gradient vectors do not have length 3")
    return(invisible(NULL))
  }
  directions <- dim(bvec)[1] # number of measured directions

  z <- matrix(0, directions, 6)

  for (d in 1:directions) {
    z[d,1] <- bvec[d,1]^2
    z[d,2] <- bvec[d,2]^2
    z[d,3] <- bvec[d,3]^2
    z[d,4] <- 2*bvec[d,1]*bvec[d,2]
    z[d,5] <- 2*bvec[d,1]*bvec[d,3]
    z[d,6] <- 2*bvec[d,2]*bvec[d,3]
    z[d,] <- bvalue*z[d,]
  }
  
  z
}

calculate.lm.dti <- function(ttt,z,res=FALSE) {
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
  
  # calculate the paramters and residuals for all voxels
  beta <- ttt %*% u %*% lambda1 %*% vt
  residuals <- ttt - beta %*% t(z)
  b <- rep(1/dy[4],length=dy[4])
  variance <- ((residuals^2 %*% b) * dim(z)[1] / (dim(z)[1]-dim(z)[2]))
  dim(variance) <- c(dy[1:3])
  dim(beta) <- c(dy[1:3],dim(z)[2])
  dim(residuals) <- c(dy[1:3],dim(z)[1])

  if (res) {
    list(dt=beta,xtx=xtx,variance=variance,residuals=residuals)
  } else {
    list(dt=beta,xtx=xtx,variance=variance)    
  }
}
  
determine.eigenvalue <- function(diff) {

  eigen <- array(0,dim=c(dim(diff)[1],dim(diff)[2],dim(diff)[3],3))
  for (i in 1:dim(diff)[1]) {
    cat(".")
    for (j in 1:dim(diff)[2]) {
      for (k in 1:dim(diff)[3]) {
        eigen[i,j,k,] <- eigen(matrix(diff[i,j,k,c(1,4,5,4,2,6,5,6,3)],c(3,3)),symmetric=TRUE,only.values=TRUE)$values
      }
    }
  }
# ediff<-apply(array(diff[,,,c(1,2,3,2,4,5,3,5,6)],c(dim(diff)[-4],3,3)),1:3,eigen)  
#  gives a list of dimension dim(diff)  containing the eigenvalues and eigenvectors
  list(eigen=eigen)
}

determine2.eigenvalue <- function(diff,only.values=TRUE) {
  n<-3
  ddiff<-prod(dim(diff)[-4])
  z<-.Fortran("eigen3", as.integer(ddiff), as.double(diff), evectors = double(ddiff*9), evalues = double(ddiff*3),
            !only.values, ierr = integer(ddiff), package="fmri")[c("evalues","evectors","ierr")]
  list(evalues=array(t(matrix(z$evalues,3,ddiff)),c(dim(diff)[-4],3)),evectors=array(t(matrix(z$evectors,9,ddiff)),c(dim(diff)[-4],3,3)),ierr=z$ierr)
}

anisotropy <- function(eigen) {
  
  fa <-
    array(0,dim=c(dim(eigen)[1],dim(eigen)[2],dim(eigen)[3]))
  ra <-
    array(0,dim=c(dim(eigen)[1],dim(eigen)[2],dim(eigen)[3]))
  trc <-
    array(0,dim=c(dim(eigen)[1],dim(eigen)[2],dim(eigen)[3]))
  
  for (i in 1:dim(eigen)[1]) {
    cat(".")
    for (j in 1:dim(eigen)[2]) {
      for (k in 1:dim(eigen)[3]) {
        trc[i,j,k] <- mean(eigen[i,j,k,])
        fa[i,j,k] <- sqrt(3*sum((eigen[i,j,k,]-trc[i,j,k])^2))/sqrt(2*(sum(eigen[i,j,k,]^2)))
        ra[i,j,k] <- sqrt(sum((eigen[i,j,k,]-trc[i,j,k])^2))/sqrt(3*trc[i,j,k])
      }
    }
  }

  list(fa=fa,ra=ra, trace=trc)
}
