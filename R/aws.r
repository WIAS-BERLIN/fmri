#########################################################################################################################
#
#    R - function  aws3D  for vector-valued  Adaptive Weights Smoothing (AWS) in 3D
#
#    reduced version for qtau==1, heteroskedastic variances only, exact value for Variance reduction
#
#    emaphazises on the propagation-separation approach
#
#    Copyright (C) 2005 Weierstrass-Institut fuer
#                       Angewandte Analysis und Stochastik (WIAS)
#
#    Author:  Joerg Polzehl
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
#  USA.
#

aws3D <- function(y, qlambda=NULL, lkern="Gaussian", skern="Plateau", weighted=TRUE,
                  sigma2=NULL, mask=NULL, hmax=NULL, ladjust=1, u=NULL, wghts=NULL,
                  h0=c(0,0,0), testprop=FALSE, res=NULL){
#
#  qlambda, corrfactor adjusted for case lkern="Gaussian",skern="Plateau" only
#
#  uses sum of weights and correction factor C(h,g) in statistical penalty
#
  #  Auxilary functions
  IQRdiff <- function(y) IQR(diff(y))/1.908

  #
  # first check arguments and initialize
  #
  # test dimension of data (vector of 3D) and define dimension related stuff
  if(is.null(res)) {
      return(warning("Please specify keep=''all'' when calling fmri.lm"))
  }
  # define qlambda, lambda
  if (is.null(qlambda)) qlambda <- switch(skern,.9945,.9988,.9995)
  if (qlambda<.9) warning("Inappropriate value of qlambda")
  if (qlambda<1) {
    lambda <- ladjust*qchisq(qlambda,1)
  } else {
    lambda <- 1e50
  }

  z <- aws::aws3Dmask(y, mask, lambda, hmax, res, sigma2, lkern, skern,
	                    weighted, u, wghts, h0, testprop)

  ## "compress" the residuals
  scale <- max(abs(range(z$residuals)))/32767
  z$residuals <- writeBin(as.integer(z$residuals/scale), raw(), 2)
	z$resscale <- scale
	z$maskOnly <- TRUE
  #   vred accounts for variance reduction with respect to uncorrelated (\check{sigma}^2) data
  class(z) <- "aws.gaussian"
  invisible(z)
}

aws3Dfull <- function(y, qlambda=NULL, lkern="Gaussian", skern="Plateau", weighted=TRUE,
                   sigma2=NULL, mask=NULL, hmax=NULL, ladjust=1, u=NULL, wghts=NULL,
                   testprop=FALSE, res=NULL, resscale=NULL, ddim=NULL) {
#
#  qlambda, corrfactor adjusted for case lkern="Gaussian",skern="Plateau" only
#
#  implements estimation of true variances from residuals (slower) than vaws3D
#
  #  Auxilary functions
  IQRdiff <- function(y) IQR(diff(y))/1.908

  #
  # first check arguments and initialize
  #
  # test dimension of data (vector of 3D) and define dimension related stuff
  if(is.null(res)) {
      return(warning("Please specify keep=''all'' when calling fmri.lm"))
  }
  d <- 3
  dy <- dim(y)
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  if (length(dy)!=3) {
    stop("y has to be 3 dimensional")
  }
  # MAE
  mae <- NULL

  # set the code for the kernel (used in lkern) and set lambda
  lkern <- switch(lkern,
                  Triangle=2,
                  Plateau=1,
                  Gaussian=3,
                  1)
  skern <- switch(skern,
                  Triangle=2,
                  Plateau=1,
                  Exponential=3,
                  2)

  # define qlambda, lambda
  if (is.null(qlambda)) qlambda <- switch(skern,.9945,.9988,.9995)
  if (qlambda<.9) warning("Inappropriate value of qlambda")
  if (qlambda<1) {
    lambda <- qchisq(qlambda,1)
  } else {
    lambda <- 1e50
  }
  if(skern%in%c(1,2)) {
    # to have similar preformance compared to skern="Exp"
    lambda <- 4/3*lambda
    if(skern==1) spmin <- .3
    spmax <- 1
  } else {
      spmin <- 0
      spmax <- 4
  }
  # set hinit and hincr if not provided
  hinit <- 1

  # define hmax
  if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points

  # re-define bandwidth for Gaussian lkern!!!!
  if (lkern==3) {
    # assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- fwhm2bw(hmax)*4
    hinit <- min(hinit,hmax)
  }
  if (qlambda == 1) hinit <- hmax

  # define hincr
  hincr <- 1.25
  hincr <- hincr^(1/3)

  if (length(dim(sigma2))!=3) stop("insufficient variance information for full analysis")
  if(is.null(mask)) mask <- array(TRUE,dy)
  mask[sigma2>=1e16] <- FALSE
#  in these points sigma2 probably contains NA's
  # deal with homoskedastic Gaussian case by extending sigma2
  if (length(sigma2)!=n) stop("sigma2 does not have length 1 or same length as y")
  dim(sigma2) <- dy
  sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes


  # Initialize  list for bi and theta
  if (is.null(wghts)) wghts <- c(1,1,1)
  hinit <- hinit/wghts[1]
  hmax <- hmax/wghts[1]
  wghts <- (wghts[2:3]/wghts[1])
  tobj <- list(bi= sigma2)
  theta <- y
  maxvol <- aws::getvofh(hmax,lkern,wghts)
  kstar <- as.integer(log(maxvol)/log(1.25))
  steps <- kstar+1

  if(qlambda < 1) k <- 1 else k <- kstar
  hakt <- hinit
  hakt0 <- hinit
  lambda0 <- lambda
  if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
  total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
  if(testprop) {
    if(is.null(u)) u <- 0
  }
  propagation <- NULL
##
##  determine number of cores to use
##
  mc.cores <- setCores(,reprt=FALSE)
  # run single steps to display intermediate results
  residuals <- readBin(res,"integer",prod(ddim),2)*resscale
  dim(residuals) <- c(ddim[4],ddim[1:3])
  vadjust <- apply(residuals,2:4,var)*sigma2
  while (k<=kstar) {
      hakt0 <- aws::gethani(1,3,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- aws::gethani(1,3,lkern,1.25^k,wghts,1e-4)
      hakt.oscale <- if(lkern==3) bw2fwhm(hakt/4) else hakt
      cat("step",k,"bandwidth",signif(hakt.oscale,3)," ")
    dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:3]
    hakt0 <- hakt
    theta0 <- theta
    bi0 <- tobj$bi/(1+log(tobj$bi/sigma2))
    #
    #   need these values to compute variances after the last iteration
    #
    tobj <- .Fortran(C_chawsv,as.double(y),
                     as.double(residuals),
                     as.double(sigma2),
                     as.integer(!mask),
                     as.integer(weighted),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(ddim[4]),
                     hakt=as.double(hakt),
                     as.double(lambda0),
                     as.double(theta0),
#                     as.integer(mc.cores),
                     bi=as.double(bi0),
                     resnew=double(prod(ddim)),
                     thnew=double(n1*n2*n3),
                     as.integer(lkern),
                     as.integer(skern),
                     as.double(spmin),
                     as.double(spmax),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(ddim[4]*mc.cores))[c("bi","thnew","hakt","resnew")]
    gc()
    theta <- array(tobj$thnew,dy)
    dim(tobj$bi) <- ddim[1:3]
    tobj$bi <- tobj$bi*vadjust
#  correcting for missing term in variance estimate from residuals, e.g. effects of design matrix
#  now contains 1/var(thetahat)
    if(testprop) {
      pobj <- .Fortran(C_chawsv,as.double(y),
                        as.double(residuals),
                        as.double(sigma2),
                        as.integer(!mask),
                        as.integer(weighted),
                        as.integer(n1),
                        as.integer(n2),
                        as.integer(n3),
                        as.integer(ddim[4]),
                        hakt=as.double(hakt),
                        as.double(1e50),
                        as.double(theta0),
  #                      as.integer(mc.cores),
                        bi=as.double(bi0),
                        resnew=double(prod(ddim)),
                        thnew=double(n1*n2*n3),
                        as.integer(lkern),
                        as.integer(skern),
                        as.double(spmin),
                        as.double(spmax),
                        double(prod(dlw)),
                        as.double(wghts),
                        double(ddim[4]*mc.cores))[c("bi","thnew","resnew")]
      ptheta <- array(pobj$thnew,dy)
      rm(pobj)
      gc()
      propagation <- c(propagation,sum(abs(theta-ptheta))/sum(abs(ptheta-u)))
      cat("Propagation with alpha=",max(propagation),"\n")
      cat("alpha values:","\n")
      print(rbind(hakt.oscale,signif(propagation[-1],3)))
    }
      if (!is.null(u)) {
      cat("bandwidth: ",signif(hakt.oscale,3),"eta==1",sum(tobj$eta==1),"   MSE: ",
          signif(mean((theta-u)^2),3),"   MAE: ",signif(mean(abs(theta-u)),3)," mean(bi)=",signif(mean(tobj$bi),3),"\n")
      mae <- c(mae,signif(mean(abs(theta-u)),3))
    } else if (max(total) >0) {
      cat(signif(total[k],2)*100,"%                 \r",sep="")
     }
    k <- k+1

    lambda0 <- if(k<=length(ladjust)) ladjust[k]*lambda else ladjust*lambda#/quot
    gc()
  }
  cat("fmri.smooth: estimate correlations","\n")
  lags <- pmin(c(5,5,3),ddim[1:3]-1)
  scorr <- aws::residualSpatialCorr(tobj$resnew,mask,lags)
    vartheta <- 1/tobj$bi
  vred <- vartheta*sigma2
  ## "compress" the residuals
  scale <- max(abs(range(tobj$resnew)))/32767
  residuals <- writeBin(as.integer(tobj$resnew/scale), raw(), 2)
  z <- list(theta=theta,ni=tobj$bi,var=vartheta,vred=vred,
            vred0=median(vred[mask]),y=y,
            hmax=tobj$hakt*switch(lkern,1,1,bw2fwhm(1/4)),mae=mae,
            alpha=propagation,scorr=scorr,scale=scale,res=residuals,mask=mask)
  #   vred accounts for variance reduction with respect to uncorrelated (\check{sigma}^2) data
  class(z) <- "aws.gaussian"
  invisible(z)
}

smooth3D <- function(y,lkern="Gaussian",weighted=FALSE,sigma2=NULL,mask=NULL,hmax=NULL,
                     wghts=NULL) {
  #
  # first check arguments and initialize
  #
  # test dimension of data (vector of 3D) and define dimension related stuff
  d <- 3
  dy <- dim(y)
  n1 <- dy[1]
  n2 <- dy[2]
  n3 <- dy[3]
  n <- n1*n2*n3
  if (length(dy)==d) {
    dim(y) <- dy <- c(dy,1)
  } else if (length(dy)!=d+1) {
    stop("y has to be 3 or 4 dimensional")
  }
  dv <- dim(y)[d+1]
    if(is.null(sigma2)) {
       weighted <- FALSE
    } else {
      if(length(sigma2)!=n) weighted <- FALSE
      sigma2 <- 1/sigma2
    }
  if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points

  # re-define bandwidth for Gaussian lkern!!!!
  lkern <- switch(lkern,
                  Triangle=2,
                  Plateau=1,
                  Gaussian=3,
                  1)
  if (lkern==3) {
    # assume  hmax was given in  FWHM  units (Gaussian kernel will be truncated at 4)
    hmax <- fwhm2bw(hmax)*4
  }
  if (is.null(wghts)) wghts <- c(1,1,1)
  if(is.null(mask)) mask <- array(TRUE,dy[1:3])
  hmax <- hmax/wghts[1]
  wghts <- (wghts[2:3]/wghts[1])
    dlw <- (2*trunc(hmax/c(1,wghts))+1)[1:d]
    ysmooth <- .Fortran(C_smooth3d,
                     as.double(y),
                     as.double(sigma2),
                     as.integer(!mask),
                     as.integer(weighted),
                     as.integer(n1),
                     as.integer(n2),
                     as.integer(n3),
                     as.integer(dv),
                     hakt=as.double(hmax),
                     thnew=double(n1*n2*n3*dv),
                     as.integer(lkern),
                     double(prod(dlw)),
                     as.double(wghts),
                     double(dv))$thnew
array(ysmooth,dy)
}
