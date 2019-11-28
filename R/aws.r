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
#  implements estimation of true variances from residuals (slower) than aws3D
#
  if(is.null(res)) {
      return(warning("Please specify keep=''all'' when calling fmri.lm"))
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

	z <- aws::aws3Dmaskfull(y, mask, lambda, hmax, res, sigma2, lkern, skern,
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
