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

segm3D <- function(y,weighted=TRUE,
                   sigma2=NULL,mask=NULL,hinit=NULL,hmax=NULL,
                   ladjust=1,graph=FALSE,wghts=NULL,
                   df=100,h0=c(0,0,0),res=NULL, resscale=NULL, 
                   ddim=NULL,delta=0,fov=NULL,alpha=.05) {
#
#
#  Auxilary functions
   IQRdiff <- function(y) IQR(diff(y))/1.908
   getbeta <- function(df,h0){
#    df in (40 : 100)
#    h0 in (0 : 1.5) (in FWHM)
#    see file sim_fmri_kritval.r in R/segmentation/fmrikrv/
      h <- fwhm2bw(sqrt(h0[1]*h0[2]))
      1.461+2150/df^2+0.3*h
}
   getkrval <- function(df,h0,ladj,n,kstar,alpha){
#
#    this delivers an upper bound for kritical values over a wide range of parameters
#    covering the typical situations in fMRI
#    n in (32^3 : 64^2*32)
#    df in (40 : 100)
#    h0 in (0 : 1.5) (in FWHM)
#    ladj in (1 : 1.6) 
#    kstar in (20:29)  corrsponding to maximal bandwidths  3 - 6 
#    internal adjustment by 0.1*kstar 
#    see file sim_fmri_kritval.r in R/segmentation/fmrikrv/

      h <- fwhm2bw(sqrt(h0[1]*h0[2]))
      explvar <- c(1,              1/df,               n,                        ladj,
                   log(ladj),      log(h+1),           kstar,                    h,
                   1/df^2,         1/df*log(ladj),     n*log(h+1),               n*log(ladj),
                   ladj*h,         log(ladj)*log(h+1), log(ladj)*h,              kstar*h,
                   log(h+1)*kstar, ladj*kstar*h,       log(ladj)*log(h+1)*kstar, log(ladj)*h*kstar)
      if(alpha>=.2){
      cat("Using critical values for alpha = 0.2 \n")
      coefs <- c( 0.1485,  48.4,    1.75e-06,    0.7107, 
                 -1.039,    7.804,  0.1764,     13.62, 
                 -1078,   -11.72,   2.231e-06,  -1.69e-06, 
                 -18.06,  -28.15,  41.95,       -0.6027, 
                 -0.3967,   0.8376, 1.073,      -1.844)
      qres <- 0.0484
      } else if(alpha>=.1){
      cat("Using critical values for alpha = 0.1 \n")
      coefs <- c( 0.3369,  52.17,   1.508e-06,   0.6634, 
                 -0.9903,   9.034,  0.1796,     12.64, 
                 -1011, -  13.68,   2.194e-06,  -1.397e-06, 
                 -18.13,  -29.41,  43.75,       -0.555, 
                 -0.4388,   0.8297, 1.106,      -1.894)
      qres <- 0.049
      } else if(alpha>=.08) {
      cat("Using critical values for alpha = 0.08 \n")
      coefs <- c( 0.4256,  53.41,   1.455e-06,   0.6179, 
                 -0.9293,   9.289,  0.1804,     12.28, 
                 -987.2,  -14.62,   2.152e-06,  -1.348e-06, 
                 -17.99,  -29.35,  43.79,       -0.5383, 
                 -0.4485,   0.8227, 1.106,      -1.895)
      qres <- 0.0487
      } else if(alpha>=.06) {
      cat("Using critical values for alpha = 0.06 \n")
      coefs <- c( 0.5601,  55.66,   1.393e-06,   0.5289, 
                 -0.805,    9.886,  0.1815,     11.02, 
                 -971.4,  -15.86,   2.141e-06,  -1.318e-06, 
                 -17.23,  -30.27,  43.81,       -0.4877, 
                 -0.4705,   0.7915, 1.142,      -1.895)
      qres <- 0.048
      } else if(alpha>=.05) {
      cat("Using critical values for alpha = 0.05 \n")
      coefs <- c( 0.6215,  56.87,   1.359e-06,   0.4963, 
                 -0.7564,  10.2,    0.1822,     10.69, 
                 -954.4,  -16.67,   2.161e-06,  -1.31e-06, 
                 -17.16,  -30.23,  43.88,       -0.4751,  
                 -0.4826,  0.7901,  1.141,      -1.901)
      qres <- 0.0474
      } else if(alpha>=.04) {
      cat("Using critical values for alpha = 0.04 \n")
      coefs <- c( 0.6452,  58.79,   1.341e-06,   0.4984, 
                 -0.7409,  10.32,   0.1832,     10.28, 
                 -947.7,  -17.74,   2.152e-06,  -1.346e-06, 
                 -16.89,  -29.4,   43.13,       -0.458, 
                 -0.4863,   0.7782, 1.104,      -1.867)
      qres <- 0.0486
      } else if(alpha>=.025){
      cat("Using critical values for alpha = 0.025 \n")
      coefs <- c( 0.8513,  60.5,    1.292e-06,   0.3913, 
                 -0.583,   10.85,   0.184,       9.247, 
                 -862.6,  -19.84,   2.1e-06,    -1.244e-06, 
                 -16.33,  -29.41,  42.8,        -0.4224, 
                 -0.5051,   0.7618, 1.1,        -1.86)
      qres <- 0.049
      } else if(alpha>=.02){
      cat("Using critical values for alpha = 0.02 \n")
      coefs <- c( 0.8986,  60.42,   1.298e-06,   0.3855, 
                 -0.5647,  11.18,   0.1848,      8.507, 
                 -794.8,  -20.65,   2.081e-06,  -1.247e-06, 
                 -15.88,  -28.92,  42.07,       -0.3935, 
                 -0.5177,  0.7447,  1.077,      -1.831)
      qres <- 0.0516
      } else if(alpha>=.01){
      cat("Using critical values for alpha = 0.01 \n")
      coefs <- c( 1.143,   60.28,   1.236e-06,   0.301, 
                 -0.4117,  11.68,   0.1857,      5.67, 
                 -589.9,  -23.35,   2.3e-06,    -1.302e-06, 
                 -13.63,  -27.54,  38.9,        -0.285, 
                 -0.537,    0.661,  1.02,       -1.712)
      qres <- 0.0557
      } else {
        warning("Significance levels smaller then 0.01 are not implemented,
         using alpha = 0.01 instead")
      coefs <- c( 1.143,   60.28,   1.236e-06,   0.301, 
                 -0.4117,  11.68,   0.1857,      5.67, 
                 -589.9,  -23.35,   2.3e-06,    -1.302e-06, 
                 -13.63,  -27.54,  38.9,        -0.285, 
                 -0.537,    0.661,  1.02,       -1.712)
      qres <- 0.0557
      }
#  add 0.9 quantile of residuals
      sum(coefs*explvar)+qres
   }
#
# first check arguments and initialize
#
   args <- match.call()
# test dimension of data (vector of 3D) and define dimension related stuff
   d <- 3
   dy <- dim(y)
   n1 <- dy[1]
   n2 <- dy[2]
   n3 <- dy[3]
   n <- n1*n2*n3
   nt <- ddim[4]
   if (length(dy)==d+1) {
      dim(y) <- dy[1:3]
   } else if (length(dy)!=d) {
      stop("y has to be 3 dimensional")
   }
# set the code for the kernel (used in lkern) and set lambda
   lkern <- 1
   skern <- 1
# define lambda
   lambda <- ladjust*16
# to have similar preformance compared to skern="Exp"
   spmin <- .25 
   spmax <- 1
   if (is.null(hinit)||hinit<1) hinit <- 1
  
# define hmax
   if (is.null(hmax)) hmax <- 5    # uses a maximum of about 520 points


# define hincr
# determine corresponding bandwidth for specified correlation
   if(is.null(h0)) h0 <- rep(0,3)

# estimate variance in the gaussian case if necessary  
# deal with homoskedastic Gaussian case by extending sigma2
   if (length(sigma2)==1) sigma2<-array(sigma2,dy[1:3]) 
   if (length(sigma2)!=n) stop("sigma2 does not have length 1 or same length as y")
   dim(sigma2) <- dy[1:3]
   if(is.null(mask)) mask <- array(TRUE,dy[1:3])
   mask[sigma2>=1e16] <- FALSE
#  in these points sigma2 probably contains NA's
   sigma2 <- 1/sigma2 #  taking the invers yields simpler formulaes 
# deal with homoskedastic Gaussian case by extending sigma2
   residuals <- readBin(res,"integer",prod(ddim),2)
  cat("\nfmri.smooth: first variance estimate","\n")
  vartheta0 <- .Fortran("ivar",as.double(residuals),
                           as.double(resscale),
                           as.logical(mask),
                           as.integer(ddim[1]),
                           as.integer(ddim[2]),
                           as.integer(ddim[3]),
                           as.integer(ddim[4]),
                           var = double(n1*n2*n3),
                           PACKAGE="fmri",DUP=FALSE)$var
   varest0 <- vartheta0/nt
   vq <- varest0*sigma2
# Initialize  list for bi and theta
   if (is.null(wghts)) wghts <- c(1,1,1)
   hinit <- hinit/wghts[1]
   hmax <- hmax/wghts[1]
   wghts <- (wghts[2:3]/wghts[1])
   tobj <- list(bi= rep(1,n))
   theta <- y
   segm <- array(0,dy[1:3])
   varest <- varest0
   fix <- array(FALSE,dy[1:3])
   maxvol <- getvofh(hmax,lkern,wghts)
   if(is.null(fov)) fov <- sum(mask)
   kstar <- as.integer(log(maxvol)/log(1.25))
   steps <- kstar+1
   k <- 1 
   hakt <- hinit
   hakt0 <- hinit
   lambda0 <- lambda
   if (hinit>1) lambda0 <- 1e50 # that removes the stochstic term for the first step
   scorr <- numeric(3)
   if(h0[1]>0) scorr[1] <-  get.corr.gauss(h0[1],2)
   if(h0[2]>0) scorr[2] <-  get.corr.gauss(h0[2],2)
   if(h0[3]>0) scorr[3] <-  get.corr.gauss(h0[3],2)
   total <- cumsum(1.25^(1:kstar))/sum(1.25^(1:kstar))
   beta <- getbeta(df,h0)
   thresh <- getkrval(df,h0,ladjust,n,kstar,alpha)
   cat("FOV",fov,"delta",delta,"thresh",thresh,"ladjust",ladjust,"lambda",lambda,"df",df,"\n")
# run single steps to display intermediate results
   residuals <- residuals*resscale
#
#   need these values to compute variances 
#
   while (k<=kstar) {
      hakt0 <- gethani(1,10,lkern,1.25^(k-1),wghts,1e-4)
      hakt <- gethani(1,10,lkern,1.25^k,wghts,1e-4)
      cat("step",k,"bandwidth",signif(hakt,3)," ")
      dlw <- (2*trunc(hakt/c(1,wghts))+1)[1:d]
#  need bandwidth in voxel for Spaialvar.gauss, h0 is in voxel
      if (any(h0>0)) lambda0 <- lambda0 * Spatialvar.gauss(bw2fwhm(hakt0)/4/c(1,wghts),h0,d)/
      Spatialvar.gauss(h0,1e-5,d)/Spatialvar.gauss(bw2fwhm(hakt0)/4/c(1,wghts),1e-5,d)
# Correction C(h0,hakt) for spatial correlation depends on h^{(k-1)}  all bandwidth-arguments in FWHM 
      hakt0 <- hakt
      theta0 <- theta
      bi0 <- tobj$bi
      tobj <- .Fortran("segm3d",
                       as.double(y),
                       fix=as.logical(fix),
                       as.double(residuals),
                       as.double(sigma2),
                       as.logical(!mask),
                       as.logical(weighted),
                       as.integer(n1),
                       as.integer(n2),
                       as.integer(n3),
                       as.integer(nt),
                       hakt=as.double(hakt),
                       as.double(lambda0),
                       as.double(theta0),
                       bi=as.double(bi0),
                       thnew=double(n1*n2*n3),
                       double(prod(dlw)),
                       as.double(wghts),
                       double(nt),#swres
                       double(n1*n2*n3),#pvalues
                       segm=as.integer(segm),
                       as.double(beta),
                       as.double(delta),
                       as.double(thresh),
                       as.integer(k),
                       as.double(fov),
                       as.double(vq),
                       as.double(varest0),
                       varest=as.double(varest),
                       PACKAGE="fmri",DUP=FALSE)[c("bi","thnew","hakt","segm","varest","fix")]
      gc()
      theta <- array(tobj$thnew,dy[1:3]) 
      fix <- array(tobj$fix,dy[1:3]) 
      segm <- array(tobj$segm,dy[1:3])
      varest <- array(tobj$varest,dy[1:3])
      dim(tobj$bi) <- dy[1:3]
      if (graph) {
         par(mfrow=c(2,2),mar=c(1,1,3,.25),mgp=c(2,1,0))
         image(y[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Observed Image  min=",signif(min(y),3)," max=",signif(max(y),3)))
         image(theta[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Reconstruction  h=",signif(hakt,3)," min=",signif(min(theta),3),"   max=",signif(max(theta),3)))
         image(segm[,,n3%/%2+1]>0,col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Segmentation  h=",signif(hakt,3)," detected=",sum(segm>0)))
         image(tobj$bi[,,n3%/%2+1],col=gray((0:255)/255),xaxt="n",yaxt="n")
         title(paste("Sum of weights: min=",signif(min(tobj$bi),3)," mean=",signif(mean(tobj$bi),3)," max=",signif(max(tobj$bi),3)))
      }
      if (max(total) >0) {
         cat(signif(total[k],2)*100,"%                 \r",sep="")
      }
      k <- k+1
#  adjust lambda for the high intrinsic correlation between  neighboring estimates 
      c1 <- (prod(h0+1))^(1/3)
      c1 <- 2.7214286 - 3.9476190*c1 + 1.6928571*c1*c1 - 0.1666667*c1*c1*c1
      x <- (prod(1.25^(k-1)/c(1,wghts)))^(1/3)
      scorrfactor <- (c1+x)/(c1*prod(h0+1)+x)
      lambda0 <- lambda*scorrfactor
      gc()
   }

  z <- list(theta=theta,ni=tobj$bi,var=varest,y=y,segm=segm,
            hmax=tobj$hakt*switch(lkern,1,1,bw2fwhm(1/4)),
            call=args,scorr=scorr,mask=mask)
  class(z) <- "aws.gaussian"
  invisible(z)
}

