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
      explvar <- c(1,    log(h+1),     h,     kstar,      sqrt(log(n)),   1/df,
                   log(ladj)/df, log(h+1)/df, log(ladj)*log(h+1), sqrt(log(n))/df)
      if(alpha>=.25){
      cat("Using critical values for alpha = 0.25 \n")
      coefs <- c(-1.759, 6.085, -3.159, 0.00257, 0.3874, -27.11, -30.04, 24.78, -1.553, 21.44)
      qres <- 0.0592  # res.quantiles(.9,.95,.99,1) 0.0300 0.0356 0.0445 0.0592
      } else if(alpha>=.2){
      cat("Using critical values for alpha = 0.2 \n")
      coefs <- c(-1.56, 6.076, -3.15, 0.003117, 0.3441, -26.21, -30.21, 25.6, -1.559, 21.73)
      qres <- 0.0634  # res.quantiles(.9,.95,.99,1) 0.0303 0.0365 0.0468 0.0634
      } else if(alpha>=.15){
      cat("Using critical values for alpha = 0.15 \n")
      coefs <- c(-1.346, 6.061, -3.127, 0.003738, 0.2988, -24.54, -30.29, 25.84, -1.556, 21.98)
      qres <- 0.0732  # res.quantiles(.9,.95,.99,1) 0.0310 0.0377 0.0476 0.0732
      } else if(alpha>=.125){
      cat("Using critical values for alpha = 0.125 \n")
      coefs <- c(-1.266, 6.075, -3.132, 0.004119, 0.2861, -19.34, -30.4, 26, -1.556, 20.9)
      qres <- 0.0691  # res.quantiles(.9,.95,.99,1) 0.0326 0.0392 0.0492 0.0691
      } else if(alpha>=.1){
      cat("Using critical values for alpha = 0.1 \n")
      coefs <- c(-1.036, 6.077, -3.126, 0.004547, 0.2316, -21.86, -30.78, 26.44, -1.571, 22.21)
      qres <- 0.0749  # res.quantiles(.9,.95,.99,1) 0.0342 0.0415 0.0538 0.0749
      } else if(alpha>=.075){
      cat("Using critical values for alpha = 0.075 \n")
      coefs <- c(-0.7676, 6.11, -3.147, 0.005197, 0.1678, -27.27, -31.17, 27.47, -1.58, 24.55)
      qres <- 0.0904  # res.quantiles(.9,.95,.99,1) 0.0373 0.0477 0.0654 0.0904
      } else if(alpha>=.05){
      cat("Using critical values for alpha = 0.05 \n")
      coefs <- c(-0.4409, 6.1, -3.129, 0.006261, 0.09223, -31.3, -30.87, 28.65, -1.652, 26.68)
      qres <- 0.0917  # res.quantiles(.9,.95,.99,1) 0.0419 0.0545 0.0764 0.0917
      } else if(alpha>=.02){
      cat("Using critical values for alpha = 0.02 \n")
      coefs <- c(-0.1448, 6.076, -3.09, 0.00857, 0.04783, -22.23, -30.81, 31.67, -1.848, 26.27)
      qres <- 0.105  # res.quantiles(.9,.95,.99,1) 0.0522 0.0663 0.0861 0.1050
      } else if(alpha>=.01){
      cat("Using critical values for alpha = 0.01 \n")
      coefs <- c(0.3636, 6.102, -3.085, 0.01067, -0.07914, -46.55, -31.02, 33.09, -1.999, 35.57)
      qres <- 0.149  # res.quantiles(.9,.95,.99,1) 0.0626 0.0784 0.1100 0.1490
      } else if(alpha>=.005){
      cat("Using critical values for alpha = 0.005 \n")
      coefs <- c(0.3803, 6.179, -3.102, 0.01388, -0.07502, -32.43, -30.45, 35.01, -2.143, 33.69)
      qres <- 0.213  # res.quantiles(.9,.95,.99,1) 0.0783 0.1040 0.1560 0.2130
      } else if(alpha>=.002){
      cat("Using critical values for alpha = 0.002 \n")
      coefs <- c(0.2181, 6.374, -3.127, 0.02009, -0.03312, -12.31, -30.59, 31.57, -2.374, 30.89)
      qres <- 0.211  # res.quantiles(.9,.95,.99,1) 0.112 0.137 0.189 0.211
      } else if(alpha>=.001){
      cat("Using critical values for alpha = 0.001 \n")
      coefs <- c(-0.0509, 6.654, -3.255, 0.02556, 0.03278, 2.36, -29.05, 29.92, -2.466, 29.13)
      qres <- 0.416  # res.quantiles(.9,.95,.99,1) 0.128 0.167 0.253 0.416
      } else {
      cat("Using critical values for alpha = 0.001 \n")
      coefs <- c(-0.0509, 6.654, -3.255, 0.02556, 0.03278, 2.36, -29.05, 29.92, -2.466, 29.13)
      qres <- 0.416  # res.quantiles(.9,.95,.99,1) 0.128 0.167 0.253 0.416
      }
#  add maximum of residuals
      print(rbind(coefs,explvar,coefs*explvar))
      sum(coefs*explvar)+qres
   }
#
#  set beta depending on df
#
   beta <- 2+30/df+285/df^2+5000/df^3
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
   thresh <- getkrval(df,h0,ladjust,fov,kstar,alpha)
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

