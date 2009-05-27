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
      if(alpha>=.1){
      explvar <- c(1,df,n,ladj,log(ladj),log(h+1),kstar,h,1/df^2,df*ladj,
                  df*log(ladj),df*h,df*kstar,n*ladj,n*log(h+1),n/df^2,ladj*log(h+1),
                  ladj*kstar,ladj/df^2,log(ladj)*log(h+1),log(ladj)*kstar, 
                  log(h+1)*kstar,kstar*h,h/df^2,ladj*kstar*h)
      coefs <- c(10.74, 0.006418, 2.506e-06, -9.16, 10.98, -1.978, -0.2392, 0.07853, 
                 218.9, -0.008211, 0.01178, -0.0008796, -8.596e-05, -1.166e-06, 2.293e-06, 
                 0.0006664, 3.949, 0.428, -182.7, -6.049, -0.5332, -0.1624, 0.08483, -200.7, -0.0243)
      qres <- 0.048
      } else if(alpha>=.05) {
      explvar <- c(1,df,n,ladj,log(ladj),log(h+1),kstar,h,1/df^2,df*ladj,
                  df*log(ladj),df*h,df*kstar,n*ladj,n*log(h+1),n/df^2,ladj*log(h+1),
                  ladj*kstar,ladj/df^2,log(ladj)*log(h+1),log(ladj)*kstar, 
                  log(h+1)*kstar,kstar*h,h/df^2,ladj*kstar*h)
      coefs <- c(10.81, 0.01105, 2.405e-06, -9.067, 11.02, -1.185, -0.2431, -0.7695, 
                 426.9, -0.01335, 0.01822, -0.0008828, -8.345e-05, -1.189e-06, 2.334e-06, 
                 0.0006246, 4.192, 0.4362, -250.8, -5.718, -0.55, -0.2066, 0.139, -198.8, -0.04031)
      qres <- 0.049
      } else if(alpha>=.025){
      explvar <- c(1,df,n,ladj,log(ladj),log(h+1),kstar,h,1/df^2,df*ladj,
                  df*log(ladj),df*h,df*kstar,n*ladj,n*log(h+1),n/df^2,ladj*log(h+1),
                  ladj*kstar,ladj/df^2,log(ladj)*log(h+1),log(ladj)*kstar, 
                  log(h+1)*kstar,kstar*h,h/df^2,ladj*kstar*h)
      coefs <- c(10.49, 0.01185, 2.138e-06, -8.508, 10.36, -1.221, -0.2212, -1.089, 
                 553.5, -0.01539, 0.02176, -0.0007787, -5.232e-05, -1.013e-06, 2.129e-06, 
                 0.0007087, 4.567, 0.415, -253.5, -5.852, -0.5289, -0.2222, 0.1647, -178.5, -0.0495)
      qres <- 0.056
      } else if(alpha>=.01){
      explvar <- c(1,df,n,ladj,log(ladj),log(h+1),kstar,h,1/df^2,df*ladj,
                  df*h,n*ladj,n*log(h+1),ladj*log(h+1),
                  ladj*kstar,ladj/df^2,log(ladj)*log(h+1),log(ladj)*kstar, 
                  log(h+1)*kstar,kstar*h,h/df^2,ladj*kstar*h)
      coefs <-c(11.2, -0.006744, 2.364e-06, -9.002, 11.19, -1.256, -0.189, -2.164,
                941.2, 0.002191, -0.001276, -1.023e-06, 2.073e-06, 5.819, 0.3824, 
                -306.5, -6.634, -0.4981, -0.2734, 0.2338, -223.8, -0.06904)
      qres <- 0.073
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
   varquot <- mean(vartheta0/nt*sigma2)
# Initialize  list for bi and theta
   if (is.null(wghts)) wghts <- c(1,1,1)
   hinit <- hinit/wghts[1]
   hmax <- hmax/wghts[1]
   wghts <- (wghts[2:3]/wghts[1])
   tobj <- list(bi= rep(1,n))
   theta <- y
   segm <- array(0,dy[1:3])
   varest <- 1/sigma2
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
#
#   need these values to compute variances after the last iteration
#
      residuals <- readBin(res,"integer",prod(ddim),2)*resscale
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
                       as.double(varquot),
                       as.double(delta),
                       as.double(thresh),
                       as.integer(k),
                       as.double(fov),
                       varest=as.double(varest),
                       PACKAGE="fmri",DUP=FALSE)[c("bi","thnew","hakt","segm","varest","fix")]
      gc()
      theta <- array(tobj$thnew,dy[1:3]) 
      fix <- array(tobj$fix,dy[1:3]) 
      segm <- array(tobj$segm,dy[1:3])
     par(mfrow=c(1,1))
      image(apply(segm>0,1:2,sum),zlim=c(0,26))
#      cat("\n",sum(signal&(segm==1)),sum(signal&(segm!=1)),sum(!signal&(segm==1)),sum(!signal&(segm!=1)),"\n")
      readline("Press enter")
      varest <- array(tobj$varest,dy[1:3])
#      varest <- array(tobj$varest/vartheta0/sigma2,dy[1:3])
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
            call=args,scorr=scorr)
  class(z) <- "aws.gaussian"
  invisible(z)
}
