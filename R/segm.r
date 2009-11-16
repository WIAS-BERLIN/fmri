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
#    ladj in (1 : 1.4) 
#    kstar in (11:29)  corrsponding to maximal bandwidths  3 - 6 
#    internal adjustment by 0.1*kstar 
#    see file sim_fmri_kritval.r in R/segmentation/fmrikrv/
      kstar <- max(kstar,11)  #  use minimal kstar
      h <- fwhm2bw(sqrt(h0[1]*h0[2]))
      explvar <- c(1,    df^(-2.75),     1/df,     1/ladj,      kstar,   kstar^(-.5),
                   sqrt(h), log(n), 1/df/ladj, kstar/df, 1/sqrt(kstar)/df, sqrt(h)/df, log(n)/df)
      a <- c(1, alpha, alpha^2, sqrt(alpha), log(alpha), exp(alpha), exp(-alpha))
      acoef <- matrix(c( 9.930e-1,  6.713e-3,  0       ,  0       ,  0       ,  0       ,  0       ,
                         7.408e+2,  0       ,  0       ,  3.321e+2, -1.561e+2,  0       ,  0       ,
                        -7.435e+2,  0       , -4.135e+2, -2.537e+1,  0       ,  3.934e+2,  3.623e+2,
                        -6.816e-1,  6.083e-1, -2.266e-1,  3.585e-2,  0       ,  0       ,  6.731e-1,
                        -9.116e-3,  0       ,  0       ,  5.295e-3, -2.893e-4,  7.569e-4,  6.941e-3,
                        -8.356e-1, -1.466e+0,  0       ,  7.695e-1, -5.010e-2,  6.394e-1,  0       ,
                        -7.673e-4,  0       ,  0       ,  1.430e-3,  0       ,  0       ,  0       ,
                        -2.461e-3,  0       ,  2.315e-3,  0       ,  0       ,  0       ,  0       ,
                         2.512e+3, -7.747e+2,  1.186e+3, -8.260e+0,  0       , -8.556e+2, -1.653e+3,
                         4.776e-1,  1.047e+0,  0       , -5.422e-1,  3.536e-2, -4.287e-1,  0       ,
                         7.860e+1,  1.839e+2,  0       , -1.031e+2,  8.786e+0, -7.154e+1,  0       ,
                         3.407e+0, -2.543e-1,  0       ,  0       ,  1.517e-1,  0       ,  0       ,
                         8.652e-1,  0       ,  0       , -2.830e-1,  1.345e-2,  0       ,  0       ),7,13)
      dimnames(acoef) <- list(c("(Intercept)","a","a2","ah","al","ae","aei"),NULL)
      ecoefs <- t(acoef)%*%a
      t(explvar)%*%ecoefs
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
                       as.double(df),
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

