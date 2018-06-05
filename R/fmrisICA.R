fmri.sICA <- function(data, mask=NULL, ncomp=20,
               alg.typ=c("parallel","deflation"),
               fun=c("logcosh","exp"), alpha=1, detrend=TRUE,
               degree=2, nuisance= NULL, smooth=TRUE,
               bw=8, unit=c("SD","FWHM")){
  if(detrend) data <- fmri.detrend(data,degree,nuisance)
  bwvoxel <- if(length(data$header$pixdim)>4)
        bw/data$header$pixdim[2:4] else bw
  if(smooth) data <- smooth.fmridata(data,bwvoxel,unit)
  cat("Computing independent components with fastICA \n")
  ttt <- extract.data(data)
  if(is.null(mask)) mask <- data$mask
  if(!is.logical(mask)) mask <- as.logical(mask)
  ddim <- data$dim[1:3]
  dim(ttt) <- c(prod(ddim),dim(ttt)[4])
  ttt0 <- ttt[mask,]
  fsica <- fastICA::fastICA(ttt0,n.comp=ncomp,
       alg.typ=alg.typ, fun=fun, alpha=alpha, method="C",
       maxit=500, tol=1e-5)
  cimgs <- array(0,c(prod(ddim),ncomp))
  cimgs[mask,] <- fsica$S
  dim(cimgs) <- c(ddim,ncomp)
  list(scomp=cimgs,X=fsica$X,k=fsica$K,W=fsica$W,A=fsica$A,
       mask=mask,pixdim=data$header$pixdim[2:4])
  }

ICAfingerprint <- function(icaobj,nbin=256,plot=FALSE){
##
##  Martino-Goebel(2006)
##
   scomp <- icaobj$scomp
   dscomp <- dim(scomp)
   ncomp <- dscomp[4]
   dim(scomp) <- c(prod(dscomp[1:3]),ncomp)
   scomp <- scomp[icaobj$mask,]
   nvox <- dim(scomp)[1]
   tcomp <- icaobj$A
   ntime <- dim(tcomp)[2]
   icafp <- matrix(0,11,ncomp)
   dimnames(icafp) <- list(c("kurt","skew","sentropy","dclust",
                        "ar1","tentropy","power1","power2",
                        "power3","power4","power5"),
                        paste0("Component",1:ncomp))
   for(i in 1:ncomp){
      sdcompi <- sd(scomp[,i])
      scomp[,i] <- scomp[,i]/sdcompi
      icafp[1,i] <- sum(scomp[,i]^4)/nvox - 3
      icafp[2,i] <- sum(scomp[,i]^3)/nvox
      z <- hist(scomp[,i],nbin,plot=FALSE)$counts
      z <- z[z>0]/sum(z)
      icafp[3,i] <- -sum(z*log(z,2))
      scompi <- icaobj$scomp[,,,i]
      scompi <- abs(scompi-mean(scompi))/sdcompi
      z <- findclusters(scompi,2.5)
      icafp[4,i] <- sum(z$size>270/prod(icaobj$pixdim))/nvox
      icafp[5,i] <- abs(sum(tcomp[i,-1]*tcomp[i,-ntime]))/sum(tcomp[i,]^2)*ntime/(ntime-1)
      z <- hist(tcomp[i,],nbin,plot=FALSE)$counts
      z <- z[z>0]/sum(z)
      icafp[6,i] <- -sum(z*log(z,2))
      z <- spectrum(tcomp[i,],plot=FALSE)
      ns <- length(z$freq)
      ind1 <- (1:ns)[z$freq<=.008]
      ind2 <- (1:ns)[z$freq>.008&z$freq<=.02]
      ind3 <- (1:ns)[z$freq>.02&z$freq<=.05]
      ind4 <- (1:ns)[z$freq>.05&z$freq<=.1]
      ind5 <- (1:ns)[z$freq>.1&z$freq<=.25]
      sumsp <- sum(z$spec)
      icafp[7,i] <- sum(z$spec[ind1])/sumsp
      icafp[8,i] <- sum(z$spec[ind2])/sumsp
      icafp[9,i] <- sum(z$spec[ind3])/sumsp
      icafp[10,i] <- sum(z$spec[ind4])/sumsp
      icafp[11,i] <- sum(z$spec[ind5])/sumsp
   }
   ## original normalizations from Martino/Goebel are invalid
   ## replaced by normalizations that are monotone and don't produce NaN's
   icafp[1,] <- log(icafp[1,]-min(icafp[1,])+1)
   icafp[1,] <- icafp[1,]/max(icafp[1,])
   icafp[2,] <- log(icafp[2,]-min(icafp[2,])+1)
   icafp[2,] <- icafp[2,]/max(icafp[2,])
   icafp[3,] <- log(icafp[3,]+1)
   icafp[3,] <- icafp[3,]/max(icafp[3,])
   icafp[6,] <- log(icafp[6,]+1)
   icafp[6,] <- icafp[6,]/max(icafp[6,])
   icafp[7:11,] <- sweep(icafp[7:11,],2,apply(icafp[7:11,],2,max),"/")
   icafp[7:11,] <- sweep(icafp[7:11,],1,apply(icafp[7:11,],1,max),"/")
   if(plot) stars(t(icafp),key.xpd=NA)
   invisible(t(icafp))
}
