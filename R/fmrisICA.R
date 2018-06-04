fmri.sICA <- function(data, mask=NULL, ncomp=20, detrend=TRUE,
               degree=2, nuisance= NULL, smooth=TRUE,
               bw=8, unit=c("SD","FWHM")){
  if(detrend) data <- fmri.detrend(data,degree,nuisance)
  if(smooth) data <- smooth.fmridata(data,bw,unit)
  cat("Computing independent components with fastICA \n")
  ttt <- extract.data(data)
  if(is.null(mask)) mask <- data$mask
  if(!is.logical(mask)) mask <- as.logical(mask)
  ddim <- data$dim[1:3]
  dim(ttt) <- c(prod(ddim),dim(ttt)[4])
  ttt0 <- ttt[mask,]
  fsica <- fastICA::fastICA(ttt0,n.comp=ncomp,
       alg.typ="parallel",method="C")
  cimgs <- array(0,c(prod(ddim),ncomp))
  cimgs[mask,] <- fsica$S
  dim(cimgs) <- c(ddim,ncomp)
  list(scomp=cimgs,X=fsica$X,k=fsica$K,W=fsica$W,A=fsica$A,mask=mask)
  }
