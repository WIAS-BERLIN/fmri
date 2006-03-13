fmri.smooth <- function(spm,hmax=4,type="",...) {
  cat("fmri.smooth: entering function\n")
  
  if (!class(spm) == "fmrispm") {
    warning("fmri.smooth: data not of class <fmrispm>. Try to proceed but strange things may happen")
  }

  if (!is.null(spm$smooth)) {
    if (spm$smooth) {
      warning("fmri.smooth: Parametric Map seems to be smoothed already!")
    }
  }
  
  variance <- spm$var
#  variance[variance < quantile(variance,0.25)] <- quantile(variance,0.25)
  variance[variance == 0] <- 1e20


  if (is.null(spm$weights)) {
    weights <- c(1,1,1)
  } else {
    weights <- spm$weights
  }
  if (is.null(spm$scorr)) {
    scorr <- 0
  } else {
    scorr <- spm$scorr
  }

  if (type == "old") {
    cat("fmri.smooth: smoothing the Statistical Paramteric Map\n")
    ttthat <- vaws3Dold(y=spm$cbeta, sigma2=variance, hmax=hmax, wghts=weights, scorr=scorr, qtau=1, vwghts = spm$vwghts, ...)
    cat("\n")
    cat("fmri.smooth: determine local smoothness\n")
    bwx <- get.bw.gauss(corr[1])
    bwy <- get.bw.gauss(corr[2])
    bwz <- get.bw.gauss(corr[3])
    bw <- get3Dh.gauss.old(ttthat$vred,c(bwx,bwy,bwz),weights)
    rxyz <- c(resel(1,bw[,1]), resel(1,bw[,2]), resel(1,bw[,3]))
    dim(rxyz) <- c(dim(bw)[1],3)
  } else {
    cat("fmri.smooth: smoothing the Statistical Paramteric Map\n")
    ttthat <- vaws3D(y=spm$cbeta, sigma2=variance, hmax=hmax, wghts=weights, scorr=scorr, vwghts = spm$vwghts, ...)
    cat("\n")
    
    cat("fmri.smooth: determine local smoothness\n")
    bw <- get3Dh.gauss(ttthat$vred,weights)
    rxyz <- c(resel(1,bw[,1]), resel(1,bw[,2]), resel(1,bw[,3]))
    dim(rxyz) <- c(dim(bw)[1],3)
  }

  cat("fmri.smooth: exiting function\n")
    
  if (dim(ttthat$theta)[4] == 1) {
    z <- list(cbeta = ttthat$theta[,,,1], var = ttthat$var, rxyz =
              rxyz, scorr = spm$scorr, weights = spm$weights, vwghts = spm$vwghts, smooth =
              TRUE, hmax = ttthat$hmax, dim = spm$dim)
  } else {
    z <- list(cbeta = ttthat$theta, var = ttthat$var, rxyz = rxyz,
              scorr = spm$scorr, weights = spm$weights, vwghts = spm$vwghts, smooth = TRUE,
              hmax = ttthat$hmax, dim = spm$dim)
  }    

  class(z) <- "fmrispm"
  z
}

fmri.detect <- function(spm, maxpvalue = 0.05, mode="plog", global=FALSE, delta=NULL) {
  cat("fmri.detect: entering function\n")

  if (!class(spm) == "fmrispm") {
    warning("fmri.detect: data not of class <fmrispm>. Try to proceed but strange things may happen")
  }

  if (length(dim(spm$cbeta)) < 4) {
    stat <- spm$cbeta/sqrt(spm$var)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.detect: calculate treshold\n")
    thresh <- threshold(maxpvalue,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm")
    cat("fmri.detect: calculate p-value\n")
    pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm")
  } else if (!is.null(delta)) {
    l1 <- sqrt(spm$vwghts[2]/spm$vwghts[1]) * delta[1]
    l2 <- sqrt(spm$vwghts[2]/spm$vwghts[1]) * delta[2]
    theta1 <- atan(l1)
    theta2 <- atan(l2)
    t1 <- spm$cbeta[,,,1]/sqrt(spm$var * spm$vwghts[1])
    t2 <- spm$cbeta[,,,2]/sqrt(spm$var * spm$vwghts[2])
    ratio <- t2/t1
    w1 <- (t1 + t2 * l1) / sqrt(1+l1^2)
    w2 <- (t1 + t2 * l2) / sqrt(1+l2^2)
    w3 <- (t1 > 0) * (l1 <= ratio) * (ratio <= l2) * sqrt(t1^2 + t2^2)
    stat <- pmax(w1,w2,w3)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.detect: calculate treshold\n")
    thresh <- 
      thresholdc(maxpvalue,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm",cone=theta2-theta1)
    cat("fmri.detect: calculate p-value\n")
    pv <-
      pvaluec(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm",cone=theta2-theta1)
  } else {
    stat <- spm$cbeta[,,,1]^2/spm$var + spm$cbeta[,,,2]^2/spm$var/spm$vwghts[2]  # Wert der Statistik
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.detect: calculate treshold\n")
    thresh <- 
      threshold(maxpvalue,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="chisq",df=2)
    cat("fmri.detect: calculate p-value\n")
    pv <-
      pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="chisq",df=2)
  }
  cat("fmri.detect: thresholding\n")
  mask <- rep(1,length=prod(spm$dim[1:3]))
  if (global) thresh <- median(thresh)
  mask[stat < thresh] <- 0
  pv[!mask] <- 1
  pv[pv <= 1e-113] <- 1e-113 # ????????
  pvlog <- -log(pv)
  pvlog[mask == 0] <- 0
  dim(pvlog) <- spm$dim[1:3]
  dim(mask) <- spm$dim[1:3]
  
  cat("fmri.detect: exiting function\n")

  z <- list(detect = pvlog)
  
  class(z) <- "fmridetect"
  z
}



#
# plot the (logical) difference between mask1 and mask2 as three
# color overlay over anatomic image:
# red for mask1 && mask2
# green for mask1 && !mask2
# blue for !mask1 && mask2
# mark position pos with lines
#
fmriplot.diff <- function(mask1, mask2, anatomic = array(0,dim=dim(mask1)), pos = c(-1,-1,-1), device="X11", file="plot.png") {
  partition <- ceiling(sqrt(dim(anatomic)[3]))

  switch (device,
          "png" = png(filename=file, width = 200*partition, height = 200*partition, pointsize=12, bg="transparent", res=NA),
          "jpeg" = jpeg(filename=file, width = 200*partition, height = 200*partition,
            quality = 100, pointsize = 12, bg = "transparent", res=NA),
          "ppm" = bitmap(file,type="ppm",height=2*partition,width=2*partition,res=64,pointsize=12),
          X11())

  oldpar <- par(mfrow=c(partition,partition), mar=c(0.25,0.25,0.25,.25), mgp=c(2,1,0))
  
  firstmask <- mask1 * (1-mask2)
  secondmask <- (1-mask1) * mask2
  bothmask <- mask1 * mask2

  alim <- range(anatomic)
  for (i in 1:dim(anatomic)[3]) {
    image(anatomic[,,i], xaxt="n", yaxt="n", zlim=alim, col=grey(0:255/255))
    if (any(bothmask[,,i])) image(bothmask[,,i], col=c(0,"#FF0000"), add=TRUE)
    if (any(firstmask[,,i])) image(firstmask[,,i], col=c(0,"#00FF00"), add=TRUE)
    if (any(secondmask[,,i])) image(secondmask[,,i], col=c(0,"#0000FF"), add=TRUE)

    if (i == pos[3]) {
      lines(c(0,1),c(pos[2],pos[2])/dim(anatomic)[2],col=2)
      lines(c(pos[1],pos[2])/dim(anatomic)[1],c(0,1),col=2)
    }     
  }
    
  par(oldpar)
  switch (device,
          "png" = dev.off(),
          "jpeg" = dev.off(),
          "ppm" = dev.off())
}

#
# plot the signal (which may be pvalue or signal itself or any other
# varying quantity) as overlay over anatomic image
# mark position pos with lines
#
fmriplot.signal <- function(signal, anatomic = array(0,dim=dim(signal)), pos = c(-1,-1,-1), device="X11", file="plot.png") {
  partition <- ceiling(sqrt(dim(anatomic)[3]))

  switch (device,
          "png" = png(filename=file, width = 200*partition, height = 200*partition, pointsize=12, bg="transparent", res=NA),
          "jpeg" = jpeg(filename=file, width = 200*partition, height = 200*partition,
            quality = 100, pointsize = 12, bg = "transparent", res=NA),
          "ppm" = bitmap(file,type="ppm",height=2*partition,width=2*partition,res=64,pointsize=12),
          X11())

  oldpar <- par(mfrow=c(partition,partition), mar=c(0.25,0.25,0.25,.25), mgp=c(2,1,0))
  
  alim <- range(anatomic)
  zlim <- range(signal)
  for (i in 1:dim(anatomic)[3]) {
    image(anatomic[,,i], xaxt="n", yaxt="n", zlim=alim, col=grey(0:255/255))
    if (any(signal[,,i])) image(signal[,,i], zlim=zlim ,col=c(0,heat.colors(512)), add=TRUE)

    if (i == pos[3]) {
      lines(c(0,1),c(pos[2],pos[2])/dim(anatomic)[2],col=2)
      lines(c(pos[1],pos[2])/dim(anatomic)[1],c(0,1),col=2)
    }     
  }
    
  par(oldpar)
  switch (device,
          "png" = dev.off(),
          "jpeg" = dev.off(),
          "ppm" = dev.off())
}

#
# plot timeseries given in fmri (2+1-dim!) use the fit to draw lines
# through data. print signal value(s) and value(s) of t statistic in
# sigvec and tvec (2+1-dim!) for any voxel. print voxel position
# starting from origin. mark active voxel given by threshold with red
# data points. 
fmriplot.timeseries <- function(fmri, fit, sigvec, tvec, thr, origin, device="X11", file="plot.png") {
  partition <- dim(fmri)[2:1]
  switch (device,
          "png" = png(filename=file, width = 200*partition[2], height = 200*partition[1], pointsize=12, bg="transparent", res=NA),
          "jpeg" = jpeg(filename=file, width = 200*partition[2], height = 200*partition[1],
            quality = 100, pointsize = 12, bg = "transparent", res=NA),
          "ppm" = bitmap(file,type="ppm",height=2*partition[1],width=2*partition[2],res=64,pointsize=12),
          X11())

  oldpar <- par(mfrow=partition, mar=c(0.25,0.25,0.25,.25), mgp=c(2,1,0))

  ylim <- range(fmri)
  for (j in 1:dim(fmri)[2]) {
    for (i in 1:dim(fmri)[1]) {
      if (tvec[i,j,1] > thr) {
        plot(fmri[i,j,],cex=0.5,xaxt="n",yaxt="n",col="red")
      } else {
        plot(fmri[i,j,],cex=0.5,xaxt="n",yaxt="n")
      }
      lines(fit[i,j,],col=3)
      text(10,max(fmri[i,j,]),pos=4,paste(origin[1]+i-1,origin[2]+j-1,sep="/"))    
      text(10,min(fmri[i,j,]),pos=4,
           paste(paste(signif(sigvec[i,j,],digits=2),signif(tvec[i,j,],digits=2),sep="/"),collapse="---"))
    }
  }
  
  
  par(oldpar)
  switch (device,
          "png" = dev.off(),
          "jpeg" = dev.off(),
          "ppm" = dev.off())
}

fmriplot.slices3d <- function(signal, anatomic =
                              array(0,dim=dim(signal)), pos =
                              c(-1,-1,-1)) {

  if (require(misc3d)) {

    mri.colors <- function (n1, n2, factor=n1/(n1+n2), from=0, to=.2) {
      colors1 <- gray((0:n1)/(n1+n2))
      colors2 <- hsv(h = seq(from,to,length=n2),
                     s = seq(from = n2/(n2+factor*n1) - 1/(2 * (n2+factor*n1)), to =
                       1/(2 * (n2+factor*n1)), length = n2),
                     v = 1,
                     gamma=1)
      list(all=c(colors1,colors2),gray=colors1,col=colors2)
    }

    # re-scale anatomic to 0 ... 0.5
    anatomic <- 0.5 * (anatomic - min(anatomic)) / diff(range(anatomic))

    # re-scale signal to 0.5 ... 1
    signal <-  0.5 + 0.5 * (signal - min(signal)) / diff(range(signal))

    anatomic[signal > 0.5] <- signal[signal > 0.5]

    slices3d(anatomic,col=mri.colors(256,256)$all)
    
  } else {
    stop("misc3d is required.")
  }
  
}

plot.fmridetect <- function(x, anatomic =
                             array(0,dim=dim(x$detect)),
                             pos = c(-1,-1,-1), type="slice",
                             device="X11", file="plot.png",...) {
  mri.colors <- function (n1, n2, factor=n1/(n1+n2), from=0, to=.2) {
    colors1 <- gray((0:n1)/(n1+n2))
    colors2 <- hsv(h = seq(from,to,length=n2),
                   s = seq(from = n2/(n2+factor*n1) - 1/(2 * (n2+factor*n1)), to =
                     1/(2 * (n2+factor*n1)), length = n2),
                   v = 1,
                   gamma=1)
    list(all=c(colors1,colors2),gray=colors1,col=colors2)
  }

  signal <- x$detect
  
  if (type == "3d") {
    if (require(misc3d)) {
                                       # re-scale anatomic to 0 ... 0.5
      anatomic <- 0.5 * (anatomic - min(anatomic)) / diff(range(anatomic))
      
                                        # re-scale signal to 0.5 ... 1
      signal <-  0.5 + 0.5 * (signal - min(signal)) / diff(range(signal))
      
      anatomic[signal > 0.5] <- signal[signal > 0.5]
      
      slices3d(anatomic,col=mri.colors(256,256)$all)
      
    } else {
      stop("misc3d is required.")
    }
  } else {
    partition <- ceiling(sqrt(dim(anatomic)[3]))

    switch (device,
            "png" = png(filename=file, width = 200*partition, height = 200*partition, pointsize=12, bg="transparent", res=NA),
            "jpeg" = jpeg(filename=file, width = 200*partition, height = 200*partition,
              quality = 100, pointsize = 12, bg = "transparent", res=NA),
            "ppm" = bitmap(file,type="ppm",height=2*partition,width=2*partition,res=64,pointsize=12),
            X11())
    
    oldpar <- par(mfrow=c(partition,partition), mar=c(0.25,0.25,0.25,.25), mgp=c(2,1,0))
    
    alim <- range(anatomic,na.rm=TRUE,finite=TRUE)
    zlim <- range(signal,na.rm=TRUE,finite=TRUE)
    for (i in 1:dim(anatomic)[3]) {
      image(anatomic[,,i], xaxt="n", yaxt="n", zlim=alim, col=mri.colors(256,256)$gray)
      if (any(signal[,,i])) image(signal[,,i], zlim=zlim ,col=c(0,mri.colors(256,256)$col), add=TRUE)
      
      if (i == pos[3]) {
        lines(c(0,1),c(pos[2],pos[2])/dim(anatomic)[2],col=2)
        lines(c(pos[1],pos[2])/dim(anatomic)[1],c(0,1),col=2)
      }     
    }
    
    par(oldpar)
    switch (device,
            "png" = dev.off(),
            "jpeg" = dev.off(),
            "ppm" = dev.off())
  }
  invisible(NULL)
}
