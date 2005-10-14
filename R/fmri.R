perform.aws <- function(beta,variance,hmax=4,hinit=1,weights=c(1,1,1),qlambda=1,lkern="Gaussian",scorr=0,...) {
  require(aws)
  variance[variance < quantile(variance,0.25)] <- quantile(variance,0.25)

  ttthat <- aws(beta, sigma2=variance, hmax=hmax, hinit=hinit,
                 qlambda=qlambda, qtau=1,wghts=weights,lkern=lkern,scorr=scorr,...)

  z <- list(hat=ttthat$theta, var=ttthat$var, ni=ttthat$ni,
  hmax=ttthat$hmax, lseq=ttthat$lseq, y = ttthat$y, mae=ttthat$mae)
  z
}

plot.pvalue <- function(stat, anatomic,rx,ry,rz, pvalue, x=-1, y=-1, z=-1, zlim=0, device="X11", file="plot.png") {
  alim <- range(anatomic)

  anatomic[anatomic<0] <- 0

  i <- dim(stat)[1]
  j <- dim(stat)[2]
  k <- dim(stat)[3]
  p <- pvalue(stat,i,j,k,rx,ry,rz)
  mask <- array(1,dim=dim(stat)[1:3])
  maximum <- threshold(max(p),i,j,k,rx,ry,rz)
  mask[stat<maximum] <- 0  
  p[!mask] <- 1
  p[p>pvalue] <- 1

  signal <- -log(p)
  signal[mask==0] <- 0
#  zlim <- max(zlim,signal)
  
  switch (device,
          "png" = png(filename=file, width = 1000, height = 1000, pointsize=12, bg="transparent", res=NA),
          "jpeg" = jpeg(filename=file, width = 1000, height = 1000,
            quality =100, pointsize=12, bg="transparent", res=NA),
          "ppm" = bitmap(file,type="ppm",height=10,width=10,res=64,pointsize=12),
          X11())

  partition <- as.integer(sqrt(dim(anatomic)[3])) + 1
  oldpar <- par(mfrow=c(partition,partition), mar=c(0,0,0,.25), mgp=c(2,1,0))
  
  for (i in 1:dim(anatomic)[3]) {
    image(anatomic[,,i], xaxt="n", yaxt="n", zlim=alim, col=grey(1:255/255))
    if (any(signal[,,i]))
#      image(signal[,,i], zlim=c(0,zlim) ,col=c(0,heat.colors(512)), add=TRUE)
      image(signal[,,i], col=c(0,heat.colors(512)), add=TRUE)
    if (i == z) {
      lines(c(0,1),c(y,y)/dim(anatomic)[2],col=2)
      lines(c(x,x)/dim(anatomic)[1],c(0,1),col=2)
    }     
  }

  image(t(matrix(seq(0,zlim,length=512),512,20)),col=c(0,heat.colors(512)))

  lines(c(0,1),rep(-log(0.005),2)/zlim,col=1)
  text(0,-log(0.005)/zlim,"0.005")

  lines(c(0,1),rep(-log(0.01),2)/zlim,col=1)
  text(0,-log(0.01)/zlim,"0.01")

  lines(c(0,1),rep(-log(0.025),2)/zlim,col=1)
  text(0,-log(0.025)/zlim,"0.025")

  lines(c(0,1),rep(-log(0.05),2)/zlim,col=1)
  text(0,-log(0.05)/zlim,"0.05")

  lines(c(0,1),rep(-log(0.1),2)/zlim,col=1)
  text(0,-log(0.1)/zlim,"0.1")

  lines(c(0,1),rep(-log(0.2),2)/zlim,col=1)
  text(0,-log(0.2)/zlim,"0.2")

  lines(c(0,1),rep(-log(0.4),2)/zlim,col=1)
  text(0,-log(0.4)/zlim,"0.4")
  
  par(oldpar)
  as.numeric(dev.cur())
}

plot.fmri <- function(signal, mask, anatomic, x=-1, y=-1, z=-1, zlim=0, device="X11", file="plot.png") {
  alim <- range(anatomic)
  zlim <- max(zlim,signal)

  signal[mask==0] <- 0
  anatomic[anatomic<0] <- 0
  
  if (device == "png") {
    png(filename=file, width = 480, height = 480, pointsize=12,
    bg="transparent", res=NA)
  } else {
    X11()
  }

  partition <- as.integer(sqrt(dim(anatomic)[3])) + 1
  oldpar <- par(mfrow=c(partition,partition), mar=c(0,0,0,.25), mgp=c(2,1,0))
  
  for (i in 1:dim(anatomic)[3]) {
    image(anatomic[,,i], xaxt="n", yaxt="n", zlim=alim, col=grey(1:255/255))
    if (any(signal[,,i]))
      image(signal[,,i], zlim=c(0,zlim) ,col=c(0,heat.colors(512)), add=TRUE)
    if (i == z) {
      lines(c(0,1),c(y,y)/dim(anatomic)[2],col=2)
      lines(c(x,x)/dim(anatomic)[1],c(0,1),col=2)
    }     
  }
  
  par(oldpar)
  as.numeric(dev.cur())
}

#
# plot the (logical) difference between mask1 and mask2 as three
# color overlay over anatomic image:
# red for mask1 && mask2
# green for mask1 && !mask2
# blue for !mask1 && mask2
# mark position pos with lines
#
plot.diff <- function(mask1, mask2, anatomic = array(0,dim=dim(mask1)), pos = c(-1,-1,-1), device="X11", file="plot.png") {
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
plot.signal <- function(signal, anatomic = array(0,dim=dim(signal)), pos = c(-1,-1,-1), device="X11", file="plot.png") {
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
plot.timeseries <- function(fmri, fit, sigvec, tvec, thr, origin, device="X11", file="plot.png") {
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
