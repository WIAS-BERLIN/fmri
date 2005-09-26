perform.aws <- function(beta,variance,hmax=4,hinit=1,weights=c(1,1,1),vweights=NULL,qlambda=1,lkern="Gaussian") {
  require(aws)
  variance[variance < quantile(variance,0.25)] <- quantile(variance,0.25)

  ttthat <- aws(beta, sigma2=variance, hmax=hmax, hinit=hinit,
                 qlambda=qlambda, qtau=1,wghts=weights,lkern=lkern)

  z <- list(hat=ttthat$theta, var=ttthat$var)
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


