read.ANALYZE <- function(prefix = "", picstart = "", numbpic = 1) {
  if (require(AnalyzeFMRI)) {

    ttt <- f.read.analyze.volume(paste(prefix, picstart, ".img", sep=""));
    dt <- dim(ttt)
    cat(".")
    header <- f.read.analyze.header(paste(prefix, picstart, ".img", sep=""));

    if (numbpic > 1) { 
      for (i in (picstart+1):(picstart+numbpic-1)) {
        a <- f.read.analyze.volume(paste(prefix, i, ".img", sep=""))
        if (sum() != 0)
          cat("Error: wrong spatial dimension in picture",i)
        ttt <- c(ttt,a);
        dt[4] <- dt[4] + dim(a)[4]
        cat(".")
      }
    }

    cat("\n")
    dim(ttt) <- c(dt)
    list(ttt=ttt,header=header)
  } else {
    cat("Error: library AnalyzeFMRI not found\n")
    list(ttt=NA,header=NA)
  }
}

read.AFNI <- function(file) {
  conhead <- file(paste(file,".HEAD",sep=""),"r")
  header <- readLines(conhead)
  close(conhead)

  types <- NULL
  args <- NULL
  counts <- NULL
  values <- NULL
  
  for (i in 1:length(header)) {
    if (regexpr("^type *= *", header[i]) == 1) {
      tmptype <- strsplit(header[i]," *= *")[[1]][2]
      types <- c(types,tmptype)
      args <- c(args,strsplit(header[i+1]," *= *")[[1]][2])
      tmpcounts <- as.numeric(strsplit(header[i+2]," *= *")[[1]][2])
      counts <- c(counts,tmpcounts)
      i <- i+3
      tmpvalue <- ""
      while ((regexpr("^$", header[i]) != 1) && (i <= length(header))) {
        tmpvalue <- paste(tmpvalue,header[i])
        i <- i+1
      }
      tmpvalue <- sub("^ +","",tmpvalue)
      if ((tmptype == "integer-attribute") || (tmptype == "float-attribute")) {
        tmpvalue <- as.numeric(strsplit(tmpvalue," +")[[1]])
      }
      values <- c(values,list(value=tmpvalue))
    }        
  }

  names(values) <- args

  dx <- values$DATASET_DIMENSIONS[1]
  dy <- values$DATASET_DIMENSIONS[2]
  dz <- values$DATASET_DIMENSIONS[3]
  dt <- values$DATASET_RANK[2]
  size <- file.info(paste(file,".BRIK",sep=""))$size/(dx*dy*dz*dt)

  if (as.integer(size) == size) {
    conbrik <- file(paste(file,".BRIK",sep=""),"rb")
    myttt<- readBin(conbrik, "int", n=dx*dy*dz*dt*size, size=size, signed=FALSE, endian="big")
    close(conbrik)
    dim(myttt) <- c(dx,dy,dz,dt)
    list(ttt=myttt,header=values)
  } else {
    cat("Error reading file: Could not detect size per voxel\n")
    list(ttt=NA,header=values)    
  }
}

write.AFNI <- function(file, ttt, label, note="", origin=c(0,0,0), delta=c(4,4,4), idcode="WIAS_noid") {
  conhead <- file(paste(file, ".HEAD", sep=""), "w")
  writeChar(AFNIheaderpart("string-attribute","HISTORY_NOTE",note),conhead,eos=NULL)
  writeChar(AFNIheaderpart("string-attribute","TYPESTRING","3DIM_HEAD_ANAT"),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","IDCODE_STRING",idcode),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","IDCODE_DATE",date()),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","SCENE_DATA",c(0,2,1,-999,-999,-999,-999,-999)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","ORIENT_SPECIFIC",c(0,3,4)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("float-attribute","ORIGIN",origin),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("float-attribute","DELTA",delta),conhead,eos=NULL)  
  minmax <- function(y) {r <- NULL;for (k in 1:dim(y)[4]) {r <- c(r,min(y[,,,k]),max(y[,,,k]))}; r}
  writeChar(AFNIheaderpart("float-attribute","BRICK_STATS",minmax(ttt)),conhead,eos=NULL)
  writeChar(AFNIheaderpart("integer-attribute","DATASET_RANK",c(3,dim(ttt)[4],0,0,0,0,0,0)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","DATASET_DIMENSIONS",c(dim(ttt)[1:3],0,0)),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("integer-attribute","BRICK_TYPES",rep(1,dim(ttt)[4])),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("float-attribute","BRICK_FLOAT_FACS",rep(0.0,dim(ttt)[4])),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","BRICK_LABS",paste(label,collapse="~")),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","BRICK_KEYWORDS",paste(rep("",length(label)),collapse="~")),conhead,eos=NULL)  
  writeChar(AFNIheaderpart("string-attribute","BYTEORDER_STRING","MSB_FIRST"),conhead,eos=NULL)  
  close(conhead)

  conbrik <- file(paste(file, ".BRIK", sep=""), "wb")
  dim(ttt) <- NULL
  writeBin(as.integer(ttt), conbrik,size=2)
  close(conbrik)
}

AFNIheaderpart <- function(type, name, value) {
  a <- "\n"
  a <- paste(a, "type = ", type, "\n", sep="")
  a <- paste(a, "name = ", name, "\n", sep="")
  if (regexpr("string",type) == 1) {
    value <- paste("'", value, "~", sep="")
    a <- paste(a, "count = ", nchar(value) - 1, "\n", sep ="")
    a <- paste(a, value, "\n", sep="")
  } else {
    a <- paste(a, "count = ", length(value), "\n", sep ="")
    j <- 0
    while (j<length(value)) {
      left <- length(value) - j
      if (left>4) left <- 5
      a <- paste(a, paste(value[(j+1):(j+left)],collapse="  "), "\n", sep="  ")
      j <- j+5
    }
  }
  a
}

write.ANALYZE <- function(ttt, name = "data", size="int", voxelsize = c(2,2,2)) {
  if (require(AnalyzeFMRI)) {
    f.write.analyze(ttt, file=name, size=size, voxelsize)
  } else {
    cat("Error: library AnalyzeFMRI not found\n")
    NA
  }
}

create.stimulus <- function(scans=1 ,onsets=c(1) ,length=1, rt=3, mean=TRUE) {
  numberofonsets <- length(onsets)
  stimulus <- rep(0, scans)
  
  for (i in 1:numberofonsets) {
    onset <- onsets[i]
    for (j in onset:(onset+length-1)) {
      stimulus[j] <- 1
    }
  }

  mygamma <- function(x, a1, a2, b1, b2, c) {
    d1 <- a1 * b1
    d2 <- a2 * b2
    c1 <- ( x/d1 )^a1
    c2 <- c * ( x/d2 )^a2
    res <- c1 * exp(-(x-d1)/b1) - c2 * exp(-(x-d2)/b2)
    res
  }

  hrf <- convolve(stimulus,mygamma(scans:1, 6, 12, 0.9/rt, 0.9/rt, 0.35))
  
  if (mean) {
    hrf - mean(hrf)
  } else {
    hrf
  }
}

create.designmatrix <- function(hrf, order=0) {
  stimuli <- dim(hrf)[2]
  scans <- dim(hrf)[1]

  z <- matrix(0, scans, stimuli+order+1)

  for (i in 1:stimuli) {
    z[,i] <- hrf[,i]
  }

  z[,stimuli+1] <- 1

  if (order != 0) {
    for (i in (stimuli+2):(stimuli+order+1)) {
      z[,i] <- (1:scans)^(i-stimuli-1)
    }
  }
  
  z
}

create.arcorrection <- function(scans, rho=0) {
  a <- array(0,dim=c(scans,scans))
  
  a[1,1] <- 1

  rho0 <- sqrt(1-rho^2)
  rho1 <- -rho/rho0
  rho2 <- 1/rho0
  
  for (i in 2:scans) {
    a[i,i-1] <- rho1
    a[i,i] <- rho2
  }

  a
}


calculate.lm <- function(ttt,z,ar,vtype="var") {
  zprime <- a %*% z
  svdresult <- svd(zprime)
  u <- svdresult$u
  v <- svdresult$v
  vt <- t(v)
  lambda1 <- diag(1/svdresult$d)
  lambda2 <- diag(1/svdresult$d^2)

  sigma <- v %*% lambda2 %*% vt

  dy <- dim(ttt)
  variance <- array(0,dim=c(dy[1:3]))
  dim(ttt) <- c(prod(dy[1:3]),dy[4])
  beta <- ttt %*% t(a) %*% u %*% lambda1 %*% vt
  residuals <- ttt - beta %*% t(z)
  b <- rep(1/dim(z)[1],dim(z)[1])
  if (vtype == "var") {
    variance <- (residuals^2 %*% b)*dim(z)[1]/(dim(z)[1]-dim(z)[2])*sigma[1,1]
  }
  dim(beta)<-c(dy[1:3],dim(z)[2])
  dim(variance) <- c(dy[1:3])
  dim(residuals) <- dy
  
  result <- list(beta=beta,var=variance,residuals=residuals)
  result
}

perform.aws <- function(beta,variance,hmax=4,hinit=1,weights=c(1,1,1),vweights=NULL,qlambda=1) {
  require(aws)
  variance[variance < quantile(variance,0.25)] <- quantile(variance,0.25)

  ttthat <- vaws(beta, sigma2=variance, hmax=hmax, hinit=hinit,
                 qlambda=qlambda, qtau=1,wghts=weights,vwghts=vweights)

  z <- list(hat=ttthat$theta, var=ttthat$var)
  z
}

calculate.threshold <- function(var1,var2,alpha=0.95,gamma=0.5) {
  delta <- var1 * qchisq(gamma,1)
  calpha <- qchisq(alpha,1 + delta/var2)
  test <- (var2 + 2 * delta)/(1 + delta/var2) * calpha
  test
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
  zlim <- max(zlim,signal)
  
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
      image(signal[,,i], zlim=c(0,zlim) ,col=c(0,heat.colors(512)), add=TRUE)
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


pvalue <- function(z,i,j,k,rx,ry,rz) {
  rho0 <- 1 - pnorm(z)
  rho1 <- (4 * log(2))^.5 / (2 * pi) * exp(-z*z/2)
  rho2 <- (4 * log(2)) / (2 * pi)^1.5 * exp(-z*z/2) * z
  rho3 <- (4 * log(2))^1.5 / (2 * pi)^2 * exp(-z*z/2) * (z*z -1)
  r0 <- 1
  r1 <- (i-1)*rx + (j-1)*ry +(k-1)*rz
  r2 <- (i-1)*(j-1)*rx*ry + (j-1)*(k-1)*ry*rz +(i-1)*(k-1)*rx*rz
  r3 <- (i-1)*(j-1)*(k-1)*rx*ry*rz

  rho0 * r0 + rho1 * r1 + rho2 * r2 + rho3 * r3
}

pvalueprime <- function(z,i,j,k,rx,ry,rz) {
  rho0prime <- dnorm(z)
  rho1prime <- -(4 * log(2))^.5 / (2 * pi) * exp(-z*z/2) * z
  rho2prime <- -(4 * log(2)) / (2 * pi)^1.5  * exp(-z*z/2) * (z*z -1)
  rho3prime <- -(4 * log(2))^1.5 / (2 * pi)^2 * exp(-z*z/2) * (z*z*z - z)
  r0 <- 1
  r1 <- (i-1)*rx + (j-1)*ry +(k-1)*rz
  r2 <- (i-1)*(j-1)*rx*ry + (j-1)*(k-1)*ry*rz +(i-1)*(k-1)*rx*rz
  r3 <- (i-1)*(j-1)*(k-1)*rx*ry*rz

  rho0prime * r0 + rho1prime * r1 + rho2prime * r2 + rho3prime * r3
}

resel <- function(voxeld, hmax, hv=0.919) hv * voxeld / hmax

threshold <- function(p,i,j,k,rx,ry,rz) {
  f <- function(x,par){
    abs(pvalue(x,par$i,par$j,par$k,par$rx,par$ry,par$rz)-par$p)
  }
  optimize(f=f,lower=2,upper=10,tol=1e-6,par=list(p=p,i=i,j=j,k=k,rx=rx,ry=ry,rz=rz))$minimum
}

gkernsm <- function(y,h=1) {
  grid <- function(d) {
    d0 <- d%/%2+1
    gd <- seq(0,1,length=d0)
    if (2*d0==d+1) gd <- c(gd,-gd[d0:2]) else gd <- c(gd,-gd[(d0-1):2])
    gd
  }
  dy <- dim(y)
  if (is.null(dy)) dy<-length(y)
  ldy <- length(dy)
  if (length(h)!=ldy) h <- rep(h[1],ldy)
  kern <- switch(ldy,dnorm(grid(dy),0,h/dy),
                 outer(dnorm(grid(dy[1]),0,h[1]/dy[1]),
                       dnorm(grid(dy[2]),0,h[2]/dy[2]),"*"),
                 outer(outer(dnorm(grid(dy[1]),0,h[1]/dy[1]),
                             dnorm(grid(dy[2]),0,h[2]/dy[2]),"*"),
                       dnorm(grid(dy[3]),0,h[3]/dy[3]),"*"))
  kern <- kern/sum(kern)
  kernsq <- sum(kern^2)
  list(gkernsm=convolve(y,kern,conj=TRUE),kernsq=kernsq)
}

correlation <- function(res,mask) {
  meanpos <- function(a) mean(a[a!=0])
  varpos <- function(a) var(a[a!=0])

  dr <- dim(res)
  if ( (length(dim(mask)) == length(dr))
      && (sum(dim(mask)[1:3] != dr[1:3]) == 0)
      && (dim(mask)[4] == 1)) {
    mask <- rep(mask,dr[4])
    dim(mask) <- dr
    x <- meanpos(res[-1,,,]*res[-dr[1],,,]*mask[-1,,,]/sqrt(varpos(res[-1,,,]*mask[-1,,,]) * varpos(res[-dr[1],,,]*mask[-dr[1],,,])))
    y <- meanpos(res[,-1,,]*res[,-dr[2],,]*mask[,-1,,]/sqrt(varpos(res[,-1,,]*mask[,-1,,]) * varpos(res[,-dr[2],,]*mask[,-dr[2],,])))
    z <- meanpos(res[,,-1,]*res[,,-dr[3],]*mask[,,-1,]/sqrt(varpos(res[,,-1,]*mask[,,-1,]) * varpos(res[,,-dr[3],]*mask[,,-dr[3],])))
    c(x,y,z)
  } else {
    cat("Error: dimension of mask and residui matrices do not match\n")    
  }
}

smoothness <- function(cor, dim, h = 0, step = 0.1) {
  field <- array(rnorm(prod(dim)),dim)
  repeat {
    h <- h + step
    z <- gkernsm(field, h)$gkernsm;
    corg <- mean(z[,-1,]*z[,-dim[2],]/sqrt(var(z[,-1,]) * var(z[,-dim[2],])))
#    cat(corg,h,"\n")
    if (corg > cor) break
  }
  h
}

bandwidth <- function(res,mask) { # second argument !!!!
  cxyz <- correlation(res,mask)
  bwx <- smoothness(cxyz[1],dim(res)[1:3],0.5,0.01)
  bwy <- smoothness(cxyz[2],dim(res)[1:3],0.5,0.01)
  bwz <- smoothness(cxyz[3],dim(res)[1:3],0.5,0.01)
  meingauss<-gkernsm(array(rnorm(prod(dim(res)[1:3])),dim=dim(res)[1:3]),c(bwx,bwy,bwz))
  list(bwcorr=1/meingauss$kernsq,bw=c(bwx,bwy,bwz))
}
