findclusters <- function(x,thresh){
   dx <- dim(x)
   tx <- as.integer(x>thresh)
   z <- .Fortran("ccluster",
                 size=as.integer(tx),
                 as.integer(dx[1]),
                 as.integer(dx[2]),
                 as.integer(dx[3]),
                 clusterid=integer(prod(dx)),
                 PACKAGE="fmri")[c("size","clusterid")]
    dim(z$size) <- dx
    dim(z$clusterid) <- dx
    z
}

simclusterthr <- function(n,distr=c("norm","t"),df=100,bw=0,kern="Gaussian",
    nsim=1000,seed=1,qgrid = pnorm(seq(1.4,4,.05)),filename="tempres"){
    require(aws)
    nq <- length(qgrid)
    nsize <- 20
    results <- array(0,c(nq,nsize))
    set.seed(seed)
    dx <- c(n,n,n)
    quant <- if(distr[1]=="t") qt(qgrid,df) else qnorm(qgrid)
    for(i in 1:nsim){
       x <- if(distr[1]=="t") rt(n^3,df) else rnorm(n^3)
       dim(x) <- dx
       if(bw>0) {
          x <- kernsm(x,bw,kern=kern,unit="FWHM")@yhat
          x <- x/sqrt(var(x))
        }
        for(j in 1:nq){
            mc <- min(20,max(findclusters(x,quant[j])$size))
           results[j,1:mc] <- results[j,1:mc]+1
        }
        cat("Time",format(Sys.time()),"after ",i,"simulations:\n")
        if(i%/%100*100==i) print(results[seq(1,nq,5),]/i)
        if(i%/%1000*1000==i){
           tmp <- list(results=results/i,qgrid=qgrid,n=n,distr=distr[1],df=df,bw=bw,kern=kern)
           save(tmp,file=filename)
        }
    }
    list(results=results/nsim,qgrid=qgrid,n=n,distr=distr[1],df=df,bw=bw,kern=kern)
}
getkv0 <- function(param,mpredf=mpredfactor,irho=1,alpha=.05,ncmin=2){
   ns <- 50000
   data <- clnorm64comp$results[,1:19,]
   kgrid <- qnorm(clnorm64comp$qgrid)
   nc <- 2:20
   kv <- (param[1]+param[2]*log((alpha+1/2/ns)/(1-alpha+1/ns)))*mpredf[nc,irho]
   kv
      }

      getalphaclust <- function(alpha,clustertable,calpha=seq(.001,.1,.001),ncmin=2){
        ## get reference alpha for clusterthreshold with clustersizes ncmin:20
         if(ncmin <2) ncmin <- 2
         if(ncmin >20) ncmin <- 20
         cta <- clustertable[,ncmin-1]
         ca0 <- min(calpha[cta<alpha])
         cta0 <- max(cta[cta<alpha])
         ca1 <- max(calpha[cta>alpha])
         cta1 <- min(cta[cta>alpha])
         ca0+ (cta1-alpha)*(ca1-ca0)/(cta1-cta0)
      }

      fmri.cluster <- function(spm, mode="basic", na.rm=FALSE, alpha=.05, ncmin=2, minimum.signal=0){
        args <- sys.call()
        args <- c(spm$call,args)
      cat("fmri.pvalue: entering function\n")

      if (!("fmrispm" %in% class(spm)) ) {
        warning("fmri.cluster: data not of class <fmrispm>. Try to proceed but strange things may happen")
      }

      if (!is.null(attr(spm, "smooth"))) {
        if (!is.null(attr(spm, "residuals"))) {
          type <- "t"
          df <- spm$df
          if(is.null(df)) df <- abs(diff(dim(attr(spm, "design"))))
        } else {
          type <- "norm"
          df <- 1000 # this is actually not needed, placeholder
        }
      } else {
        type <- "t"
        df <- spm$df
      }
        corr <- mean(spm$scorr)
        if(is.null(corr)) corr <- 0
        stat <- (spm$cbeta-minimum.signal)/sqrt(spm$var)
        dim(stat) <- prod(spm$dim[1:3])

      #  clustersizes to use
        clusters <- ncmin:20
        data(clustertable)
        alphaclust <- getalphaclust(alpha,clustertable,ncmin)
        load(system.file("extdata/cluster.rda",package="fmri"))
        irho <- as.integer(corr/0.05)+1
        if(irho>13) {
           error("to much spatial correlation")
        }
        # correct for size of multiplicity
        n2 <- sum(spm$mask)
        n1 <- 64^3
        alpha <- 1-(1-alpha)^(n2/n1)
        # this reflects the use of n2 instead of 64^3 voxel (simulation)
        # get critical values
        kv <- getkv0(param,mpredf=mpredfactor,irho=irho,alpha=alpha,ncmin=2)
      # this gives a vector of kritical values corresponding to cluster sizes
      # now adjust for distribution
        if(type=="t") kv <- qf(pnorm(kv),df)
        detected <- array(0,dim(spm$dim))
        for(ic in 1:length(clusters)){
           ttt <- findclusters(stat,kv[ic])
           detected[ttt>clusters[ic]] <- 1
        }
        detected <- detected*spm$mask
        detected
      }
