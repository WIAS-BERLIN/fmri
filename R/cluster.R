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
       if(bw>0) x <- kernsm(x,bw,kern=kern,unit="FWHM")@yhat
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
