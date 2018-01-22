searchlight <- function(radius){
   rad <- as.integer(radius)
   nr <- 2*rad+1
   indices <- rbind(rep(-rad:rad,nr*nr),
              rep(-rad:rad,rep(nr,nr)),
              rep(-rad:rad,rep(nr*nr,nr)))
   indices[,apply(indices^2,2,sum)<=radius^2]
}

searchlightdistr <- function(df,nregion,kind=0,nsim=10000){
rn <- matrix(rt(10000*nregion,df),nregion,10000)
rsl <- apply(if(kind==0) abs(rn) else rn^2, 2, mean)
list(p=seq(1/(2*nsim),(2*nsim-1)/(2*nsim),length=nsim),kvalue=sort(rsl),df=df,kind=kind,nsim=nsim)
}

fmri.searchlight <- function(spm, alpha=.05, radius, minimum.signal=0, kind=c("abs","squared")){
  args <- sys.call()
  args <- c(spm$call,args)

cat("fmri.serchlight: entering function\n")

if (!("fmrispm" %in% class(spm)) ) {
  warning("fmri.searchlight: data not of class <fmrispm>. Try to proceed but strange things may happen")
}
stat <- spm$beta
stat[stat>0] <- pmax(0,stat[stat>0]-minimum.signal)
stat[stat<0] <- pmin(0,stat[stat<0]+minimum.signal)
stat <- stat/sqrt(spm$var)
dim(stat) <- dimspm <- spm$dim[1:3]
sregion <- searchlight(radius)
nregion <- dim(sregion,2)
kind <- if(kind=="abs") 0 else 1
stat <- if(kind==0) abs(stat) else stat^2
##
##  get statistics over searchlights
##
stat <- .Fortran("slight",
                as.double(stat),
                as.logical(spm$mask),
                as.integer(dimspm[1]),
                as.integer(dimspm[2]),
                as.integer(dimspm[3]),
                as.integer(sregion),
                as.integer(nregion),
                stat=double(prod(dimspm)),
                PACKAGE="fmri")$stat
##
##  approximate distribution of searchlight statistics
##  assumes t-distributed spm
##
distr <- searchlightdistr(spm$df,nregion,kind)
nvoxel <- sum(mask)
pv <- array(1,dimspm)
pv[mask] <- .Fortran("getslpv",
               as.double(stat[mask]),
               as.integer(nvoxel),
               as.double(distr$p),
               as.double(distr$kvalue),
               as.integer(distr$nsim),
               pvalue = double(nvoxel),
               PACKAGE="fmri")$pvalue
ind <- fdr(pv[mask], alpha)
thresh <- min(stat[mask][ind])
mask[stat < thresh] <- FALSE
    mask <- mask & spm$mask
    pv[!mask] <- 1
    z <- list(pvalue = pv, weights = spm$weights, dim = spm$dim,
            hrf = spm$hrf)
        class(z) <- c("fmridata", "fmripvalue")
        z$roixa <- spm$roixa
        z$roixe <- spm$roixe
        z$roiya <- spm$roiya
        z$roiye <- spm$roiye
        z$roiza <- spm$roiza
        z$roize <- spm$roize
        z$roit <- spm$roit
        z$header <- spm$header
        z$format <- spm$format
        z$dim0 <- spm$dim0
        z$args <- z$args
        attr(z, "file") <- attr(spm, "file")
        attr(z, "white") <- attr(spm, "white")
        attr(z, "design") <- attr(spm, "design")
        attr(z, "smooth") <- "Not smoothed"
        attr(z, "mode") <- paste("Threshold: searchlight radius=",radius, "\n")
        z
    }
