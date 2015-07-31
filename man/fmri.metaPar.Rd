\name{fmri.metaPar}
\alias{fmri.metaPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fmri.metaPar(G, V, X = NULL, model = NULL, method = "REML", weighted = TRUE, mask = NULL, cluster = 2, white = 1, wghts = c(1, 1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
%%     ~~Describe \code{G} here~~
}
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{weighted}{
%%     ~~Describe \code{weighted} here~~
}
  \item{mask}{
%%     ~~Describe \code{mask} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
}
  \item{white}{
%%     ~~Describe \code{white} here~~
}
  \item{wghts}{
%%     ~~Describe \code{wghts} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (G, V, X = NULL, model = NULL, method = "REML", weighted = TRUE, 
    mask = NULL, cluster = 2, white = 1, wghts = c(1, 1, 1)) 
{
    t1 <- Sys.time()
    cat("fmri.meta: entering function\n")
    if (length(dim(G)) != 4) 
        stop("fmri.meta: The first argument does not seem to be a 4D-Array.\n")
    if (length(dim(V)) != 4) 
        stop("fmri.meta: The second argument does not seem to be a 4D-Array.\n")
    if (!identical(dim(G), dim(V))) 
        stop("fmri.meta: dimensionality of the two arrays are not the same.\n")
    if ((!is.null(X)) && (!is.data.frame(X))) 
        stop("fmri.meta: The design matrix X is not a data frame.")
    if (is.null(mask)) {
        mask <- array(TRUE, dim = dim(G)[1:3])
    }
    if (!identical(dim(mask), dim(G)[1:3])) 
        stop("fmri.meta: mask dimensionality does not match the data.\n")
    if (sum(mask) == 0) 
        stop("fmri.meta: your head mask is empty, nothing to do.\n")
    dim.spm <- dim(G)[1:3]
    dx <- dim(G)[1]
    dy <- dim(G)[2]
    dz <- dim(G)[3]
    dk <- dim(G)[4]
    total0 <- dx * dy * dz
    if (!is.null(attr(G, "xind"))) {
        xind <- attr(G, "xind")
    }
    else {
        xind <- NULL
    }
    if (!is.null(attr(G, "yind"))) {
        yind <- attr(G, "yind")
    }
    else {
        yind <- NULL
    }
    if (!is.null(attr(G, "zind"))) {
        zind <- attr(G, "zind")
    }
    else {
        zind <- NULL
    }
    dim(G) <- c(total0, dk)
    dim(V) <- c(total0, dk)
    total <- sum(mask)
    GV <- array(0, dim = c(total, 2 * dk))
    GV[, 1:dk] <- G[mask, ]
    GV[, (dk + 1):(2 * dk)] <- V[mask, ]
    rm(G, V)
    gc()
    cbeta <- array(0, dim = dim.spm)
    var <- array(1, dim = dim.spm)
    resid <- array(1, dim = c(dk, total0))
    if (is.null(X)) {
        mods <- NULL
    }
    else if (is.null(model)) {
        mods <- NULL
    }
    else {
        mods = model
    }
    funcD <- function(y, mods, X, method, weighted, dk, dot = FALSE) {
        g <- y[1:dk]
        v <- y[(dk + 1):(2 * dk)]
        if (dot) 
            cat(".")
        fit <- rma.uni(g, v, mods = mods, data = X, method = method, 
            weighted = weighted, control = list(stepadj = 0.5, 
                threshold = 0.001, maxiter = 1000))
        cbeta <- fit$b[1]
        var <- fit$vb[1, 1]
        resid <- g - cbeta
        df <- dk - dim(fit$X)[2]
        result <- c(cbeta, var, df, resid)
    }
    cat("fmri.meta: calculating MEMA model\n")
    if (cluster == 1) {
        calc.time <- round(total * 1/11400, 0)
        cat(paste("Time of calculation takes about", calc.time, 
            "min.\n", collapse = " "))
        param <- apply(GV, 1, funcD, mods, X, method, weighted, 
            dk, dot = FALSE)
    }
    if (cluster > 1) {
        cl <- makeCluster(cluster)
        clusterEvalQ(cl, c(library(metafor), "funcD"))
        param <- parApply(cl, GV, 1, funcD, mods, X, method, 
            weighted, dk)
        stopCluster(cl)
    }
    cbeta[mask] <- param[1, ]
    var[mask] <- param[2, ]
    df <- param[3, 1]
    resid[, mask] <- param[4:(dk + 3), ]
    cat("fmri.meta: finished\n")
    cat("fmri.meta: calculating spatial correlation\n")
    lags <- c(5, 5, 3)
    corr <- .Fortran("mcorr", as.double(resid), as.logical(mask), 
        as.integer(dx), as.integer(dy), as.integer(dz), as.integer(dk), 
        scorr = double(prod(lags)), as.integer(lags[1]), as.integer(lags[2]), 
        as.integer(lags[3]), PACKAGE = "fmri", DUP = FALSE)$scorr
    dim(corr) <- lags
    qscale <- range(resid)
    scale <- max(abs(qscale))/32767
    resid <- writeBin(as.integer(resid/scale), raw(), 2)
    if (sum(!is.finite(corr))) {
        warning("fmri.meta:  infinite spatial correlations were found.\n              The analyzed region maybe too small. No bandwidths determinable.")
        bw <- NULL
        rxyz <- NULL
    }
    else {
        bw <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B", 
            lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 6), 
            lag = lags, data = corr)$par
        bw[bw <= 0.25] <- 0
        dim(bw) <- c(1, 3)
        rxyz <- c(resel(1, bw[1]), resel(1, bw[2]), resel(1, 
            bw[3]))
        dim(rxyz) <- c(1, 3)
    }
    cat("fmri.meta: finished\n")
    cat("fmri.meta: exiting function\n")
    result <- list(cbeta = cbeta, var = var, mask = mask, res = resid, 
        resscale = scale, rxyz = rxyz, scorr = corr, weights = wghts, 
        dim = c(dim.spm, dk), bw = bw, df = df)
    if (!is.null(xind)) {
        result$roixa <- min(xind)
        result$roixe <- max(xind)
    }
    if (!is.null(yind)) {
        result$roiya <- min(yind)
        result$roiye <- max(yind)
    }
    if (!is.null(zind)) {
        result$roiza <- min(zind)
        result$roize <- max(zind)
    }
    result$subjects <- dk
    result$method <- method
    result$weighted <- weighted
    if (is.null(model)) {
        result$model <- as.character("~ 1")
    }
    else {
        result$model <- paste(c("", as.character(model)), collapse = " ")
    }
    result$cluster <- cluster
    class(result) <- c("fmridata", "fmrispm")
    intrcpt <- rep(1, dk)
    attr(result, "design") <- data.frame(cbind(intrcpt, X))
    attr(result$cbeta, "note") <- "intercept or first coefficient"
    attr(result, "estimator") <- "two-stage process"
    attr(result, "white") <- white
    attr(result, "residuals") <- !is.null(scale)
    t2 <- Sys.time()
    cat("time elapsed", difftime(t2, t1, units = "mins"), "mins\n")
    invisible(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
