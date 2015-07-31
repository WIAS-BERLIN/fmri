\name{fmri.lmePar}
\alias{fmri.lmePar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fmri.lmePar(bold, z, fixed = NULL, random = NULL, mask = NULL, ac = 0.3, vtype = "individual", cluster = 2, wghts = c(1, 1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bold}{
%%     ~~Describe \code{bold} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{fixed}{
%%     ~~Describe \code{fixed} here~~
}
  \item{random}{
%%     ~~Describe \code{random} here~~
}
  \item{mask}{
%%     ~~Describe \code{mask} here~~
}
  \item{ac}{
%%     ~~Describe \code{ac} here~~
}
  \item{vtype}{
%%     ~~Describe \code{vtype} here~~
}
  \item{cluster}{
%%     ~~Describe \code{cluster} here~~
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
function (bold, z, fixed = NULL, random = NULL, mask = NULL, 
    ac = 0.3, vtype = "individual", cluster = 2, wghts = c(1, 
        1, 1)) 
{
    t1 <- Sys.time()
    cat("fmri.lme: entering function\n")
    dim.spm <- dim(bold)[1:3]
    dx <- dim(bold)[1]
    dy <- dim(bold)[2]
    dz <- dim(bold)[3]
    dNt <- dim(bold)[4]
    total0 <- dx * dy * dz
    if (length(dim(bold)) != 4) 
        stop("fmri.lme: Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
    if (!is.data.frame(z)) 
        stop("fmri.lme: Something is wrong with the design matrix.")
    if (dim(bold)[4] != dim(z)[1]) 
        stop("fmri.lme: design matrix dimensionality does not match functional data.")
    if ((!length(ac) == 1) && (!identical(dim(ac), dim.spm))) 
        stop("fmri.lme: AR1-array dimensionality does not match functional data.")
    if (is.null(mask)) {
        mask <- array(TRUE, dim = dim.spm)
    }
    if (!identical(dim(mask), dim.spm)) 
        stop("fmri.lme: mask dimensionality does not match functional data.")
    if (sum(mask) == 0) 
        stop("fmri.meta: your head mask is empty, nothing to do.")
    if (!is.null(attr(bold, "xind"))) {
        xind <- attr(bold, "xind")
    }
    else {
        xind <- NULL
    }
    if (!is.null(attr(bold, "yind"))) {
        yind <- attr(bold, "yind")
    }
    else {
        yind <- NULL
    }
    if (!is.null(attr(bold, "zind"))) {
        zind <- attr(bold, "zind")
    }
    else {
        zind <- NULL
    }
    if (length(ac) == 1) {
        white <- 3
    }
    else {
        white <- 1
    }
    dim(bold) <- c(total0, dNt)
    total <- sum(mask)
    if ((cluster == 1) && (total > 26000)) {
        stop("Time for this calculation takes longer than 5 hours. Please use parallelizing option.")
    }
    else if ((cluster == 1) && (total > 16000)) {
        stop("Time for this calculation takes longer than 3 hours. Please use parallelizing option.")
    }
    bold2 <- array(0, dim = c(total, dNt + 1))
    bold2[, 1] <- if (length(ac) == total0) 
        ac[mask]
    else rep(ac[1], total)
    bold2[, -1] <- bold[mask, ]
    rm(ac, bold)
    gc()
    dim <- c(dim.spm, dNt)
    cbeta <- array(0, dim = dim.spm)
    var <- array(1, dim = dim.spm)
    resid <- array(1, dim = c(dNt, total0))
    arfactor <- array(0, dim = dim.spm)
    df <- 100
    dim2 <- NULL
    cbeta2 <- NULL
    var2 <- NULL
    resid2 <- NULL
    df2 <- NULL
    gr <- nlevels(z$group)
    runs <- nlevels(z$run)
    if (gr == 1) {
        if (is.null(fixed)) {
            fixed <- ~0 + hrf + session + drift1:session + drift2:session
        }
        fe.model <- formula(paste(c("y", as.character(fixed)), 
            collapse = " "))
        if (is.null(random)) {
            if (runs == 1) {
                random <- ~0 + hrf | subj
                ar1 <- ~1 | subj/session
            }
            if (runs > 1) {
                if (subj == 1) {
                  random <- ~0 + hrf | session
                  ar1 <- ~1 | session
                }
                else {
                  random <- ~0 + hrf | subj/session
                  ar1 <- ~1 | subj/session
                }
            }
        }
        else if (subj == 1) {
            ar1 <- ~1 | session
        }
        else {
            ar1 <- ~1 | subj/session
        }
        if (vtype == "individual") {
            weights <- varIdent(form = ~1 | subj)
        }
        if (vtype == "equal") {
            weights <- NULL
        }
        funcA <- function(y, z, fe.model, random, ar1, weights, 
            dot = FALSE) {
            z$y <- y[-1]
            if (dot) 
                cat(".")
            lmm <- lme(fixed = fe.model, random = random, correlation = corAR1(value = y[1], 
                form = ar1, fixed = TRUE), weights = weights, 
                method = "REML", control = lmeControl(rel.tol = 1e-06, 
                  returnObject = TRUE), data = z)
            cbeta <- fixef(lmm)[1]
            var <- vcov(lmm)[1, 1]
            resid <- resid(lmm)
            arfactor <- y[1]
            df <- lmm$fixDF$X[1]
            result <- c(cbeta, var, arfactor, df, resid)
        }
        cat("fmri.lme: calculating MFX model\n")
        if (cluster == 1) {
            if (total > 900) {
                calc.time <- round(total * 23/2000, 0)
                warning(paste("Time of calculation takes about", 
                  calc.time, "min, better use parallelizing option.", 
                  collapse = " "))
            }
            param <- apply(bold2, 1, funcA, z, fe.model, random, 
                ar1, weights, dot = TRUE)
        }
        if (cluster > 1) {
            cl <- makeCluster(cluster)
            clusterEvalQ(cl, c(library(nlme), "funcA"))
            param <- parApply(cl, bold2, 1, funcA, z, fe.model, 
                random, ar1, weights)
            stopCluster(cl)
        }
        cbeta[mask] <- param[1, ]
        var[mask] <- param[2, ]
        arfactor[mask] <- param[3, ]
        df <- param[4, 1]
        resid[, mask] <- param[5:(dNt + 4), ]
        random.model <- paste(c("", as.character(random)), collapse = " ")
        cat("fmri.lme: finished\n")
        cat("fmri.lme: calculating spatial correlation\n")
        lags <- c(5, 5, 3)
        corr <- .Fortran("mcorr", as.double(resid), as.logical(mask), 
            as.integer(dx), as.integer(dy), as.integer(dz), as.integer(dNt), 
            scorr = double(prod(lags)), as.integer(lags[1]), 
            as.integer(lags[2]), as.integer(lags[3]), PACKAGE = "fmri", 
            DUP = FALSE)$scorr
        dim(corr) <- lags
        qscale <- range(resid)
        scale <- max(abs(qscale))/32767
        resid <- writeBin(as.integer(resid/scale), raw(), 2)
        if (sum(!is.finite(corr))) {
            warning("fmri.lme: infinite spatial correlations were found. \n              The analyzed region maybe too small. No bandwidths determinable.")
            bw <- NULL
            rxyz <- NULL
        }
        else {
            bw <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B", 
                lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 
                  6), lag = lags, data = corr)$par
            bw[bw <= 0.25] <- 0
            dim(bw) <- c(1, 3)
            rxyz <- c(resel(1, bw[1]), resel(1, bw[2]), resel(1, 
                bw[3]))
            dim(rxyz) <- c(1, 3)
        }
        scale2 <- NULL
        corr2 <- NULL
        bw2 <- NULL
        rxyz2 <- NULL
        dim2 <- NULL
        cat("fmri.lme: finished\n")
    }
    if (gr == 2) {
        cbeta2 <- cbeta
        var2 <- var
        gr.label <- levels(z$group)
        dNt1 <- length(z$group[z$group == gr.label[1]])
        dNt2 <- length(z$group[z$group == gr.label[2]])
        resid <- array(1, dim = c(dNt1, total0))
        resid2 <- array(1, dim = c(dNt2, total0))
        df2 <- df
        if (is.null(fixed)) {
            fixed <- ~0 + hrf:group + session + drift1:session + 
                drift2:session
        }
        fe.model <- formula(paste(c("y", as.character(fixed)), 
            collapse = " "))
        if (is.null(random)) {
            if (runs == 1) {
                random <- list(subj = pdDiag(~0 + hrf:group))
            }
            else if (runs > 1) {
                stop("fmri.lme: Sorry, repeated measures for two groups are not predefined.")
            }
        }
        if (vtype == "individual") {
            weights <- varIdent(form = ~1 | subj)
        }
        if (vtype == "equal") {
            weights <- NULL
        }
        cbeta.label <- c(paste(c("hrf:group", gr.label[1]), collapse = ""), 
            paste(c("hrf:group", gr.label[2]), collapse = ""))
        funcB <- function(y, z, fe.model, random, weights, cbeta.label, 
            dot = FALSE) {
            z$y <- y[-1]
            if (dot) 
                cat(".")
            lmm <- lme(fixed = fe.model, random = random, correlation = corAR1(value = y[1], 
                form = ~1 | subj/session, fixed = TRUE), weights = weights, 
                method = "REML", control = lmeControl(rel.tol = 1e-06, 
                  returnObject = TRUE), data = z)
            cbeta1 <- fixef(lmm)[cbeta.label[1]]
            cbeta2 <- fixef(lmm)[cbeta.label[2]]
            var1 <- vcov(lmm)[cbeta.label[1], cbeta.label[1]]
            var2 <- vcov(lmm)[cbeta.label[2], cbeta.label[2]]
            resid1 <- resid(lmm)[z$group == levels(z$group)[1]]
            resid2 <- resid(lmm)[z$group == levels(z$group)[2]]
            arfactor <- y[1]
            df1 <- lmm$fixDF$X[cbeta.label[1]]
            df2 <- lmm$fixDF$X[cbeta.label[2]]
            result <- c(cbeta1, cbeta2, var1, var2, df1, df2, 
                arfactor, resid1, resid2)
        }
        cat("fmri.lme: calculating MFX model\n")
        if (cluster == 1) {
            if (total > 900) {
                calc.time <- round(total * 23/2000, 0)
                warning(paste("Time of calculation takes about", 
                  calc.time, "min, better use parallelizing option.", 
                  collapse = " "))
            }
            param <- apply(bold2, 1, funcB, z, fe.model, random, 
                weights, cbeta.label, dot = TRUE)
        }
        if (cluster > 1) {
            cl <- makeCluster(cluster)
            clusterEvalQ(cl, c(library(nlme), "funcB"))
            param <- parApply(cl, bold2, 1, funcB, z, fe.model, 
                random, weights, cbeta.label)
            stopCluster(cl)
        }
        cbeta[mask] <- param[1, ]
        cbeta2[mask] <- param[2, ]
        var[mask] <- param[3, ]
        var2[mask] <- param[4, ]
        df <- param[5, 1]
        df2 <- param[6, 1]
        arfactor[mask] <- param[7, ]
        resid[, mask] <- param[8:(dNt1 + 7), ]
        resid2[, mask] <- param[(8 + dNt1):(dNt2 + dNt1 + 7), 
            ]
        dim <- c(dim.spm, dNt1)
        dim2 <- c(dim.spm, dNt2)
        attr(cbeta, "coeff") <- cbeta.label[1]
        attr(var, "var") <- cbeta.label[1]
        attr(dim, "group label") <- gr.label[1]
        attr(cbeta2, "coeff") <- cbeta.label[2]
        attr(var2, "var") <- cbeta.label[2]
        attr(dim2, "group label") <- gr.label[2]
        random.model <- list(formula = paste(c(attr(random$subj, 
            "formula"), " | subj"), collapse = ""), classes = class(random$subj))
        cat("fmri.lme: finished\n")
        cat("fmri.lme: calculating spatial correlation\n")
        lags <- c(5, 5, 3)
        corr <- .Fortran("mcorr", as.double(resid), as.logical(mask), 
            as.integer(dx), as.integer(dy), as.integer(dz), as.integer(dNt1), 
            scorr = double(prod(lags)), as.integer(lags[1]), 
            as.integer(lags[2]), as.integer(lags[3]), PACKAGE = "fmri", 
            DUP = FALSE)$scorr
        dim(corr1) <- lags
        corr2 <- .Fortran("mcorr", as.double(resid2), as.logical(mask), 
            as.integer(dx), as.integer(dy), as.integer(dz), as.integer(dNt2), 
            scorr = double(prod(lags)), as.integer(lags[1]), 
            as.integer(lags[2]), as.integer(lags[3]), PACKAGE = "fmri", 
            DUP = FALSE)$scorr
        dim(corr2) <- lags
        qscale <- range(resid)
        scale <- max(abs(qscale))/32767
        resid <- writeBin(as.integer(resid/scale), raw(), 2)
        attr(resid, "resid") <- paste(c("res:group", gr.label[1]), 
            collapse = "")
        qscale2 <- range(resid2)
        scale2 <- max(abs(qscale2))/32767
        resid2 <- writeBin(as.integer(resid2/scale2), raw(), 
            2)
        attr(resid2, "resid") <- paste(c("res:group", gr.label[2]), 
            collapse = "")
        if (sum(!is.finite(corr))) {
            warning("fmri.lme:  infinite spatial correlations were found.\n            The analyzed region maybe too small. No bandwidths determinable.")
            bw <- NULL
            rxyz <- NULL
        }
        else {
            bw <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B", 
                lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 
                  6), lag = lags, data = corr)$par
            bw[bw <= 0.25] <- 0
            dim(bw) <- c(1, 3)
            rxyz <- c(resel(1, bw[1]), resel(1, bw[2]), resel(1, 
                bw[3]))
            dim(rxyz) <- c(1, 3)
        }
        if (sum(!is.finite(corr2))) {
            bw2 <- NULL
            rxyz2 <- NULL
        }
        else {
            bw2 <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B", 
                lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 
                  6), lag = lags, data = corr2)$par
            bw2[bw2 <= 0.25] <- 0
            dim(bw2) <- c(1, 3)
            rxyz2 <- c(resel(1, bw2[1]), resel(1, bw2[2]), resel(1, 
                bw2[3]))
            dim(rxyz2) <- c(1, 3)
        }
        attr(scale, "resid") <- paste(c("res:group", gr.label[1]), 
            collapse = "")
        attr(corr, "spatial correlations:group") <- gr.label[1]
        attr(bw, "bandwidths:group") <- gr.label[1]
        attr(rxyz, "smoothness in resel space:group") <- gr.label[1]
        attr(scale2, "resid") <- paste(c("res:group", gr.label[2]), 
            collapse = "")
        attr(corr2, "spatial correlations:group") <- gr.label[2]
        attr(bw2, "bandwidths:group") <- gr.label[2]
        attr(rxyz2, "smoothness in resel space:group") <- gr.label[2]
        cat("fmri.lme: finished\n")
    }
    else if (gr > 2) {
        stop("fmri.lme: Sorry, more than two groups are not planned in this version.")
    }
    cat("fmri.lme: exiting function\n")
    result <- list(cbeta = cbeta, var = var, cbeta2 = cbeta2, 
        var2 = var2, mask = mask, res = resid, resscale = scale, 
        res2 = resid2, resscale2 = scale2, arfactor = arfactor, 
        rxyz = rxyz, scorr = corr, rxyz2 = rxyz2, scorr2 = corr2, 
        weights = wghts, dim = dim, dim2 = dim2, bw = bw, bw2 = bw2, 
        df = df, df2 = df2)
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
    result$subjects <- max(as.numeric(z$subj))
    result$subj.runs <- max(as.numeric(z$run))
    result$sessions <- max(as.numeric(z$session))
    result$groups <- gr
    result$fixedModel <- paste(c("", as.character(fixed)), collapse = " ")
    result$randomModel <- random.model
    result$VarModel <- vtype
    result$cluster <- cluster
    class(result) <- c("fmridata", "fmrispm")
    attr(result, "design") <- z
    attr(result, "estimator") <- "one-stage process"
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
