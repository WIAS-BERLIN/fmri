## Group-Designmatrix, fmri-package is needed

fmri.designG <- function (hrf, subj = 1, runs = 1)
{
  if (is.null(dim(hrf))) dim(hrf) <- c(length(hrf), 1)
  stimuli <- dim(hrf)[2]
  scans <- dim(hrf)[1]
  id <- c(rep(1, scans*runs))
  if (subj > 1) {
    for (s in 2:subj) {id<-c(id,c(rep(s*1,scans*runs)))}
  }
  id <- factor(id)
  run <- c(rep(1,scans))
  if (runs > 1) {
    for (r in 2:runs) {run<-c(run,c(rep(r*1,scans)))}
  }
  run <- factor(run)
  run.g <- rep(run,subj)
  scan <- c(1:scans)
  scan.g <- rep(scan,subj*runs)
  session <- c(rep(1, scans))
  s <- 1
  while (s < subj*runs) {
    s <- s+1
    session <- c(session,c(rep(s,scans)))
  }
  session <- factor(session)
  hrf.list="hrf"
  if (stimuli > 1) {
    for (s in 2:stimuli) {
      list2<-c(paste("hrf",s,sep=""))
      hrf.list<-c(hrf.list,list2)
    }
  }
  x <- fmri.design(hrf, order=2) # fmri-package used
  x.g <- rep(1,subj*runs) %x% x
  x.group <- data.frame(id, run.g, scan.g, session, x.g)
  colnames(x.group)<-c("subj","run","scan","session",hrf.list,"drift0","drift1","drift2")
  x.group
}

## Output-Matrix
# subj: consecutive subject number, 1 to subj
# run:  consecutive run number within subj, 1 to run, 
#       for experiments with repeated measurents
# scan: consecutive scan number, 1 to T, for every single experiment
# session: consecutive experiment number, 1 to (subj*run)
# hrf: expected BOLD-response
# drift0, drift1, drift2: polynomial drift terms, orthogonal to the stimuli



## Estimation of the group map, nlme- and parallel-packages are needed

# Input data:
# bold  = Multi-Subject-Boldsignal-Array (Bold time series)
# z     = Group-Designmatrix, typically output from fmri.designG(.)
# fixed = fixed effects, one-sided formula, 
#       e.g. fixed <- ~ 0 + subj+drift1:subj+drift2:subj
# random = random effects, one-sided formula with grouping factor(s)
#       e.g. random <- ~ 0 + hrf|subj
# ac    = AR1-coefficent(s), single value (global) or voxel specific array (local)
# mask  = head mask, if available
# vtype = homoscedastic ("eqaul") or heteroscedastic/subject specific
#        ("individual") residuals
# cluster = number of cores for parallel estimations

fmri.lmePar <- function (bold, z, fixed, random, ac = 0.3, 
                      mask = NULL, vtype ="individual", cluster = 2)
{
  t1 <- Sys.time()
  if (length(dim(bold)) != 4) {
    stop("Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
  }
  if (dim(bold)[4] != dim(z)[1]) {
    stop("the design matrix do not match to the bold array")
  }
  dim.spm <- dim(bold)[1:3]
  dx <- dim(bold)[1]
  dy <- dim(bold)[2]
  dz <- dim(bold)[3]
  dNt <- dim(bold)[4]
  total0 <- dx*dy*dz
  if ((!length(ac==1))||(!dim(ac)==dim.spm)) {
    stop("Error: dimension of the array with the AR1-coefficients 
         do not match to the data set")
  }
  
  # Anwenden der Kopfmaske auf die Daten
  dim(bold) <- c(total0,dNt)
  if (is.null(mask)) {
    mask <- array(TRUE, dim=dim.spm)
  }
  total <- sum(mask)
  bold2 <- array(0, dim=c(total, dNt+1))
  bold2[, 1] <- if(length(ac)==total0) ac[mask] else rep(ac[1],total)
  bold2[, -1] <- bold[mask,]
  rm(ac, bold)
  gc()
  
  # leere Ausgabe-Arrays anlegen    
  cbeta <- array(0, dim = dim.spm)
  var <- array(1, dim = dim.spm)
  resid <- array(1, dim = c(dNt, total0))
  arfactor <- array(0, dim = dim.spm)
  df <- 100
  
  # Modellsyntax konfigurieren
  fe.model <- formula(paste(c("y",as.character(fixed)),collapse = " "))
  if (vtype=="individual") {
    weights <- varIdent(form =~ 1|subj)
  }
  if (vtype=="equal") {
    weights <- NULL
  }
  
  # innere Hilfsfunktion: Schaetzung der Modellparameter im Voxel i
  funB <- function(y, z, fe.model, random, weights) {
    z$y <- y[-1]
    lmm <-  lme(fixed = fe.model,
                random = random,           
                correlation = corAR1(value=y[1], form = ~1|subj, fixed=TRUE),                   
                weights = weights, 
                method ="REML",
                control = lmeControl(rel.tol=1e-6, returnObject = TRUE),
                data = z)
    cbeta <- fixef(lmm)[1]
    var <- vcov(lmm)[1,1]
    resid <- resid(lmm)
    arfactor <- y[1]
    df <- lmm$fixDF$X[1]
    result <- c(cbeta,var,arfactor,df,resid)
  }
  
  ## voxelweise Auswertung
  cl <- makeCluster(cluster)
  clusterEvalQ(cl, c(library(nlme),"funB"))
  param <- parApply(cl, bold2, 1, funB, z, fe.model, random, weights)
  stopCluster(cl)
  #param <- parRapply(cl, bold2, funB, z, fe.model, random, weights)
  # gibt keinen Array, sondern eine Liste zurueck, macht deswegen:
  # Error in param[1, ] : incorrect number of dimensions
  cbeta[mask] <- param[1, ]
  var[mask] <- param[2, ]
  arfactor[mask] <- param[3, ]
  df <- param[4,1]
  resid[, mask] <- param[5:(dNt+4), ]
  qscale <- range(resid)
  scale <- max(abs(qscale))/32767
  resid <- writeBin(as.integer(resid/scale), raw(), 2)
  result <- list(cbeta = cbeta, var = var,
                 res = resid, resscale=scale,
                 arfactor = arfactor, df = df, 
                 mask = mask, dim = c(dim.spm,dNt))
  class(result) <- c("fmridata", "fmrispm")
  
  t2 <- Sys.time()
  cat("time elapsed",difftime(t2,t1,units="mins"),"mins\n")
  invisible(result)
  }


# same, but no parallelizing
fmri.lme <- function (bold, z, fixed, random, ac = 0.3, 
                      mask = NULL, vtype ="individual")
{
  t1 <- Sys.time()
  if (length(dim(bold)) != 4) {
    stop("Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
  }
  if (dim(bold)[4] != dim(z)[1]) {
    stop("the design matrix do not match to the bold array")
  }
  dim.spm <- dim(bold)[1:3]
  dx <- dim(bold)[1]
  dy <- dim(bold)[2]
  dz <- dim(bold)[3]
  dNt <- dim(bold)[4]
  total0 <- dx*dy*dz
  if ((!length(ac==1))||(!dim(ac)==dim.spm)) {
    stop("Error: dimension of the array with the AR1-coefficients 
         do not match to the data set")
  }
  
  # Anwenden der Kopfmaske auf die Daten
  dim(bold) <- c(total0,dNt)
  if (is.null(mask)) {
    mask <- array(TRUE, dim=dim.spm)
  }
  total <- sum(mask)
  bold2 <- array(0, dim=c(total, dNt+1))
  bold2[, 1] <- if(length(ac)==total0) ac[mask] else rep(ac[1],total)
  bold2[, -1] <- bold[mask,]
  rm(ac, bold)
  gc()
  
  # leere Ausgabe-Arrays anlegen    
  cbeta <- array(0, dim = dim.spm)
  var <- array(1, dim = dim.spm)
  resid <- array(1, dim = c(dNt, total0))
  arfactor <- array(0, dim = dim.spm)
  df <- 100
  
  # Modellsyntax konfigurieren
  fe.model <- formula(paste(c("y",as.character(fixed)),collapse = " "))
  if (vtype=="individual") {
    weights <- varIdent(form =~ 1|subj)
  }
  if (vtype=="equal") {
    weights <- NULL
  }
  
  # innere Hilfsfunktion: Schaetzung der Modellparameter im Voxel i
  funB <- function(y, z, fe.model, random, weights) {
    z$y <- y[-1]
    #i <<- i+1
    cat(".")
    #setTxtProgressBar(pb, i)
    #cat("fmri.lme : calculating voxel: ", i, "/", total,"     \r", sep="")
    lmm <-  lme(fixed = fe.model,
                random = random,           
                correlation = corAR1(value=y[1], form = ~1|subj, fixed=TRUE),                   
                weights = weights, 
                method ="REML",
                control = lmeControl(rel.tol=1e-6, returnObject = TRUE),
                data = z)
    cbeta <- fixef(lmm)[1]
    var <- vcov(lmm)[1,1]
    resid <- resid(lmm)
    arfactor <- y[1]
    df <- lmm$fixDF$X[1]
    result <- c(cbeta,var,arfactor,df,resid)
  }
  
  # voxelweise Auswertung
  #pb <- txtProgressBar(0, total-1, style = 3)
  #i <- 0 # Laufzeitindikator
  param <- apply(bold2, 1, funB, z, fe.model, random, weights)  
  cbeta[mask] <- param[1, ]
  var[mask] <- param[2, ]
  arfactor[mask] <- param[3, ]
  df <- param[4,1]
  resid[, mask] <- param[5:(dNt+4), ]
  #close(pb)
  qscale <- range(resid)
  scale <- max(abs(qscale))/32767
  resid <- writeBin(as.integer(resid/scale), raw(), 2)
  result <- list(cbeta = cbeta, var = var,
                 res = resid, resscale=scale,
                 arfactor = arfactor, df = df, 
                 mask = mask, dim = c(dim.spm,dNt))
  class(result) <- c("fmridata", "fmrispm")
  t2 <- Sys.time()
  cat("time elapsed",difftime(t2,t1,units="mins"),"mins\n")
  invisible(result)
  }
