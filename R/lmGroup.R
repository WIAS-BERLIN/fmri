## Construct one-stage Group-Design-Matrix, fmri-package is needed

fmri.designG <- function (hrf, 
                          subj = 1, 
                          runs = 1, 
                          group = NULL,
                          age = NULL,
                          sex = NULL,
                          iq = NULL)
{
  if (is.null(dim(hrf))) 
    dim(hrf) <- c(length(hrf), 1)
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
  
  if (!is.null(group)) {
    if (length(group) != subj) {
      stop("fmri.designG: Length of group vector: ",length(group), 
           " does not match the number of subjects: " ,subj, ".")
    } else {
      group.g <- NULL
      for (s in 1:subj) {group.g <- c(group.g,c(rep(group[s],scans*runs)))}
    }
  } else {
    group.g <- c(rep(1,scans*runs*subj))
  }
  group.g <- factor(group.g)
  
  hrf.list="hrf"
  if (stimuli > 1) {
    for (s in 2:stimuli) {
      list2<-c(paste("hrf",s,sep=""))
      hrf.list<-c(hrf.list,list2)
    }
  }
  x <- fmri.design(hrf, order=2) # fmri-package used
  x.g <- rep(1,subj*runs) %x% x
  
  ## add population traits
  if (!is.null(age)) {
    if (length(age) != subj) {
      stop("fmri.designG: Length of group vector: ",length(age), 
           " does not match the number of subjects: " ,subj, ".")
    } else {
      age.g <- NULL
      for (s in 1:subj) {age.g <- c(age.g,c(rep(age[s],scans*runs)))}
    }
  } else {
    age.g <- c(rep(NA,scans*runs*subj))
  }
  age.gc <- age.g - mean(age.g) # centered
  
  if (!is.null(sex)) {
    if (length(sex) != subj) {
      stop("fmri.designG: Length of group vector: ",length(sex), 
           " does not match the number of subjects: " ,subj, ".")
    } else {
      sex.g <- NULL
      for (s in 1:subj) {sex.g <- c(sex.g,c(rep(sex[s],scans*runs)))}
    }
  } else {
    sex.g <- c(rep(NA,scans*runs*subj))
  }
  sex.gc <- sex.g - mean(sex.g) # centered
  sex.g <- factor(sex.g)
  
  if (!is.null(iq)) {
    if (length(iq) != subj) {
      stop("fmri.designG: Length of group vector: ",length(iq), 
           " does not match the number of subjects: " ,subj, ".")
    } else {
      iq.g <- NULL
      for (s in 1:subj) {iq.g <- c(iq.g,c(rep(iq[s],scans*runs)))}
    }
  } else {
    iq.g <- c(rep(NA,scans*runs*subj))
  }
  iq.gc <- iq.g - mean(iq.g) # centered
  
  x.group <- data.frame(id, run.g, scan.g, session, group.g, x.g, 
                        age.g, age.gc, sex.g, sex.gc, iq.g, iq.gc)
  colnames(x.group)<-c("subj","run","scan","session","group",
                       hrf.list,"drift0","drift1","drift2",
                       "age","age_m0","sex","sex_m0", "iq","iq_m0")  
  if (runs > 1) 
    warning("If not all runs within each subject were completed, please delete this rows manually!")
  x.group
}

# Input data:
# hrf = one or more expected BOLD response(s) for one session, 
#       typically output from fmri.stimulus(.)
# subj = number of subjects in the study
# runs = number of repeated measures within subjects
# group = optional group vector, length == subj
# age = optional vector with ages for all subjects, length(age) == subj
# sex = optional vector with gender for all subjects, length(sex) == subj
# iq  = optional vector with IQ or other test scores for all subjects, length(iq) == subj

## Output-Matrix (data.frame)
# subj: consecutive subject number, 1 to subj
# run:  consecutive run number within subj, 1 to run, 
#       for experiments with repeated measures
# scan: consecutive scan number, 1 to T, for every single experiment
# session: consecutive experiment number, 1 to (subj*run)
# group: grouping variable, default: 1 group only
# hrf: expected BOLD-response
# drift0, drift1, drift2: polynomial drift terms, orthogonal to the stimuli
# age: multiplied age information, if unknown: NAs
# age_m0: like above, but mean centered by subtracting the overall mean
# sex: multiplied sex information, if unknown: NAs
# sex_m0: like above, but mean centered by subtracting the overall mean
# iq:  multiplied IQ information, if unknown: NAs
# iq_m0: like above, but mean centered by subtracting the overall mean


## Estimation of the group map, nlme- and parallel-packages are needed

# Input data:
# bold  = Multi-Subject-Boldsignal-Array (Bold time series),
#         attention: order of cumulative data sets is importent! 
#         (all runs/sessions of first subject, all runs of second subject, ...)
# z     = Group-Design-Matrix, typically output from fmri.designG(.)
# fixed = fixed effects, one-sided formula,
#         no intercept because intercept will be a linear combination of 
#         the session factor-variable, modeling session-specific intercepts
#       e.g. fixed <- ~ 0 + hrf + session + drift1:session + drift2:session 
#            (default in case of one group of subjects)
#       e.g. fixed <- ~ 0 + hrf:group + session + drift1:session + drift2:session
#            (default in case of two groups of subjects, estimates one SPM per group)
#
# random = random effects, one-sided formula with grouping factor(s)
#       e.g. random <- ~ 0 + hrf|subj 
#            (default in case of one group without repeated measures)
#       e.g. random <- ~ 0 + hrf|subj/session 
#            (default in case of one group with repeated measures)
#       e.g. random <- list(subj = pdDiag(~ 0 + hrf:group)) 
#            (default in case of two groups without repeated measures) 
# ac    = AR1-coefficent(s), single value (global) or voxel specific array (local)
# mask  = head mask, if available
# vtype = homoscedastic ("eqaul") or heteroscedastic/subject-specific
#        ("individual") residuals
# cluster = number of CPU cores for parallel voxel calculating,
#         cluster = 1 used the simple unparallelized method (very slow);
#         if you do not know your CPUs, try: detectCores() from parallel-package
# wghts = isotropic voxel (MNI-space: wghts = c(1,1,1), default) or other vector

fmri.lmePar <- function (bold,
                         z,
                         fixed = NULL,
                         random = NULL,
                         mask = NULL,
                         ac = 0.3, 
                         vtype ="individual",
                         cluster = 2,
                         wghts = c(1,1,1)) {
  
  t1 <- Sys.time()
  cat("fmri.lme: entering function\n")
  
  ## set dimensions
  dim.spm <- dim(bold)[1:3]
  dx <- dim(bold)[1]
  dy <- dim(bold)[2]
  dz <- dim(bold)[3]
  dNt <- dim(bold)[4]
  total0 <- dx*dy*dz
  
  ## test dimensionality of objects and design matrix
  if (length(dim(bold)) != 4)
    stop("fmri.lme: Hmmmm, this does not seem to be a fMRI time series. I better stop executing! Sorry!\n")
  if (!is.data.frame(z))
    stop("fmri.lme: Something is wrong with the design matrix.")
  if (dim(bold)[4] != dim(z)[1])
    stop("fmri.lme: design matrix dimensionality does not match functional data.")
  if ((!length(ac)==1) && (!identical(dim(ac),dim.spm)))
    stop("fmri.lme: AR1-array dimensionality does not match functional data.")
  if (is.null(mask)) {
    mask <- array(TRUE, dim=dim.spm)
  }
  if (!identical(dim(mask),dim.spm))
    stop("fmri.lme: mask dimensionality does not match functional data.")
  if (sum(mask)==0)
    stop("fmri.meta: your head mask is empty, nothing to do.")
  
  ## keep roi-attributes, bold object will be removed
  if (!is.null(attr(bold, "xind"))) {
    xind <- attr(bold, "xind")
  } else { 
    xind <- NULL
  }
  if (!is.null(attr(bold, "yind"))) {
    yind <- attr(bold, "yind")
  } else {
    yind <- NULL
  }
  if (!is.null(attr(bold, "zind"))) {
    zind <- attr(bold, "zind")
  } else {
    zind <- NULL
  }
  
  ## keep ac-type, ac object will be removed
  if (length(ac)==1) {
    white <- 3 # actype == "ac"     ... AC only 
  } else {
    white <- 1 # actype == "smooth" ... calc AC, smooth AC, calc prewhitened model 
  } 
  
  ## apply mask and merge ac-array with functional data 
  dim(bold) <- c(total0,dNt)
  total <- sum(mask)
  
  ## insert: early termination criterion
  if ((cluster==1) && (total > 26000)) {
    stop("Time for this calculation takes longer than 5 hours. Please use parallelizing option.")
  } else if ((cluster==1) && (total > 16000)) {
    stop("Time for this calculation takes longer than 3 hours. Please use parallelizing option.") 
  }
  
  bold2 <- array(0, dim=c(total, dNt+1))
  bold2[, 1] <- if(length(ac)==total0) ac[mask] else rep(ac[1],total)
  bold2[, -1] <- bold[mask,]
  rm(ac, bold)
  gc()
  
  ## define some arrays capturing the estimated values
  dim <- c(dim.spm,dNt)
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
  
  gr <- nlevels(z$group) # get number of groups
  runs <- nlevels(z$run) # get number of repeated measures
  
  ## one group only
  if (gr==1) {
    
    ## configure model syntax
    if (is.null(fixed)) {
      fixed <- ~ 0 + hrf + session + drift1:session + drift2:session
    }
    fe.model <- formula(paste(c("y",as.character(fixed)),collapse = " "))
    # session-specific intercepts strictly necessary, mostly high variable;
    # subject-specific intercepts not sufficient, even in case of repeated
    # measurements model missfit were found under this restiction (dames: LR-tests)
        
    if (is.null(random)) {
      if (runs == 1) {
        random <- ~ 0 + hrf|subj
        ar1 <- ~ 1|subj/session
      }
      if (runs > 1) {
        if (subj==1) {
          random <- ~ 0 + hrf|session
          ar1 <- ~ 1|session
        } else {
        random <- ~ 0 + hrf|subj/session
        ar1 <- ~ 1|subj/session
        }
      }
    } else if (subj==1) {
      ar1 <- ~ 1|session
    } else {
      ar1 <- ~ 1|subj/session
    } 
          
    # ... but subject-specific error terms still be valid, whereas
    # homoscedasticity assumption mostly not acceptable (dames: LR-tests) 
    if (vtype=="individual") {
      weights <- varIdent(form =~ 1|subj)
    }
    if (vtype=="equal") {
      weights <- NULL
    } 
    
    ## internal auxiliary function: estimation of the parameters in one voxel i
    funcA <- function(y, z, fe.model, random, ar1, weights, dot = FALSE) {
      z$y <- y[-1]
      if (dot) cat(".")
      lmm <-  lme(fixed = fe.model,
                  random = random,           
                  correlation = corAR1(value=y[1], form = ar1, fixed=TRUE),                   
                  weights = weights, 
                  method ="REML",
                  control = lmeControl(rel.tol=1e-6, returnObject = TRUE),
                  data = z)
      cbeta <- fixef(lmm)[1] # alternatively: fixef(lmm)["hrf"]
      var <- vcov(lmm)[1,1]  # alternatively: vcov(lmm)["hrf","hrf"]
      resid <- resid(lmm)
      arfactor <- y[1]
      df <- lmm$fixDF$X[1]   # alternatively: lmm$fixDF$X["hrf"]
      result <- c(cbeta,var,arfactor,df,resid)
    }
    
    cat("fmri.lme: calculating MFX model\n")
    
    ## voxel-by-voxel analysis
    # A) without parallelizing 
    if (cluster==1) {
      if (total > 900) {
        calc.time <- round(total*23/2000,0)  # ca. 23 min per 2000 voxels (dames)
        warning(paste("Time of calculation takes about",calc.time,
                      "min, better use parallelizing option.", collapse = " "))
      }
      param <- apply(bold2, 1, funcA, z, fe.model, random, ar1, weights, dot = TRUE)
    } 
    
    # B) with parallelizing
    if (cluster > 1) {
      cl <- makeCluster(cluster)
      clusterEvalQ(cl, c(library(nlme),"funcA"))
      param <- parApply(cl, bold2, 1, funcA, z, fe.model, random, ar1, weights)
      stopCluster(cl)
      # param <- parRapply(cl, bold2, funcA, z, fe.model, random, ar1, weights)
      # gives back no array but a list, therefore getting:
      # Error in param[1, ] : incorrect number of dimensions
    }
    
    cbeta[mask] <- param[1, ]
    var[mask] <- param[2, ]
    arfactor[mask] <- param[3, ]
    df <- param[4,1]
    resid[, mask] <- param[5:(dNt+4), ]
    random.model <- paste(c("",as.character(random)),collapse = " ")
    cat("fmri.lme: finished\n")
    
    cat("fmri.lme: calculating spatial correlation\n")
    lags <- c(5, 5, 3)
    corr <- .Fortran("mcorr", as.double(resid), as.logical(mask), 
                     as.integer(dx), as.integer(dy), as.integer(dz), 
                     as.integer(dNt), scorr = double(prod(lags)), as.integer(lags[1]), 
                     as.integer(lags[2]), as.integer(lags[3]), PACKAGE = "fmri", 
                     DUP = FALSE)$scorr
    dim(corr) <- lags
    # attention: NaNs in corr
    # if analyzed brain region is too small (e.g. 1 or 2 Slices)
    
    ## "compress" the residuals
    qscale <- range(resid)
    scale <- max(abs(qscale))/32767
    resid <- writeBin(as.integer(resid/scale), raw(), 2)
    
    ## determine local smoothness  
    if (sum(!is.finite(corr))) {
      warning("fmri.lme: infinite spatial correlations were found. 
              The analyzed region maybe too small. No bandwidths determinable.")
      bw <- NULL
      rxyz <- NULL
    } else {
      bw <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B", 
                  lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 6), 
                  lag = lags, data = corr)$par
      bw[bw <= 0.25] <- 0
      dim(bw) <- c(1, 3)
      rxyz <- c(resel(1, bw[1]), resel(1, bw[2]), resel(1, bw[3]))
      dim(rxyz) <- c(1, 3)
    }
    scale2 <- NULL
    corr2 <- NULL
    bw2 <- NULL
    rxyz2 <- NULL
    dim2 <- NULL
    cat("fmri.lme: finished\n")
  }
  
  ## two independent groups
  if (gr==2) {
    cbeta2 <- cbeta
    var2 <- var
    gr.label <- levels(z$group)
    dNt1 <- length(z$group[z$group==gr.label[1]])
    dNt2 <- length(z$group[z$group==gr.label[2]])
    resid <- array(1, dim = c(dNt1, total0)) 
    resid2 <- array(1, dim = c(dNt2, total0)) 
    df2 <- df
    
    ## configure model syntax
    if (is.null(fixed)) {
      fixed <- ~ 0 + hrf:group + session + drift1:session + drift2:session
    } 
    fe.model <- formula(paste(c("y",as.character(fixed)),collapse = " "))

    if (is.null(random)) {
      if (runs == 1) {
        random <- list(subj = pdDiag(~ 0 + hrf:group)) # independent groups
      } else if (runs > 1) {
        stop("fmri.lme: Sorry, repeated measures for two groups are not predefined.")
      }
    }
    if (vtype=="individual") {
      weights <- varIdent(form =~ 1|subj)
    }
    if (vtype=="equal") {
      weights <- NULL
    }
    
    # nlme:lme intern labeling 
    cbeta.label <- c(paste(c("hrf:group",gr.label[1]), collapse = ""),
                     paste(c("hrf:group",gr.label[2]), collapse = ""))  
    
    ## internal auxiliary function: estimation of the parameters in one voxel i
    funcB <- function(y, z, fe.model, random, weights, cbeta.label, dot = FALSE) 
    {
      z$y <- y[-1]
      if (dot) cat(".")
      lmm <-  lme(fixed = fe.model,
                  random = random,           
                  correlation = corAR1(value=y[1], form = ~ 1|subj/session, fixed=TRUE),                   
                  weights = weights, 
                  method ="REML",
                  control = lmeControl(rel.tol=1e-6, returnObject = TRUE),
                  data = z)
      
      cbeta1 <- fixef(lmm)[cbeta.label[1]]
      cbeta2 <- fixef(lmm)[cbeta.label[2]]
      var1 <- vcov(lmm)[cbeta.label[1],cbeta.label[1]]
      var2 <- vcov(lmm)[cbeta.label[2],cbeta.label[2]]
      resid1 <- resid(lmm)[z$group==levels(z$group)[1]]
      resid2 <- resid(lmm)[z$group==levels(z$group)[2]]
      arfactor <- y[1]
      df1 <- lmm$fixDF$X[cbeta.label[1]]
      df2 <- lmm$fixDF$X[cbeta.label[2]]
      result <- c(cbeta1,cbeta2,var1,var2,df1,df2,arfactor,resid1,resid2)
    }
    
    cat("fmri.lme: calculating MFX model\n")
    
    ## voxel-by-voxel analysis
    # A) without parallelizing 
    if (cluster==1) {
      if (total > 900) {
        calc.time <- round(total*23/2000,0)  # ca. 23 min per 2000 voxels (dames)
        warning(paste("Time of calculation takes about",calc.time,
                      "min, better use parallelizing option.", collapse = " "))
      }
      param <- apply(bold2, 1, funcB, z, fe.model, random, weights, 
                     cbeta.label, dot = TRUE)
    } 
    
    # B) with parallelizing
    if (cluster > 1) {
      cl <- makeCluster(cluster)
      clusterEvalQ(cl, c(library(nlme),"funcB"))
      param <- parApply(cl, bold2, 1, funcB, z, fe.model, random, weights,
                        cbeta.label)
      stopCluster(cl)
    }
    
    cbeta[mask] <- param[1, ]
    cbeta2[mask] <- param[2, ]
    var[mask] <- param[3, ]
    var2[mask] <- param[4, ]
    df <- param[5,1]
    df2 <- param[6,1]
    arfactor[mask] <- param[7, ]
    resid[, mask] <- param[8:(dNt1+7), ]
    resid2[, mask] <- param[(8+dNt1):(dNt2+dNt1+7), ]
    dim <- c(dim.spm,dNt1)
    dim2 <- c(dim.spm,dNt2)
    
    attr(cbeta, "coeff") <- cbeta.label[1]
    attr(var, "var") <- cbeta.label[1]
    attr(dim, "group label") <- gr.label[1]
    attr(cbeta2, "coeff") <- cbeta.label[2]
    attr(var2, "var") <- cbeta.label[2]   
    attr(dim2, "group label") <- gr.label[2]
    random.model <- list(formula=paste(c(attr(random$subj, "formula"), " | subj"),
                                 collapse = ""), 
                   classes=class(random$subj))
    cat("fmri.lme: finished\n")
    
    cat("fmri.lme: calculating spatial correlation\n")
    lags <- c(5, 5, 3)
    corr <- .Fortran("mcorr", as.double(resid), as.logical(mask), 
                      as.integer(dx), as.integer(dy), as.integer(dz), 
                      as.integer(dNt1), scorr = double(prod(lags)), as.integer(lags[1]), 
                      as.integer(lags[2]), as.integer(lags[3]), PACKAGE = "fmri", 
                      DUP = FALSE)$scorr
    dim(corr1) <- lags
    
    corr2 <- .Fortran("mcorr", as.double(resid2), as.logical(mask), 
                      as.integer(dx), as.integer(dy), as.integer(dz), 
                      as.integer(dNt2), scorr = double(prod(lags)), as.integer(lags[1]), 
                      as.integer(lags[2]), as.integer(lags[3]), PACKAGE = "fmri", 
                      DUP = FALSE)$scorr
    dim(corr2) <- lags
    # attention: NaNs in corr
    # if analyzed brain region is too small (e.g. 1 or 2 Slices)
    
    ## "compress" the residuals
    qscale <- range(resid)
    scale <- max(abs(qscale))/32767
    resid <- writeBin(as.integer(resid/scale), raw(), 2)
    attr(resid, "resid") <- paste(c("res:group",gr.label[1]), collapse = "")
    
    qscale2 <- range(resid2)
    scale2 <- max(abs(qscale2))/32767
    resid2 <- writeBin(as.integer(resid2/scale2), raw(), 2)
    attr(resid2, "resid") <- paste(c("res:group",gr.label[2]), collapse = "")  
    
    ## determine local smoothness  
    if (sum(!is.finite(corr))) {
      warning("fmri.lme:  infinite spatial correlations were found.
            The analyzed region maybe too small. No bandwidths determinable.")
      bw <- NULL
      rxyz <- NULL
    } else {
      bw <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B", 
                  lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 6), 
                  lag = lags, data = corr)$par
      bw[bw <= 0.25] <- 0
      dim(bw) <- c(1, 3)
      rxyz <- c(resel(1, bw[1]), resel(1, bw[2]), resel(1, bw[3]))
      dim(rxyz) <- c(1, 3)
    }
    
    if (sum(!is.finite(corr2))) {
      bw2 <- NULL
      rxyz2 <- NULL
    } else {
      bw2 <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B", 
                   lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 6), 
                   lag = lags, data = corr2)$par
      bw2[bw2 <= 0.25] <- 0
      dim(bw2) <- c(1, 3)
      rxyz2 <- c(resel(1, bw2[1]), resel(1, bw2[2]), resel(1, bw2[3]))
      dim(rxyz2) <- c(1, 3)
    }
    attr(scale, "resid") <- paste(c("res:group",gr.label[1]), collapse = "")
    attr(corr, "spatial correlations:group") <- gr.label[1]
    attr(bw, "bandwidths:group") <- gr.label[1]
    attr(rxyz, "smoothness in resel space:group") <- gr.label[1]
    attr(scale2, "resid") <- paste(c("res:group",gr.label[2]), collapse = "")
    attr(corr2, "spatial correlations:group") <- gr.label[2]
    attr(bw2, "bandwidths:group") <- gr.label[2]
    attr(rxyz2, "smoothness in resel space:group") <- gr.label[2]    
    cat("fmri.lme: finished\n")

  } else if (gr > 2) {
    stop("fmri.lme: Sorry, more than two groups are not planned in this version.") 
  }

  cat("fmri.lme: exiting function\n")
  result <- list(cbeta = cbeta, 
                 var = var,
                 cbeta2 = cbeta2,
                 var2 = var2,
                 mask = mask,
                 res = resid, 
                 resscale = scale,
                 res2 = resid2,
                 resscale2 = scale2,
                 arfactor = arfactor, 
                 rxyz = rxyz, 
                 scorr = corr,
                 rxyz2 = rxyz2, 
                 scorr2 = corr2,
                 weights = wghts, 
                 dim = dim,
                 dim2 = dim2, 
                 bw = bw,
                 bw2 = bw2, 
                 df = df,
                 df2 = df2)
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
  ## save further group study objects
  result$subjects <- max(as.numeric(z$subj))
  result$subj.runs <- max(as.numeric(z$run))
  result$sessions <- max(as.numeric(z$session))
  result$groups <- gr
  result$fixedModel <- paste(c("",as.character(fixed)),collapse = " ")
  result$randomModel <- random.model
  result$VarModel <- vtype
  result$cluster <- cluster
  #result$vwghts <- 1 # unused argument since fmri v1.6-x 
  
  class(result) <- c("fmridata", "fmrispm")
  attr(result, "design") <- z
  attr(result, "estimator") <- "one-stage process"
  attr(result, "white") <- white
  attr(result, "residuals") <- !is.null(scale)
  
  t2 <- Sys.time()
  cat("time elapsed",difftime(t2,t1,units="mins"),"mins\n")
  invisible(result)
}
