.First.lib <- function(lib, pkg) {
  if(version$major==0)
    stop("This version for R 1.00 or later")
  library.dynam("fmri", pkg, lib)
}
