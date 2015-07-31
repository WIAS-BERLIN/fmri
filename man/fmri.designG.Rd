\name{fmri.designG}
\alias{fmri.designG}
\title{
  Construct one-stage Group-Design-Matrix
}
\description{
  Construct one-stage Group-Design-Matrix 
}
\usage{
fmri.designG(hrf, subj = 1, runs = 1, group = NULL, age = NULL, sex = NULL, iq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hrf}{one or more expected BOLD response(s) for one session, 
             typically output from fmri.stimulus(.)}
  \item{subj}{number of subjects in the study}
  \item{runs}{number of repeated measures within subjects}
  \item{group}{optional group vector, length == subj}
  \item{age}{optional vector with ages for all subjects, length(age) == subj}
  \item{sex}{optional vector with gender for all subjects, length(sex) == subj}
  \item{iq}{optional vector with IQ or other test scores for all subjects, length(iq) == subj}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  The value of this function is a dataframe with the following variables
    \item{subj}{consecutive subject number, 1 to subj}
    \item{run}{consecutive run number within subj, 1 to run, 
       for experiments with repeated measures}
    \item{scan}{consecutive scan number, 1 to T, for every single experiment}
    \item{session}{consecutive experiment number, 1 to (subj*run)}
    \item{group}{grouping variable, default: 1 group only}
    \item{hrf}{expected BOLD-response}
    \item{drift0, drift1, drift2}{polynomial drift terms, orthogonal to the stimuli}
    \item{age}{multiplied age information, if unknown: NAs}
    \item{age_m0}{like above, but mean centered by subtracting the overall mean}
    \item{sex}{multiplied sex information, if unknown: NAs}
    \item{sex_m0}{like above, but mean centered by subtracting the overall mean}
    \item{iq}{multiplied IQ information, if unknown: NAs}
    \item{iq_m0}{like above, but mean centered by subtracting the overall mean}
}
\author{
  Sibylle Dames
}
\seealso{
  \code{\link{fmri.stimulus}}, \code{\link{fmri.design}}, \code{\link{fmri.lmePar}}
}
\examples{
subj <- 6
runs <- 1
scans <- 121
times <- c(12, 48, 84, 120, 156, 192, 228, 264) 
duration <- 24          
tr <- 2.5

## Group-Designmatrix
hrf <- fmri.stimulus(scans, times, duration, tr, times = TRUE)
x.group <- fmri.designG(hrf, subj = subj, runs = runs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ design }
