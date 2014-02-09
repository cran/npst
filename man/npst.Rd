\name{npst}
\alias{npst}
\title{
  npst
}
\description{
  Package 'npst' generalizes Hewitt's (1971) test for seasonality and
  Rogerson's (1996) extension based on Monte-Carlo simulation.
}
\usage{
npst(indata=NULL, long=12, peak=6, repts=100000,
  whole.distribution=FALSE, siglevels=c(0.001,0.01,0.05,0.1),
  PARALLEL=FALSE, nodes=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indata}{
    A numeric vector whose elements are the empirical number of events
    (e.g. deaths). The length of the data is typically 12 (=months), 52
    or 53 (weeks), or 365 or 366 (days). Not providing 'indata'  is also
    okay (slightly different output then).
  }
  \item{long}{
    The basic length of the data analyzed, i.e. if we have monthly data,
    it would be 12 (hence it is an integer scalar). If 'indata' are
    provided, argument 'long' is  calculated based on argument 'indata'.
  }
  \item{peak}{
    Length of peak period (integer scalar). For instance, if we assume
    that the 'peak season' lasts six months for monthly data, 'peak'
    would be six (=default value).
  }
  \item{repts}{
    How many Monte Carlo simulation runs should be conducted (integer
    scalar)?
  }
  \item{whole.distribution}{
    Argument 'whole.distribution' indicates whether the whole
    distribution should be returned (=TRUE) or only the critical values
    (=FALSE) (Boolean Scalar).
  }
  \item{siglevels}{
    For which significance levels should the corresponding required rank
    sums be returned. Default settings are the 'typical' signficance
    levels of 0.001,0.01, 0.05, and 0.1 (numeric vector).
  }
  \item{PARALLEL}{
    If TRUE, multi-core Monte Carlo Simulation; otherwise single-core
    Simulation (Boolean Scalar). 
  }
  \item{nodes}{
    Specify on how many nodes the estimation should run
    (default=1). Only active if argument PARALLEL=TRUE.
  }
  
}
\value{
  \item{maximum.rank.sum}{The maximum rank sum theoretically possible
    with the given data (integer scalar)}
  \item{observed}{The observed maximum rank sum (with the given data)
    (integer scalar)}
  \item{observed.p.value}{What is the p-value corresponding to the
    observed maximum rank sum (numeric scalar)}
  \item{critical values}{What are the required rank sums for the entered
    significance levels (numeric data.frame)?}
  \item{distribution}{ONLY IF whole.distribution=TRUE: A numeric
    data.frame specifying all possible rank-sums and their associated
    p-values.}
}\references{
  EDWARDS, J.H. (1961):
  The recognition and estimation of cyclic trends.
  Annals of Human Genetics 25:83-86

  HEWITT, David and MILNER, Jean and CSIMA, Adele and PAKULA (1971):
  On Edwards' Criterion of Seasonality and a Non-Parametric Alternative.
  British Journal of Preventive Social Medicine 25:174-176

  ROGERSON, Peter A. (1996):
  A Generalization of Hewitt's Test for Seasonality.
  International Journal of Epidemiology 25:644-648

  WALTER, S.D. (1980):
  Exact significance levels for {H}ewitt's test for seasonality
  Journal of Epidemiology and Community Health 34:147--149
}

\author{
Roland Rau
}
\examples{
npst(indata=c(1:6, 12:7), peak=6, repts=100000,
  whole.distribution=FALSE, siglevels=c(0.01, 0.05, 0.1),
  PARALLEL=FALSE)
}
\keyword{nonparametric}
