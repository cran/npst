\name{npst-package}
\alias{npst-package}
\docType{package}
\title{
npst
}
\description{
  Package 'npst' generalizes Hewitt's (1971) test for seasonality and
  Rogerson's (1996) extension based on Monte-Carlo simulation.
}
\details{
\tabular{ll}{
Package: \tab npst\cr
Type: \tab Package\cr
Version: \tab 1.6\cr
Date: \tab 2011-08-26\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Roland Rau

Maintainer: Roland Rau <roland.rau@gmail.com>
}
\references{
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
\keyword{ package }
\examples{
npst(indata=c(1:6, 12:7), peak=6, repts=100000,
  whole.distribution=FALSE, siglevels=c(0.01, 0.05, 0.1),
  PARALLEL=FALSE)
}
