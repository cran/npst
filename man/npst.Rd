\name{npst}
\alias{npst}
\title{
  NPST - A Generalized Nonparametric Seasonality Test
}
\description{
  The function 'npst' generalizes the nonparametric seasonality tests of
  Hewitt et al. (1971) and of Rogerson (1996). At its core, the function
  npst calls a Fortran (90/95) subroutine to estimate via Monte Carlo
  simulation the empirical distribution function of the test statistic
  for the given input arguments.
}
\usage{
npst(indata = NULL, long = 12, peak = 6, repts = 1e+05,
whole.distribution = FALSE, siglevels = c(0.001, 0.01, 0.05, 0.1), seed
= NULL) 
}
\arguments{
  \item{indata}{
    \code{indata} needs to be a numeric vector whose elements are the
    empirical number of events (e.g. deaths). The length of the data is
    typically 12 (=months), 52 or 53 (weeks), or 365 or 366 (days). Not
    providing 'indata' is also okay (slightly different output then, see section'value').
  }
  \item{long}{
    \code{long} is a numeric scalar denoting the basic length of the data
    analyzed, i.e. if we have monthly data, it would be 12. If  \code{indata}
    are provided, argument \code{long} is calculated based on argument
    \code{indata}
  }
  \item{peak}{
    \code{peak} is a numeric scalar denoting the length of the peak
    period. For instance, if we assume that the 'peak season' lasts six
    months for monthly data, \code{peak} would be six (=default value)
  }
  \item{repts}{
    \code{repts}, a numeric scalar, denotes how many Monte Carlo
    simulation runs should be conducted. I would recommend at least 10 Mio.
  }
  \item{whole.distribution}{
    \code{whole.distribution}, a logical scalar,  indicates whether the
    whole distribution should be returned (=TRUE) or only the critical
    values (=FALSE) for the given significance levels (see argument \code{siglevels}).
  }
  \item{siglevels}{
    \code{siglevels}, a numeric vector, denotes for which significance
    levels should the corresponding required rank sums be returned. Default
    settings are the 'typical' signficance levels of 0.001,0.01, 0.05, and
    0.1
  }
  \item{seed}{
    With argument \code{seed}, a numeric vector, the user can supply a
    starting seed for the random number generator. Useful if results
    should be replicated. If no seed is provided, a random seed will be
    generated.
  }
}
\details{
  Please see the forthcoming(?) paper by the author of this package where
  he writes in detail about this function and how to use it in
  conjunction with an explanation of the background.
}
\value{
  \item{maximum.rank.sum}{The theoretical maximum rank sum; i.e. 57 in the
    "standard" case of a peak period of 6 months with monthly data (12+11+10+9+8+7=57)}
  \item{observed}{What is the maximum rank sum observed in the data?
    Will be only returned if \code{indata} have been provided.}
  \item{observed.p.value}{What is the 'p-Value' corresponding to the
    observed maximum rank sum (\code{observed}) (only if 'indata' have been provided)}
  \item{critical.values}{Which rank-sums are the threshold values around
    the significance levels (provided as input)?}
  \item{distribution}{A dataframe with two columns. Column 1 denotes
    *all* possible ranksums. Column 2 denotes the respective values of
    the empirical distribution function (only if input argument
    whole.distribution=TRUE; default is FALSE).}
}
\references{
  D. Hewitt, J. Milner, A. Csima, A. Pakula (1971): On Edwards'
  Criterion of Seasonality and a Non-Parametric Alternative. British
  Journal of Preventive Social Medicine 25:174--176
      
  P.A. Rogerson (1996): A Generalization of Hewitt's Test for
  Seasonality. International Journal of Epidemiology 25:644--648

  S.D. Walter (1980): Exact significance levels for Hewitt's test for
  seasonality. Journal of Epidemiology and Community Health
  34:147--149
}
\author{
Roland Rau
}

\examples{
npst()
}
\keyword{statistical test}
\keyword{seasonality}
