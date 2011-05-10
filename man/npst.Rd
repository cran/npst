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
npst(long = 12, peak = 6, repts = 1e+05, outputrows = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{long}{
    Argument \code{long} specifies how long the base period is. For
    instance, if we analyze months in a year, \code{long}'s value would
    be 12 (which is the default value to replicate Hewitt's results).
}
  \item{peak}{
    Argument \code{peak} specifies how long the peak period is. For
    instance, if we analyze months in a year, and we want to know
    whether there is peak period of 6 months, \code{peak}'s value would
    be 6 (which is the default value to replicate Hewitt's results).
}
  \item{repts}{
    Argument \code{repts} specifies how many Monte Carlo repetitions
    should be performed. The default value is 100,000. I recommend,
    however, at least 1,000,000 -- better even 10,000,000
}
  \item{outputrows}{
    
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
Roland Rau
}
\examples{
}
\keyword{nonparametric}

