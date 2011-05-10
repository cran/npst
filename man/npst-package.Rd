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
Version: \tab 1.2\cr
Date: \tab 2011-05-10\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Roland Rau

Maintainer: Roland Rau <roland.rau@gmail.com>
}
\references{
@article{edwards,
	author={Edwards, J. H. },
	title={The recognition and estimation of cyclic trends},
	journal={Annals of Human Genetics},
	volume={25},
	year={1961},
	pages={83--86}
}

@article{hewitt,
	author={Hewitt, David and Milner, Jean and Csima, Adele and Pakula, Andrew},
	title={On {E}dwards' Criterion of Seasonality and a Non-Parametric Alternative},
	journal={British Journal of Preventive Social Medicine},
	volume={25},
	year={1971},
	pages={174--176}
}

@article{rogerson,
	author={Rogerson, Peter A.},
	title={A {G}eneralization of {H}ewitt's {T}est for {S}easonality},
	journal={International Journal of Epidemiology},
	volume={25},
	year={1996},
	pages={644--648}
}


@article{walter,
	author={Walter, S.D.},
	title={Exact significance levels for {H}ewitt's test for seasonality},
	journal={Journal of Epidemiology and Community Health},
	volume={34},
	year={1980},
	pages={147--149}
}
  
  
}
\keyword{ package }
\seealso{
}
\examples{
###  replicate Walter's (1980) exact results for Hewitt's test for
###  seasonality via Monte Carlo Simulation
npst(long=12, peak=6, repts=100000, outputrows=5)
}
