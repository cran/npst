## ################################################################################
## Program : NPST
## Language: R
## This code is the wrapper function to call the shared library
## "npstfortran". The final result derives the empirical
## distribution of a generalization of Hewitt's and Rogerson's
## nonparametric seasonality tests.
##
## Copyright (C) 2004-2010  Roland Rau
## 
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## ################################################################################

.First.lib <- function(lib, pkg)
  library.dynam("npst", pkg, lib)
     

npst <- function(long=12, peak=6, repts=100000, outputrows=10) {
  

  # what is the maximum rank sum?
  maxranksum <- sum(long:(long-peak+1))
  
  # initialize vector where results are stored with zeros.
  resvec <- rep(0, maxranksum)

  # initialize base data and repeat the values.
  basedata <- 1:long
  basedata2 <- rep(1:long,2)

  # calling external Fortran (90/95) routine to do the heavy work.
  outsourcedloop <- .Fortran("npst", PACKAGE="npst",
                             lng=as.integer(long),
                             pk=as.integer(peak),
                             rpts=as.integer(repts),
                             mxrksum=as.integer(maxranksum),
                             resvec=as.integer(resvec),
                             basedata=as.integer(basedata),
                             basedata2=as.integer(basedata2))

  # Calculate the p-values, i.e. the cumulative distribution function.
  pvals <- cumsum(outsourcedloop$resvec[order(maxranksum:1)])/repts
  maxisums <- maxranksum:1
  returndf <- data.frame(Max.Rank.Sum=maxisums, p.values=pvals)
  return(returndf[1:outputrows,])
}
