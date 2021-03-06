% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NPSimulation.R
\name{Calc4GroupNPStats}
\alias{Calc4GroupNPStats}
\title{Calc4GroupNPStats}
\usage{
Calc4GroupNPStats(x1, x2, x3, x4, sigfig = -1, alpha = 0.05)
}
\arguments{
\item{x1}{is the data associated with treatment A in one block 1}

\item{x2}{is the data associated with treatment B in block 1}

\item{x3}{is the data associated with treatment A in block 2}

\item{x4}{is the data associated with treatment B in block 2}

\item{sigfig}{is the number of significant digits in the data. If >0 the datav will be appropriately truncated}

\item{alpha}{is the significance level for all statistical tests}
}
\value{
The function returns the point biserial version of Kendall's tau and its variance, Cliff's d and its variance, the probability of superiority, phat, and its variance, for the 4 group experiment experiment.
}
\description{
This function does a non-parametric analysis of a randomized blocks experiment assuming 2 blocks and 2 treatment conditions.
}
\examples{
set.seed(123)
x=list()
x[[1]]=stats::rnorm(10,0,1)
x[[2]]=stats::rnorm(10,0.8,1)
x[[3]]=stats::rnorm(10,0.5,1)
x[[4]]=stats::rnorm(10,1.3,1)
Calc4GroupNPStats(x[[1]],x[[2]],x[[3]],x[[4]],sigfig=-1,alpha=0.05)
# A tibble: 1 x 17
#       N  phat phat.var phat.df phat.test phat.pvalue phat.sig      d   vard d.sig    cor   sqse
# <int> <dbl>    <dbl>   <dbl>     <dbl>       <dbl> <lgl>     <dbl>  <dbl> <lgl>  <dbl>  <dbl>
#   1    40  0.17  0.00497    31.0     -4.68   0.0000532 TRUE     -0.660 0.0206 TRUE  -0.347 0.0132
# … with 5 more variables: ctvar <dbl>, n1 <int>, n2 <int>, sigCVt <lgl>, sigCVn <lgl>
}
\author{
Barbara Kitchenham and Lech Madeyski
}
