% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{feedlot}
\alias{feedlot}
\title{Feedlot data}
\format{A data frame with 67 observations and 4 variables:
\describe{
  \item{\code{herd}}{a factor with levels \code{9} \code{16} \code{3}
    \code{32} \code{24} \code{31} \code{19} \code{36} \code{34} \code{35}
    \code{33}, designating the herd that a feeder calf came from.}
  \item{\code{diet}}{a factor with levels \code{Low} \code{Medium}
    \code{High}: the energy level of the diet given the animal.}
  \item{\code{swt}}{a numeric vector: the weight of the animal at slaughter.}
  \item{\code{ewt}}{a numeric vector: the weight of the animal at entry to the feedlot.}
}}
\source{
Urquhart NS (1982) Adjustment in covariates when one factor affects
  the covariate. \emph{Biometrics} 38, 651-660.
}
\usage{
feedlot
}
\description{
This is an unbalanced analysis-of-covariance example, where one covariate is
affected by a factor. Feeder calves from various herds enter a feedlot, where
they are fed one of three diets. The weight of the animal at entry is the
covariate, and the weight at slaughter is the response.
}
\details{
The data arise from a Western Regional Research Project conducted at New
Mexico State University. Calves born in 1975 in commercial herds entered a
feedlot as yearlings. Both diets and herds are of interest as factors. The
covariate, \code{ewt}, is thought to be dependent on \code{herd} due to
different genetic backgrounds, breeding history, etc. The levels of
\code{herd} ordered to similarity of genetic background.

Note: There are some empty cells in the cross-classification of 
\code{herd} and \code{diet}.
}
\examples{
feedlot.lm <- lm(swt ~ ewt + herd*diet, data = feedlot)

# Obtain EMMs with a separate reference value of ewt for each 
# herd. This reproduces the last part of Table 2 in the reference
emmeans(feedlot.lm,  ~ diet | herd,  cov.reduce = ewt ~ herd)

}
\keyword{datasets}
