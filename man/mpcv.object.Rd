\name{mpcv.object}
\alias{mpcv.object}
\title{MPCV (multivariate process capability vector) object}
\description{Structure of mpcv object}
\value{
An mpcv.object is a list containing the following elements:
\item{CpV}{the percentage value of the capability measurement component.}
\item{PS}{the percentage value of the process shift component.}
\item{PSvar}{the name of the variable which influences the process shift the most.}
\item{PD}{the percentage value of the process distance component.}
\item{PDvar}{the name of the variable which has the most negative influence on the value of the process distance component.}
\item{coef.lo}{a named vector of leading coefficients of "lower" one-sided models; for the independent variable \code{NA} is returned.}
\item{coef.up}{a named vector of leading coefficients of "upper" one-sided models; for the independent variable \code{NA} is returned.}
}
\author{
Krzysztof Ciupke, \code{<krzysztof.ciupke at polsl.pl>}
}
\seealso{\code{\link{mpcv}}}
