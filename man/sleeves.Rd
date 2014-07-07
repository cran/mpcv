\name{sleeves}
\alias{sleeves}
\docType{data}
\title{Dataset describing diameters of cylindrical sleeves}
\description{Dataset containing measurements of three identifiable diameters of cylindrical sleeves referred to as \code{A}, \code{B} and \code{C}.}
\usage{data("sleeves")}
\format{
  A list containing
  \describe{
    \item{\code{x}}{a data.frame with 28 observations and three quality characteristics \code{A}, \code{B} and \code{C},}
	\item{\code{USL}}{the vector of the upper specification limits,}
	\item{\code{LSL}}{the vector of the lower specification limits,}
	\item{\code{Target}}{the vector of the target.}
  }
}
\references{
Raissi S. (2009) \emph{Multivariate process capability indices on the presence of priority for quality characteristics}, Journal of Industrial Engineering International,  Vol. 5, No. 9, 27-36.
}
\examples{
	data("sleeves")
	x <- sleeves$x
	LSL <- sleeves$LSL
	USL <- sleeves$USL
	Target <- sleeves$Target
}

\keyword{datasets}