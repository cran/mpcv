\name{automotive}
\alias{automotive}
\docType{data}
\title{Automotive bivariate dataset}
\description{The dataset describing the problem of automatic screwing the car wheels. Two characteristics are observed: the torque \code{T} of tightening a screw, and the rotation angle \code{A} of the screw until the necessary value of the torque is acquired.}
\usage{data("automotive")}
\format{
  A list containing:
  \describe{
    \item{\code{x}}{a matrix with 47 observations and two quality characteristics \code{T} and \code{A},}
	\item{\code{USL}}{the vector of the upper specification limits,}
	\item{\code{LSL}}{the vector of the lower specification limits,}
	\item{\code{Target}}{the vector of the target.}
  }
}
\references{
Ciupke K. (2014)
}
\examples{
data("automotive")
x <- automotive$x
LSL <- automotive$LSL
USL <- automotive$USL
Target <- automotive$Target
}

\keyword{datasets}