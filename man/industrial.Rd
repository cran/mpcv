\name{industrial}
\alias{industrial}
\docType{data}
\title{Industrial bivariate dataset}
\description{The dataset represents the measurements of an industrial product with two quality characteristics: Brinell hardness \code{H} and tensile strength \code{S}.}
\usage{data("industrial")}
\format{
  A list containing:
  \describe{
    \item{\code{x}}{a matrix with 25 observations and two quality characteristics: Brinell hardness \code{H} and tensile strength \code{S},}
	\item{\code{USL}}{the vector of the upper specification limits,}
	\item{\code{LSL}}{the vector of the lower specification limits,}
	\item{\code{Target}}{the vector of the target.}
  }
}
\references{
Sultan T. (1986) \emph{An acceptance chart for raw materials of two correlated properties}, Quality Assurance, Vol. 12, No. 3, 70-72. Quorum Books.
}

\examples{
data("industrial")
x <- industrial$x
LSL <- industrial$LSL
USL <- industrial$USL
Target <- industrial$Target
}

\keyword{datasets}