\name{coef.mpcv}
\alias{coef.mpcv}
\title{Extracts leading coefficients of the one-sided models from the mpcv object}
\description{\code{coef} is a generic function which extracts the leading coefficients of the model of a process region; the coefficients are extracted from the object of class \code{"mpcv"}.}
\usage{
\method{coef}{mpcv}(object, ...)
}
\arguments{
  \item{object}{object of a class \code{"mpcv"}.}
  \item{...}{other arguments, currently not used.}
}
  
\details{
Shape of the process region is mainly defined by the leading coefficients of the process models (here one-sided models ale applied). To keep a similar shape of the process region in the future process analysis, minimal values of the leading coefficients are required. This function allows to extract the leading coefficients from the \code{object}.
}

\value{
A matrix with 2 rows (named \code{coef.lo} and \code{coef.up}) and number of columns corresponding the number of variables.
}

\note{
For the viariable used as the independent one in the \code{\link{mpcv}} function the \code{NA} value is returned.
}

\references{
Ciupke K. (2014) \emph{Multivariate Process Capability Vector Based on One-Sided Model}, Quality and Reliability Engineering International, John Wiley & Sons.
}
\seealso{
 \code{\link[mpcv]{mpcv}}
}

\author{
Krzysztof Ciupke, \code{<krzysztof.ciupke at polsl.pl>}
}


\examples{
data(industrial)   
x <- industrial$x
LSL <- industrial$LSL
USL <- industrial$USL
Target<- industrial$Target
res.ind <- mpcv(x, LSL=LSL, USL=USL, Target=Target, alpha=0.025)
coef(res.ind)
}
\keyword{multivariate}


