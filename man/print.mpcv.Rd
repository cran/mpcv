\name{print}
\alias{print.mpcv}
\title{Prints an mpcv object}
\description{This function prints a \code{mpcv} object. It is a method for the generic function \code{print} of class \code{"mpcv"}.}
\usage{
\method{print}{mpcv}(x, ...)
}
\arguments{
  \item{x}{object of class \code{"mpcv"}.}
  \item{...}{other arguments, all currently ignored.}
}
 
\details{
This function prints values of the three components: \code{CpV}, \code{PS} and \code{PD} with names of variable which influence the components values the most.
}
\author{
Krzysztof Ciupke, \code{<krzysztof.ciupke at polsl.pl>}
}
\seealso{
 \code{\link[mpcv]{mpcv}}
}
\examples{
data(industrial)   
x <- industrial$x
LSL <- industrial$LSL
USL <- industrial$USL
Target<- industrial$Target
res.ind <- mpcv(x, LSL=LSL, USL=USL, Target=Target, alpha=0.025)
print(res.ind)
}
\keyword{multivariate}



