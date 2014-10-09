\name{plot.mpcv}
\alias{plot.mpcv}
\title{Plots an mpcv object.}
\description{Plots a \code{mpcv} object on the current graphics device. This function is a method for the generic function \code{plot}, for objects of class \code{"mpcv"}.}
\usage{
\method{plot}{mpcv}(x, ...)
}
\arguments{
  \item{x}{object of a class \code{"mpcv"}.}
  \item{...}{other argument, currently no used.}
}
  
\details{
This function presents graphically the dataset, the given tolerance region (specification limits), the identified process region, the target and the marginal median for all pairs of a process quality characteristics (variables).
}

%%\note{
%%Informacje dodatkowe
%%}

%% ~Make other sections like Warning with \section{Warning }{....}

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
res.ind <- mpcv(x, LSL=LSL, USL=USL, Target=Target, alpha=0.025, coef.lo=c(NA,.005))
plot(res.ind)
}

\keyword{multivariate}


