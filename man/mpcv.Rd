\name{mpcv}
\alias{mpcv}
\title{
 Multivariate process capability vector
}
\description{
 Performs the multivariate process capability analysis using three component multivariate process capability vector (mpcv).
}
\usage{
mpcv(x, indepvar = 1, LSL, USL, Target, alpha = 0.0027, distance, n.integr = 100, 
     coef.up, coef.lo)
}

\arguments{
  \item{x}{
a numeric matrix containing the data (quality characteristics).
}
  \item{indepvar}{
a number or a name of the independent variable needed for building one-sided models.
}
  \item{LSL}{
a vector of lower specification limits defined for each variable.
}
  \item{USL}{
a vector of upper specification limits defined for each variable.
}
  \item{Target}{
a vector of target of the process defined for each variable.
}
  \item{alpha}{
the proportion of nonconforming products.
}
  \item{distance}{
the distance measure to be used for removing the nonconforming elements. This must be one of \code{"mahalanobis"} (default), \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} or \code{"minkowski"}. Any unambiguous substring can be given.
  
}
  \item{n.integr}{
a number of integration intervals
}
  \item{coef.up}{
a vector of minimal values of leading coefficients for "upper" one-sided models. Value given for \code{indepvar} is omitted (could be \code{NA}).
}
  \item{coef.lo}{
a vector of minimal values of leading coefficients for "lower" one-sided models. Value for \code{indepvar} is omitted (could be \code{NA}).
}
}
\details{
If the parameter \code{Target} is not specified, then \code{Target <- LSL + (USL - LSL)/2}.

Using the applied methodology, the shape of a process region is mainly defined by the leading coefficients of the models. 
To obtain a certain shape of a process region (e.g. similar to the previous one) there is possible to specify minimal values of the leading coefficients \code{coef.lo} and \code{coef.up} of the models.
By default all the minimal values of the coefficients equal zero.

Except the \code{"mahalanobis"} distance, the available distance measures are listed in \link[stats]{dist}.
}
\value{
An mpcv object. See \code{\link{mpcv.object}} for details.
}
\references{
Ciupke K. (2014) \emph{Multivariate Process Capability Vector Based on One-Sided Model}, Quality and Reliability Engineering International, John Wiley & Sons.
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

data(automotive)
x <- automotive$x
LSL <- automotive$LSL
USL <- automotive$USL
Target<- automotive$Target
res.aut <- mpcv(x, indepvar="T",  LSL=LSL, USL=USL, Target=Target)

data(sleeves)
x <- sleeves$x
LSL <- sleeves$LSL
USL <- sleeves$USL
Target<- sleeves$Target
res.sle <- mpcv(x, indepvar=3, LSL=LSL, USL=USL, Target=Target, alpha=.02)
}

\keyword{multivariate}
\keyword{models}



