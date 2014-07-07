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
     graphic = FALSE, coef.up, coef.lo, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric matrix or data.frame containing the data (quality characteristics).
}
  \item{indepvar}{
a number or a name of the independent variable for building one-sided models.
}
  \item{LSL}{
a vector of lower specification limits defined for each variable.
}
  \item{USL}{
a vector of upper specification limits.
}
  \item{Target}{
vector of target of the process.
}
  \item{alpha}{
the proportion of nonconforming products.
}
  \item{distance}{
the distance measure to be used for removing the nonconforming elements. This must be one of \code{"mahalanobis"} (default), \code{"euclidean", "maximum", "manhattan", "canberra", "binary"} or \code{"minkowski"}. Any unambiguous substring can be given.
  
}
\item{n.integr}{
a number of integration intervals
}
  \item{graphic}{
logical value indicating whether the obtained results (process region) should be presented graphically (default \code{FALSE}).
}
  \item{coef.up}{
vector of minimal values of leading coefficients for "upper" one-sided models. Value given for \code{indepvar} is omitted (could be \code{NA}).
}
  \item{coef.lo}{
vector of minimal values of leading coefficients for "lower" one-sided models. Value for \code{indepvar} is omitted (could be \code{NA}).
}
   \item{verbose}{
logical value indicating whether the values of \code{Cpv}, \code{PS} and \code{PD} should be displayed (default \code{FALSE}).
}
}
\details{
If the parameter \code{Target} is not specified, then \code{Target <- LSL + (USL - LSL)/2}.

Shape of the process region is mainly defined by the leading coefficients. 
To obtain a similar shape of the process region in the future process analysis, there is possible to specify minimal values of the leading coefficients \code{coef.lo} and \code{coef.up}.

Except the \code{"mahalanobis"} distance, the available distance measures are listed in \link[stats]{dist}.
}
\value{
A list containing the following elements:
 \item{CpV}{the percentage value of the capability measurement component.}
 \item{PS}{the percentage value of the process shift component.}
 \item{PSvar}{the name of the variable which influences the process shift the most}
 \item{PD}{percentage value of the process distance component.}
 \item{PDvar}{the name of the variable which has the most negative influence on the value of the process distance component}
 \item{coef.lo}{a named vector of leading coefficients of "lower" one-sided models; for the independent variable \code{NA} is returned.}
 \item{coef.up}{a named vector of leading coefficients of "upper" one-sided models; for the independent variable \code{NA} is returned.}
}
\references{
Ciupke K. (2014) \emph{Multivariate Process Capability Vector Based on One-Sided Model}, Quality and Reliability Engineering International, John Wiley & Sons.
}
\author{
Krzysztof Ciupke
}

%%\note{
%%Informacje dodatkowe
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
 %%\code{"\link[MPCI]{mpci}"}
%% }
\examples{
%%\dontrun{    
	
data(industrial)   
x <- industrial$x
LSL <- industrial$LSL
USL <- industrial$USL
Target<- industrial$Target
res <- mpcv(x, LSL=LSL, USL=USL, Target=Target, alpha=0.025, graphic=TRUE)

data(automotive)
x <- automotive$x
LSL <- automotive$LSL
USL <- automotive$USL
Target<- automotive$Target
res <- mpcv(x, indepvar="T",  LSL=LSL, USL=USL, Target=Target, verbose=TRUE)

data(sleeves)
x <- sleeves$x
LSL <- sleeves$LSL
USL <- sleeves$USL
Target<- sleeves$Target
res <- mpcv(x, indepvar=3, LSL=LSL, USL=USL, Target=Target, alpha=.02, graphic=TRUE, verbose=TRUE)

%%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{models}



