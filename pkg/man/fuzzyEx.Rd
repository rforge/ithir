\name{fuzzyEx}
\docType{methods}
\alias{fuzzyEx}
\title{Derivation of fuzzy membership to classes}
\description{This is a simple function that complements the outputs from the Fuzme software developed by Minasny and McBratney (2002) and specifically upon outputs of the fuzzy kmeans with extragrades algorithm from McBratney and de Gruijter (1992). The function, given some inputs of multivariate data, together with a centroid table of the same multivariate information, will estimate the membership or belongingness of each row to each class. As the estimation is based on the fuzzy kmeans with extragrades algorithm, the membership are always derived for 1 + n classes i.e. a membership to each candidate class as defined in centroidal terms of the centroid table and a membership to the to extragrade class.}

\usage{fuzzyEx(data,centroid,cv,expon,alfa)
}
\arguments{
  \item{data}{A data frame where the first row is a row or observation identifier. The remaining columns hold data relating to the centroid table.}
  \item{centroid}{A matrix of the class centroids} 
  \item{cv}{A variance-covariance matrix used for the estimation of Mahalinobis distance.}
  \item{expon}{numeric. The fuzzy exponent value}
  \item{alfa}{numeric. a value indicating the level of membership to extragrade class. This value is an output from fuzme.}
}
\value{Re tuns a matrix with the same number ow rows in the \code{data} input and the number of columns equal to 2 + n where n is the number of classes. One of the extra columns contains the memberships to the extragrade class. The last remaining column information about which class the row has the highest membership to.}

\note{The distance measure for evaluating the difference between class centroids and observation is the Mahalinobis distance}

\author{Brendan Malone}

\references{
\itemize{
  \item McBratney, A.B., de Gruijter, J.J., (1992) \href{http://dx.doi.org/10.1111/j.1365-2389.1992.tb00127.x}{Continuum approach to soil classification by modified fuzzy k-means with extragrades}. Journal of Soil Science, 43:159-175.  
  \item Minasny, B., McBratney, A.B., (2002) \href{http://sydney.edu.au/agriculture/pal/software/fuzme.shtml}{FuzMe version 3.0}. Australian Centre for Precision Agriculture, The University of Sydney, Australia.
}
}


\examples{
library(raster)
data(edgeTarget_C)
target
bbRaster(target)
}
\keyword{methods}
