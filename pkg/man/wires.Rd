\name{wires algorithm}
\docType{methods}
\alias{wires}
\title{WIRES (wide-ranging exploratory survey)}
\description{This function is the implementation of the WIRES algorithm as described in Stockmann et al. (2015). WIRES stands for Wide-ranging exploratory survey. The algorithm was designed for the purposes of regional scale proximal sensing (soil, geophysical etc.) surveys where aerial survey is prohibitive. Due to the regional extent of survey, tightly-spaced drive lines is also prohibitive. The idea is to first place sampling points throughout the survey area using an established sampling design. Then transects are establshed sequentially by selecting a point at random then connecting it, and to the subsequent points until the transect terminates at the peripheral of a given survey area. With a user-defined number of transects (not greater than the number of sampling points), the effect will be for the study area, an arrangment of transects (wires), meandering and interesecting throughout.}

\usage{wires(points, bb, nos.transects)}
\arguments{
  \item{points}{dataframe; the first two columns of which are the coordinate locations of the sampling points.} 
  \item{bb}{matrix; specifically a 4  x 2 matrix with each row indicating a coordinate pair of the bounding box of the survey area for which the survey is to be completed. Must bound completely the sample points.}     
  \item{nos.transects}{numeric; the number of transects that are needed. Number of transects can not be greater than the number of sampling points. }
}
\value{A new folder is created in the working directory called \code{wireOuts} and will contain the following outputs:
\enumerate{
\item{A shapefile of a Delaunay triangualtion made from the given sample points called \code{"sampleTriangulation.shp"}}
\item{Shapefile objects in both point and line format of each created transect. For example \code{"transect_1_LINE.shp"} and \code{"transect_1_POINT.shp"} are the line and point shapefile objects of the first created transect.}
\item{A text file summarising the length of each transect called \code{"TransectDistances.txt"}}  
}}

\note{The WIRES algorithm follows the general format:
\enumerate{
\item {Make a configuration of \code{x} number of points across the land area to be surveyed. The selection of points could follow any type of sample design. For example, simple random or stratified simple random designs, or even conditioned Latin Hypercube sampling etc.}
\item {Form a Delaunay triangulation of these points.}
\item {Select a point at random without replacement. Call this point \code{A}. Points on the convex hull of the configuration are however not permitted to be sampled, as these form the termination sites of a formed transect.}
\item {Determine the nearest point to \code{A}. Call the nearest point \code{B}.}
\item {Form a triangle \code{ABC} with the immediate neighbouring points of \code{B}.}
\item {Point \code{C} is the neighbouring point of \code{B} which maximises the angle at \code{B}. The angles are measured using measured distances and cosine rules.}
\item {Transect segment is drawn following the triangulation lines from point \code{A}, through \code{B}, to point \code{C}. Point \code{B} is now labelled new point \code{A}, and point \code{C} is  now labelled new point \code{B}.}
\item {Go back to step 5, and repeat until point \code{C} becomes one of the termination sites, which marks the end of a transect. }
\item {Go back to step 3 to generate more transects.}
}
There is a maximum to how many transects can be created for a given number of point locations. This number is generally the number of sampling points minus the number of convex hull points. The fuction will give warning if too manny transects have been requested. 
After all transects have been formed, a filtering step is performed to remove redundant transect (those which are exactly similar to another transect). Subsequently the final number of transects created may not equal the number requested. This is not a problem however becasue in practice, only a selection of the generated transects may be used anyhows.}

\author{Brendan Malone}
\references{
\itemize{
  \item Stockmann, U., Malone, B.P., McBratney, A.B., Minasny, B., (2015) \href{http://dx.doi.org/10.1016/j.geoderma.2014.10.005}{Landscape-scale exploratory radiometric mapping using proximal
soil sensing}. Geoderma, 239-240: 115-129.}
}

\examples{
library(ithir)
library(raster)

#random points
data(hvPoints250)

#bounding box
data(hvGrid25m)
bb<- bbRaster(hvGrid25m)

#wires algorithm
wires(points = hvPoints250, bb = bb, nos.transects=10)
}
\keyword{methods}
