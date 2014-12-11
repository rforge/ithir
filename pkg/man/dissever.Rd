\name{dissever}
\docType{methods}
\alias{dissever}
\alias{dissever,RasterStack, RasterLayer-method}
\title{Spatial downscaling using generalised additive models}
\description{Algorithm performs spatial downscaling of coarse grid mapping to fine grid mapping using predictive covariates and a GAM (generalised additive model)}
\usage{
\S4method{dissever}{RasterStack, RasterLayer}(r2,c.grid, ss= 0.5, cores= 8)
}
\arguments{
  \item{r2}{object of class \code{"RasterStack"}: Raster stack of predictive covariates} 
  \item{c.grid}{object of class \code{"RasterLayer"}: The map that needs to be downscaled}   
  \item{ss}{numeric: The proportion of the map that is sampled for fitting the GAM model (a number between 0 and 1)}
  \item{cores}{numeric: number of computer cores to call upon for running the algorithm}
  \item{thresh}{numeric: dissever iterations will proceed until this value is met}
}
\value{Returns a data frame, and sequentially generates gridded maps as the algoithm goes through the initialisation and iteration phases:
  \describe{
   \item{diogRMSE}{data.frame: The RMSE is a measure of the difference between the map to be downscaled and the aggregated predictions of the downscaled mapping. Produced upon every iteration}
}
}
\note{dissever is a downscaling algorithm that seeks to map a target variable at fine resoultion given an available coarse resolution map of the same target variable, together with fine scaled predictive covariates that will inform the spatial prediction. The agorithm uses generalised additive modelling to fit the relationship betwween target variable and predictive covariates. The algorithm has two phases. 1. Initialisation - a coarse grid to fine grid resample followed by sampling and initial model fit. 2. Iteration. Adjustments are made to the predictions iteratively in order to try and ensure that the coarse grid map is linearly related to the fine grid predictions i.e. there is a mass balance property to be attained. Iterations proceed until the threshold is met. The threshold is based on the change in RMSE over three iterations. If the change is less than the threshold the algorithm will terminate. It is always assumed that the map extent of the coarse grid map to be downscaled is from the same area and same extent as the fine gridded covariates. A number of grids are produced from the algorithm, all of which are deposited into the disseverOuts folder which will be automatically created in your working directory after starting the downscaling
This version of dissever has not facilitated the ability to incorporate the prediction uncertainties of the coarse map}

\author{Brendan Malone}
\references{
\itemize{
  \item Malone, B.P, McBratney, A.B., Minasny, B., Wheeler, I., (2011) \href{http://dx.doi.org/10.1016/j.cageo.2011.08.021}{A general method for downscaling earth resource information}. Computers & Geosciences, 41: 119-125. 
}
}

\examples{
library(gam)
library(snow)
library(rgdal)
library(raster)
setwd("C:/temp")

#Map to be downscaled
data(edgeTarget_C)

#Predictive Grids
data(edgeGrids)

#Run dissever
dissever(r2=edgeGrids, c.grid=target, ss=0.75, cores= 2)
}
\keyword{methods}
