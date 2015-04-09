\name{Hunter Valley, NSW 100m digital elevation model}
\docType{data}
\alias{Hunter Valley DEM}
\title{Hunter Valley DEM}
\description{A \code{dataframe} that is easily converted to a grid raster. It has 3 columns, with the first two being spatial coordinates, and the third being ground elevation information. The coordinates are a regular grid point pattern with 100m spacing, and when converted to raster, resolve to be a digital elevation model for the Lower Hunter Valley region, NSW, Australia. The CRS of the coordinates is WGS 84 UTM Zone 56. }
\usage{data(HV_dem)}
\format{
\code{HV_dem} is a large \code{dataframe} with spatial coordinates, and elevation information labelled in respective columns. 
}
\details{This \code{dataframe} stores information that depicts a regular raster grid}

\references{
\itemize{
\item This data is copyright of the Soil Security Lab, The University of Sydney. If any part of this data is to be used in any publication or report, please provide a citation:
Soil Security Laboratory, 2015. Use R for Digital Soil Mapping Manual. The University of Sydney, Sydney, Australia
}
}

\examples{
library(ithir)
library(raster)


data(HV_dem)
map<- rasterFromXYZ(HV_dem)
plot(map) #display map


}
\keyword{datasets}
