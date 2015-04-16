\name{Edgeroi soil carbon data}
\docType{data}
\alias{edgeroi_splineCarbon}
\title{Soil point data from the Edgeroi District, NSW, Australia.}
\description{A soil information \code{dataframe} derived from 341 obsevations at various locations of the Edgeroi District, NSW, Australia (30.11S 149.66E). The soil data in this data set is soil carbon density and was estimated via a pedotransfer function that included measured soil attributes soil carbon concentration and soil texture information. Mass preserving splines were fitted to the original data in order to harmonise the depth intervals in accordances to specifications of the GlobalSoilMap project (Arrouays et al. 2014). In terms of the sampling locations, 210 were sampled on a systematic, equilateral triangular grid with a spacing of 2.8 km between sites (McGarry et al., 1989). The further 131 soil pro
files are distributed more irregularly or on transects. Locations are attributed to a site name as designated by McGarry et al. (1989) and have a recorded coordinate location. The CRS of the points is WGS84 UTM Zone 55.}
\usage{data(edgeroi_splineCarbon)}
\format{
\code{edgeroi_splineCarbon} is a 341 row \code{dataframe} with identifier, locational and soil carbon density information for the following depth intervals: 0-5cm, 5-15cm, 15-30cm, 30-60cm, 60-100cm and 100-200cm. A column which specifies the soil depth at each location irrespective of maximum spline fitted depth. 
}
\details{This data frame is a typical soil information table}

\references{
\itemize{
\item Arrouays, D., McKenzie, N., Hempel, J., Richer de Forges, A., and McBratney, A. (eds) (2014). GlobalSoilMap: Basis of the Global Spatial Soil Information System. CRC Press.
\item McGarry, D., Ward, W.T., McBratney, A.B. (1989) Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set. (2 vols) (CSIRO Division of Soils: Adelaide).
}
}

\examples{
library(ithir)
library(sp)

data(edgeroi_splineCarbon)

#plot the point locations
coordinates(edgeroi_splineCarbon)<- ~ east + north
plot(edgeroi_splineCarbon)

}
\keyword{datasets}
