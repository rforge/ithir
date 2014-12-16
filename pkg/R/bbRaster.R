# Purpose        : returns the bounding box of a raster layer
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 

bbRaster<-function(obj){
  xxx<-extent(obj)
  d.mat<- matrix(NA, nrow=4, ncol=2)
  d.mat[1,]<- c(xxx[1],xxx[3])
  d.mat[2,]<- c(xxx[1],xxx[4])
  d.mat[3,]<- c(xxx[2],xxx[3])
  d.mat[4,]<- c(xxx[2],xxx[4])
  return(d.mat)}
