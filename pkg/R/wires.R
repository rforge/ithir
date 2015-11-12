# Purpose        : Implementation of the WIRES algorithm (Wide-ranging exploratory survey) as described in Stockmann et al. (2015)
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 

#########
wires<- function(points, bb, nos.transects){
  
  #internal function
  # Selection of the uPoint 
  maxPoint<- function(polyDat, tPoint, sPoint){
    # select the polygons situated at tPoint
    np<-polyDat[polyDat$X==tPoint$X &polyDat$Y==tPoint$Y,] # select the polygons situated at tPoint
    pn<-polyDat[polyDat$ID %in% np[,1],] # All the possible neigbours of tPoint
    pn<-pn[!duplicated(pn[,2:3]),] # remove the duplicated entries
    pn<-pn[pn$X !=sPoint$X | pn$Y != sPoint$Y,] # remove the sPoint entry
    pn<-pn[pn$X !=tPoint$X | pn$Y != tPoint$Y,] # remove the tPoint entry
    
    #####plot diagnostics for the above pn object
    #points(pn[1,2],pn[1,3], pch=18, col="blue")
    #points(pn[2,2],pn[2,3], pch=18, col="red")
    #points(pn[3,2],pn[3,3], pch=18, col="blue")
    #points(pn[4,2],pn[4,3], pch=18, col="blue")
    #####
    
    # distance of sPoint to all entries of pn
    pn$su<-spDistsN1((as.matrix(pn[,2:3])), as.matrix(sPoint[,2:3]), longlat = FALSE)  # the distance from sPoint and possible vertex
    pn$st<-spDistsN1((as.matrix(tPoint[,2:3])), as.matrix(sPoint[,2:3]), longlat = FALSE)  # the distance from sPoint and possible vertex
    pn$tu<-spDistsN1((as.matrix(pn[,2:3])), as.matrix(tPoint[,2:3]), longlat = FALSE)  # the distance from sPoint and possible vertex
    
    # Derive angle
    #cos C    = (a2 + b2 ??? c2)/2ab 
    pn$deg<-0
    pn$rati<- 0
    for (j in 1: nrow(pn)){ 
      pn$rati[j]<- ((pn[j,]$st)^2 + (pn[j,]$tu)^2 - (pn[j,]$su)^2)/(2*(pn[j,]$st)*(pn[j,]$tu))
      pn$deg[j]<-acos(pmin(pmax(pn$rati[j],-1.0),1.0))*(180/pi)} 
    uuPoint<-pn[with(pn, order(deg)), ][nrow(pn),1:3] #get the nearest neighbours to sPoint
    return(uuPoint)}
  
  message("WIRES: Setting up the sampling frame")
  dir.create("wireOuts",showWarnings = F)  
  setwd(paste(getwd(),"/wireOuts",sep=""))
  W <- ripras(bb, shape="rectangle") # boundary window of the study area
  X <- as.ppp(points, W=W) # point pattern class
  Y <- delaunay(X) # perforn the delaunay triangulation
  #plot(Y, main="Delaunay Triangulation") # plot the delaunay triangulation
  
  
  # Write the delaunay triangulation image to a shape file
  Z <- as(Y, "SpatialPolygons")
  IDs <- sapply(slot(Z, "polygons"), function(x) slot(x, "ID"))
  df <- data.frame(rep(0, length(IDs)), row.names=IDs)
  SPDF <- SpatialPolygonsDataFrame(Z, df)
  tf <- tempfile()
  writePolyShape(SPDF, "sampleTriangulation")
  
  
  # prepartion for the sampling
  # PUT ALL THE POLYGON IDS AND ASSCIATED VERTICES INTO A SINGLE DATAFRAME
  # How big is the matrix we want to put all the coordinates and associated polygon labels
  ss<- 0
  for(i in 1: nrow(SPDF@data)){
    oo<-length(SPDF@polygons[[i]]@Polygons[[1]]@coords[,1])-1
    ss<-ss+oo}
  
  # Put all the polgon vertices into a single frame
  polyMat<- matrix(NA, nrow=ss, ncol=3) # matrix for all the ploygon vertices
  tt=1
  for(i in 1: nrow(SPDF@data)){
    pp<-SPDF@polygons[[i]]@Polygons[[1]]@coords
    polyMat[tt:(tt+2),2:3]<-pp[1:(nrow(pp)-1),]
    polyMat[tt:(tt+2),1]<-i
    tt=tt+(length(SPDF@polygons[[i]]@Polygons[[1]]@coords[,1])-1)}
  polyDat<- as.data.frame(polyMat)
  names(polyDat)<- c("ID", "X", "Y")
  
  ##### Convex Hull around points and remove them from the sampling dataframe (the points that make up the convex hull)
  rand.tr<-tri.mesh(points$X,points$Y)
  rand.ch<-convex.hull(rand.tr, plot.it=F) #convex hull
  CHpoints<-cbind(rand.ch$x,rand.ch$y) # Convex Hull Points
  CHpointsR<- rbind(CHpoints, CHpoints[1,])
  npt<- nrow(points)-nrow(CHpointsR) #number of possible transects that can be created
  transList<- list(CHpointsR)
  rr<-which(polyDat$X %in% CHpoints[,1] & polyDat$Y %in% CHpoints[,2])
  PolyDatr<- polyDat[-rr,] # vertex dataframe with convex hull points removed
  
  
  #plotting
  #plot(Y,col="white", main="sampling frame")
  #lines(c(rand.ch$x,rand.ch$x[1]), c(rand.ch$y,rand.ch$y[1]),col="red",lwd=1)
  
  ####Begin the sample...
  # From PolyDatr selct a point at random
  # some criteria regarding that you dont want points within some minimum distance of a convex hull point
  # Located the nearest point
  # From the nearest point located another 
  if (nos.transects > npt)
  {stop(paste("WIRES: Reduce the number of transects to less than or equal to:", npt,sep= " "))}
  
  
  message("WIRES: generating transects")
  for (i in 1:nos.transects){ # number of transects wanted
    tr<- i + 1
    transectPath<- matrix(NA, nrow=1000, ncol=3) # matrix for all the ploygon vertices
    rsp<-sample(nrow(PolyDatr), size=1) 
    sPoint<-PolyDatr[rsp,] # a random point selection FIRST POINT
    ww<-which(PolyDatr$X %in% sPoint[,2] & PolyDatr$Y %in% sPoint[,3])
    PolyDatr<- PolyDatr[-ww,] # vertex dataframe with sPoint points removed
    transectPath[1,]<- as.matrix(sPoint)
    #points(sPoint$X, sPoint$Y, pch=16, col="blue")
    
    #get the nearest neighbours to sPoint ie find tPoint
    nn<- spDistsN1((as.matrix(PolyDatr[,2:3])), as.matrix(sPoint[,2:3]), longlat = FALSE)  # the distance from function
    dfDATA<-cbind(PolyDatr,nn)
    dd<-dfDATA[with(dfDATA, order(nn)), ][1:100,] #get the nearest neighbours to sPoint
    tPoint<-dd[dd$nn != 0,][1,1:3] # SECOND POINT
    transectPath[2,]<- as.matrix(tPoint)
    #points(tPoint$X, tPoint$Y, pch=16, col="red")
    logt<- 0.5*(sum(CHpoints[,1]==tPoint[,2])+sum(CHpoints[,2]==tPoint[,3])) # see whether tPoint is a CHpoint
    if (logt >= 1)
    {i=i-1;stop}else {uPoint<- maxPoint(polyDat, tPoint, sPoint) # get the uPoint
                      logu<- 0.5*(sum(CHpoints[,1]==uPoint[,2])+sum(CHpoints[,2]==uPoint[,3])) # see whether uPoint is a CHpoint
                      if (logu >= 1)
                      {i=i-1;stop}else {#points(uPoint$X, uPoint$Y, pch=16, col="blue")
                        transectPath[3,]<- as.matrix(uPoint)
                        
                        # intertive section
                        dod<-4
                        repeat{
                          sPoint<- tPoint
                          tPoint<- uPoint
                          uPoint<- maxPoint(polyDat, tPoint, sPoint) # get the uPoint
                          logu<- 0.5*(sum(CHpoints[,1]==uPoint[,2])+sum(CHpoints[,2]==uPoint[,3])) # see whether uPoint is a CHpoint
                          if (logu >= 1)
                          {#points(uPoint$X, uPoint$Y, pch=16, col="red")
                           transectPath[dod,]<- as.matrix(uPoint)
                           transList[[tr]]<- transectPath[1:dod,]
                           #lines(transectPath[1:dod,2],transectPath[1:dod,3],col="green",lwd=1 )
                           break} else {#points(uPoint$X, uPoint$Y, pch=16, col="blue")
                             transectPath[dod,]<- as.matrix(uPoint)}
                          dod<-dod+1}}}}
  
  
  #Filter out overlapping transects
  message("WIRES: filtering out overlapping transects if any")
  nul.dat<- as.numeric(summary(transList)[,1])
  nul.dat.in<-which(nul.dat!=0)
  transList2<- transList[nul.dat.in]
  transList2<- transList2[2:length(transList2)] #remove hull points
  
  
  simDist<- matrix(NA, nrow=length(transList2), ncol=length(transList2)) # matrix to put transect distances
  for (i in 1:length(transList2)){
    sub1<-transList2[[i]] 
    for (j in 1:length(transList2)){
      comb1<- rbind(sub1,transList2[[j]] )
      simDist[j,i]<-sum(as.numeric(duplicated(comb1)))/nrow(sub1)}}
  diag(simDist)<- 0
  simDist.maxs<- colMaxs(simDist)
  selcted.transects<- which(!simDist.maxs>=0.95) # if 2 opposing transects are 95% similar
  transList3<- transList2[selcted.transects]
  
  
  
  
  # Make transects into shapefile
  # Make transect vertices into shapefile
  message("WIRES: Saving transects to shapefile objects")
  transectDist<- matrix(NA, nrow=length(transList3), ncol=2) # matrix to put transect distances
  for (i in 1:length(transList3)) {
    transectDist[i,1]<- i
    sdat<- transList3[[i]]
    if (length(sdat[,1])==0)
    {print(paste("null", i, sep="_"))} else {
      sdat<- as.data.frame(sdat)
      coordinates(sdat)<-  ~V2 + V3
      writePointsShape(sdat, paste("transect", i,"POINT" ,sep="_"))
      L = SpatialLines(list(Lines(list(Line(coordinates(sdat))),paste("transect", i ,sep="_"))))
      SLDF = SpatialLinesDataFrame(L, data.frame(Z = c("Line"), row.names = paste("transect", i ,sep="_")))
      transectDist[i,2]<-SpatialLinesLengths(SLDF)/1000
      writeLinesShape(SLDF, paste("transect", i,"LINE" ,sep="_"))}}
  
  transectDist<- as.data.frame(transectDist)
  names(transectDist)<- c("transect", "distance")
  write.table(transectDist,file="TransectDistances.txt",sep=",", col.names=T,row.names=F)
  
  #Save all transects to 1 shapefile
  l.list <- vector("list", length(transList3))
  for (i in seq_along(l.list)) {
    l.list[[i]] <- Lines(list(Line(coordinates(transList3[[i]][,2:3]))), ID= paste("transect",as.character(i),sep="_" ))}
  xxx<- SpatialLines(l.list)
  SLDF.c = SpatialLinesDataFrame(xxx,data=transectDist, match.ID=F)
  writeLinesShape(SLDF.c, paste("transects","ALL" ,sep="_"))
  message(paste(paste("WIRES: There were", length(transList3), sep=" "), "transects generated and saved", sep=" "))}


#END FUNCTION
