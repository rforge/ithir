# Purpose        : Dissever: algorithm for spatial downscaling using GAMS 
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 
# Note           : Algorithm described in: [DOI: 10.1016/j.cageo.2011.08.021]



if(!isGeneric("dissever")){
  setGeneric("dissever", function(r2,c.grid, ...){standardGeneric("dissever")})}

setMethod('dissever', signature(r2 = "RasterStack",c.grid = "RasterLayer"), 
          function(r2,c.grid, ss= 0.5, cores= 8, thresh = 0.01){
            
            message("Dissever initialisation")
            diogRMSE<- matrix(NA,ncol=3,nrow=100)  
            beginCluster(cores)
            
            dir.create("disseverOuts",showWarnings = F)
            
            #get the cell numbers of the grid
            c.gridCELL<- writeRaster(rasterFromXYZ(cbind(xyFromCell(c.grid, cell=seq(1,ncell(c.grid)), spatial=FALSE),seq(1:ncell(c.grid)))),filename="disseverOuts/cellNumbers_coarse.tif",format="GTiff",progress=F,overwrite=T)
            
            
            #resample coarse grids to fine grid
            c.grid_ds<- resample(c.grid, r2[[1]],method="ngb", filename="disseverOuts/Initialisation_fineRes.tif",format="GTiff",progress=F,overwrite=T)
            c.gridCELL_ds<- resample(c.gridCELL, r2[[1]],method="ngb", filename="disseverOuts/cellNumbers_fineRes.tif",format="GTiff",progress=F,overwrite=T) 
            
            #Make a stack of the coarse grid prediction, cell numbers and covariates
            r1<- stack(c.grid_ds,c.gridCELL_ds)  
            
            
            #Model Frame
            gm1<- paste(names(r1)[1],paste(paste("~s(",names(r2)[1],sep=""),")",sep=""),sep=" ")
            for ( i in 2:length(names(r2))){
              gm1<-  paste(gm1,paste(paste("+s(",names(r2)[i],sep=""),")",sep=""),sep= " ")}
            
            #take a sample for modelling
            r3<- stack(r1,r2)
            sr<-sampleRandom(r3,ceiling(ncell(r3[[1]])*ss)) # sample random grid cells
            sr<- as.data.frame(sr)
            
            #Fit model (Initialisation)
            gm1 <- as.formula(unclass(gm1))
            fitModel<-gam(gm1,data=sr) 
            
            #predict model
            map1 <- clusterR(r2, predict, args=list(fitModel, type="response"),filename="disseverOuts/gamMap1.tif",format="GTiff",progress=F,overwrite=T)
            names(map1)<- "map1"
            
            #coarse grid table
            c.dat<- getValues(c.grid)
            c.datXY<- xyFromCell(c.grid, cell=seq(1,ncell(c.grid)), spatial=FALSE)
            c.dat<- cbind(seq(1,ncell(c.grid)),c.datXY,c.dat)
            c.dat<- as.data.frame(c.dat)
            names(c.dat)[1]<- "cell"
            c.dat_ref<- c.dat[which(complete.cases(c.dat)),] #complete cases  
            
            message("Dissever iteration")
            for (zz in 1:100){
              
              #stack predictions and cell numbers (fine grid)
              r4<- stack(c.gridCELL_ds,map1)
              r4.dat<-getValues(r4)
              r4.dat<- as.data.frame(r4.dat)
              r4.dat<- r4.dat[which(complete.cases(r4.dat)),] #complete cases
              downFit<-aggregate(r4.dat$map1,list(group=r4.dat$cellNumbers_coarse),mean) # average within coarse grid
              
              
              #merge
              xx<- merge(c.dat_ref, downFit, by.x = "cell", by.y = "group") 
              names(xx)<- c("cell", "X", "Y", "C_grid", "F_grid")
              
              #RMSE
              diogRMSE[zz,2]<- sqrt(mean((xx$C_grid -xx$F_grid)^2))
              varCalc<- ((1)^2/(nrow(xx)*(nrow(xx)-1)))* sum((xx$C_grid -xx$F_grid)^2)
              seVar<- sqrt(varCalc)
              tCrit<- qt(1-(0.05/2), df=nrow(xx)-1)
              cI<- seVar*tCrit
              upperCi<- (mean((xx$C_grid -xx$F_grid)^2))+ cI
              lowerCi<- (mean((xx$C_grid -xx$F_grid)^2))- cI
              diogRMSE[zz,1]<- sqrt(lowerCi)
              diogRMSE[zz,3]<- sqrt(upperCi)
              
              
              #calculate adjustment factor
              xx$AF<- -99999
              xx$AF[which(complete.cases(xx))]<- xx$C_grid/xx$F_grid 
              xx[xx == -99999] <- NA
              
              
              #make adjustment factor raster
              fs<- paste("disseverOuts/",paste("iter_",zz,sep=""),sep="")
              AF.grid<- writeRaster(rasterFromXYZ(xx[,c(2,3,6)]),filename=paste(fs,"AF1_coarse.tif",sep="_"),format="GTiff",progress=F,overwrite=T)
              
              #fine grid
              AF.grid_ds<- resample(AF.grid, r2[[1]],method="ngb", filename=paste(fs,"AF1_fine.tif",sep="_"),format="GTiff",progress=F,overwrite=T)
              #plot(AF.grid_ds)
              
              
              #do the adjustment
              r5<- stack(map1,AF.grid_ds)
              f1 <- function(x) (x[[1]]*x[[2]]) 
              upd.test <- clusterR(r5, fun=f1, filename=paste(fs,"update_fine.tif",sep="_"),format="GTiff",progress=F,overwrite=T)
              names(upd.test)<- names(c.grid_ds)
              #plot(upd.test)
              
              ###DO the modelling
              #take a sample for modelling
              r1<- stack(upd.test,c.gridCELL_ds)  
              r3<- stack(r1,r2)
              sr<-as.data.frame(sampleRandom(r3,ceiling(ncell(r3[[1]])*ss))) # sample random grid cells
              
              #Fit model (iteration)
              fitModel<-gam(gm1,data=sr) 
              
              #predict model
              map2 <- clusterR(r2, predict, args=list(fitModel, type="response"),filename=paste(fs,"pred_fine.tif",sep="_"),format="GTiff",progress=F,overwrite=T)
              
              #criteria                           
              message("RMSE = ",round(diogRMSE[zz,2],3))
              if (zz >= 5) {
                FF<- mean(abs(diogRMSE[zz-2,2]-diogRMSE[zz-1,2])+ abs(diogRMSE[zz-1,2]-diogRMSE[zz,2]) + abs(diogRMSE[zz-2,2]-diogRMSE[zz,2]))
                if (FF <= thresh) {break}} 
              map1<- map2
              names(map1)<- "map1"}
            diogRMSE<- diogRMSE[which(complete.cases(diogRMSE)),]
            diogRMSE<- as.data.frame(diogRMSE)
            names(diogRMSE)<- c("lowerCI_RMSE", "RMSE", "upperCI_RMSE")
            endCluster()
            #plot(map2)
            return(diogRMSE)})

#end script
