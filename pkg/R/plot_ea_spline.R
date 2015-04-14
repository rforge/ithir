# Purpose        : Plot raw soil profile data
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 
# Contributions  : 
# Status         : working
# Note           : 


##### Function for plotting outputs from ea_spline
#plotting outputs from ea_spline
plot_ea_spline<- function(splineOuts, d = t(c(0,5,15,30,60,100,200)), maxd, type = 1, label = "", plot.which = 1){
  #type 1 (raw data and spline fit)
  if (type==1)
  { #plot the observed
    d1<- splineOuts$obs.preds 
    d1<- d1[d1$FID==plot.which,]
    vals<- d1[,4]
    depths<- d1[,2:3]
    matX<- matrix(NA, nrow=nrow(d1), ncol= 4)
    matY<- matrix(NA, nrow=nrow(d1), ncol= 4)
    for (i in 1: nrow(d1)){
      matX[i,]<- c(vals[i]-vals[i], vals[i], vals[i], vals[i]-vals[i]) 
      matY[i,]<- c(depths[i,1], depths[i,1], depths[i,2], depths[i,2]) }
    raw.plot<-plot(matX[1,], matY[1,],type="n",ylim=c(maxd,0),xlim=c(min(vals), (max(vals)*1.1)),main=paste("soil profile:",d1[1,1], sep=" "), ylab="depth", xlab= label, lty=2, lwd=1.5, xaxs="i", col="black", font.lab=2,cex.lab=1.5)
    for (i in 1: nrow(d1)){
      polygon (matX[i,],matY[i,], lty=1, lwd=2, border="black") }
    #plot the fitted spline
    lines(splineOuts$var.1cm[,plot.which],seq(1,d[length(d)]),lwd=2,col="red" )}
  
  #type 2 (raw data and spline fitted averages)
  if (type==2)
  {#plot the observed
    d1<- splineOuts$obs.preds 
    d1<- d1[d1$FID==plot.which,]
    vals<- d1[,4]
    depths<- d1[,2:3]
    matX<- matrix(NA, nrow=nrow(d1), ncol= 4)
    matY<- matrix(NA, nrow=nrow(d1), ncol= 4)
    for (i in 1: nrow(d1)){
      matX[i,]<- c(vals[i]-vals[i], vals[i], vals[i], vals[i]-vals[i]) 
      matY[i,]<- c(depths[i,1], depths[i,1], depths[i,2], depths[i,2]) }
    raw.plot<-plot(matX[1,], matY[1,],type="n",ylim=c(maxd,0),xlim=c(min(vals), (max(vals)*1.1)),main=paste("soil profile:",d1[1,1], sep=" "), ylab="depth", xlab= label, lty=2, lwd=1.5, xaxs="i", col="black", font.lab=2,cex.lab=1.5)
    
    d2<- as.matrix(splineOuts$harmonised[plot.which,2:length(d)])
    matX1<- matrix(NA, nrow=ncol(d2), ncol= 4)
    matY1<- matrix(NA, nrow=ncol(d2), ncol= 4)
    for (i in 1: ncol(d2)){
      matX1[i,]<- c(d2[i]-d2[i], d2[i], d2[i], d2[i]-d2[i]) 
      matY1[i,]<- c(d[1,i], d[1,i], d[1,i+1], d[1,i+1]) }
    #plot the spline averages
    for (i in 1: ncol(d2)){
      polygon (matX1[i,],matY1[i,], lty=1, lwd=1,col="green", border="green") }
    
    for (i in 1: nrow(d1)){
      polygon (matX[i,],matY[i,], lty=1, lwd=2, border="black") }}
  
  #type 3 (raw data and spline fitted averages)
  if (type==3)
  {#plot the observed
    d1<- splineOuts$obs.preds 
    d1<- d1[d1$FID==plot.which,]
    vals<- d1[,4]
    depths<- d1[,2:3]
    matX<- matrix(NA, nrow=nrow(d1), ncol= 4)
    matY<- matrix(NA, nrow=nrow(d1), ncol= 4)
    for (i in 1: nrow(d1)){
      matX[i,]<- c(vals[i]-vals[i], vals[i], vals[i], vals[i]-vals[i]) 
      matY[i,]<- c(depths[i,1], depths[i,1], depths[i,2], depths[i,2]) }
    raw.plot<-plot(matX[1,], matY[1,],type="n",ylim=c(maxd,0),xlim=c(min(vals), (max(vals)*1.1)),main=paste("soil profile:",d1[1,1], sep=" "), ylab="depth", xlab= label, lty=2, lwd=1.5, xaxs="i", col="black", font.lab=2,cex.lab=1.5)
    
    d2<- as.matrix(splineOuts$harmonised[plot.which,2:length(d)])
    matX1<- matrix(NA, nrow=ncol(d2), ncol= 4)
    matY1<- matrix(NA, nrow=ncol(d2), ncol= 4)
    for (i in 1: ncol(d2)){
      matX1[i,]<- c(d2[i]-d2[i], d2[i], d2[i], d2[i]-d2[i]) 
      matY1[i,]<- c(d[1,i], d[1,i], d[1,i+1], d[1,i+1]) }
    
    #plot the spline averages
    for (i in 1: ncol(d2)){
      polygon (matX1[i,],matY1[i,], lty=1, lwd=1,col="green", border="green") }
    #plot the observations
    for (i in 1: nrow(d1)){
      polygon (matX[i,],matY[i,], lty=1, lwd=2, border="black") }
    #plot the spline
    lines(splineOuts$var.1cm[,plot.which],seq(1,d[length(d)]),lwd=2,col="red" )}}

#END OF SCRIPT
