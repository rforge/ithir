# Purpose        : Plot raw soil profile data
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 
# Contributions  : 
# Status         : working
# Note           : 


##### Function for plotting raw soil profile data values
plot_soilProfile <- function(data, vals, depths, label= "") {
  matX<- matrix(NA, nrow=nrow(data), ncol= 4)
  matY<- matrix(NA, nrow=nrow(data), ncol= 4)
  for (i in 1: nrow(data)){
    matX[i,]<- c(vals[i]-vals[i], vals[i], vals[i], vals[i]-vals[i]) 
    matY[i,]<- c(depths[i,1], depths[i,1], depths[i,2], depths[i,2]) }
  raw.plot<-plot(matX[1,], matY[1,],type="n",ylim=c(max(depths),0),xlim=c(min(vals), (max(vals)*1.1)),main=paste("soil profile:",data[1,1], sep=" "), ylab="depth", xlab= label, lty=2, lwd=3, xaxs="i", col="black", font.lab=2,cex.lab=1.5)
  for (i in 1: nrow(data)){
    polygon (matX[i,],matY[i,], lty=1, lwd=3, border="black") }}

#END OF SCRIPT
