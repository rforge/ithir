# Purpose        : predict cluster memberships given input data
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 
# Note:          : This function is really only suitable when dealing with outputs from fuzme
# Note:          : software and specifiaclly when the fuzzy kmeans with extagrades algorithm is used.


fuzzyEx<- function(data,centroid,cv,expon,alfa){
  classNo<-nrow(centroid)+1
  output<-matrix(nrow=nrow(data),ncol=(1+classNo))
  temp<-matrix(ncol=nrow(centroid),nrow=3)
  for(k in 1:nrow(data))
  {
    for(i in 1:nrow(centroid))
    {
      temp[1,i]<-sqrt(mahalanobis(data[k,],centroid[i,],cv))
      temp[2,i]<-temp[1,i]^(-2/(expon-1))
      temp[3,i]<-temp[1,i]^-2
    }
    dSum<-sum(temp[2,])
    dSum2<-sum(temp[3,])
    for(i in 1:nrow(centroid))
    {
      output[k,i]<-temp[1,i]^(-2/(expon-1))/(dSum+(((1-alfa)/alfa)*dSum2)^(-1/(expon-1)))
    }
    output[k,classNo]<-1-(sum(output[k,c(1:nrow(centroid))]))
    output[k,ncol(output)]<-which.max(output[k,c(1:(1+classNo))])
  }
  return(output)
}

#END Script