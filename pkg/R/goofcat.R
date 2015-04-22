# Purpose        : Goodness of fit statistics for categorical varibles
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 
# Contributions  : 
# Status         : working
# Note           : 
#

goofcat<- function(observed = NULL, predicted = NULL, conf.mat, imp=FALSE){
  
  if (imp==TRUE){
    if(class(conf.mat)!="matrix"){
      stop("Entered data is NOT a matrix")}
    if(nrow(conf.mat)!= ncol(conf.mat)) {
      stop("Entered data is NOT a confusion matrix")}
    
    else {OA<- ceiling(sum(diag(conf.mat))/sum(colSums(conf.mat)) * 100)
      PA<- ceiling(diag(conf.mat)/colSums(conf.mat) * 100)
      UA<- ceiling(diag(conf.mat)/rowSums(conf.mat) * 100)
      
      PE_mat <- matrix(NA, ncol = 1, nrow = length(rowSums(conf.mat)))
      for (i in 1:length(rowSums(conf.mat))) {
      PE_mat[i, 1] <- (rowSums(conf.mat)[i]/sum(colSums(conf.mat))) * (colSums(conf.mat)[i]/sum(colSums(conf.mat)))}
      KS <- (sum(diag(conf.mat))/sum(colSums(conf.mat)) - sum(PE_mat))/(1 - sum(PE_mat))}}
  
  if (imp==FALSE) {obsMat<- table(observed,observed)
     df<- data.frame(observed, predicted)
     names(df)<- c("observed", "predicted")
     #make a confusion matrix
     cfuM<- function(df,obsMat){
       c.Mat<- as.matrix(obsMat)
       snames1<- c(colnames(c.Mat))
       for (i in 1:nrow(c.Mat)){
         for (j in 1:nrow(c.Mat)){
           c.Mat[j,i]<- nrow(subset(df, df$observed ==snames1[i]  & df$predicted ==snames1[j]))}}
       fmat<- matrix(NA, nrow=nrow(c.Mat), ncol=ncol(c.Mat))
       rownames(fmat)<- rownames(c.Mat)
       colnames(fmat)<- colnames(c.Mat)
       for (i in 1:nrow(c.Mat)){
         fmat[i,]<- c(c.Mat[,i])}
       return(fmat)}
     conf.mat<- cfuM(df, obsMat)
     
     OA<- ceiling(sum(diag(conf.mat))/sum(colSums(conf.mat)) * 100)
     PA<- ceiling(diag(conf.mat)/colSums(conf.mat) * 100)
     UA<- ceiling(diag(conf.mat)/rowSums(conf.mat) * 100)
     
     PE_mat <- matrix(NA, ncol = 1, nrow = length(rowSums(conf.mat)))
     for (i in 1:length(rowSums(conf.mat))) {
       PE_mat[i, 1] <- (rowSums(conf.mat)[i]/sum(colSums(conf.mat))) * (colSums(conf.mat)[i]/sum(colSums(conf.mat)))}
     KS <- (sum(diag(conf.mat))/sum(colSums(conf.mat)) - sum(PE_mat))/(1 - sum(PE_mat))}
  retval<- list(conf.mat,OA, PA, UA, KS)
  names(retval)<- c("confusion_matrix", "overall_accuracy", "producers_accuracy", "users_accuracy", "kappa")
  return(retval)}
  
  # end script

  
