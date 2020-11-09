mat.sqrt<-function(A){
  ei<-eigen(A);
  d<-ei$values;
  d<-(d+abs(d))/2;
  d2<-sqrt(d);
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors);
  return(ans)
}

normalize <- function(v){
  if(sum(v!=0) > 0)
    return(v/sqrt(sum(v^2)))
  if(sum(v !=0) ==0 )
    return(rep(0,length(v)))
}


force_sp <- function(df1){
  
  a1 <- dim(df1)[2];df2 <- df1;pos <- NULL
  for(i in 1:a1){
    df2[pos,i] <- 0;
    pos <- c(which(df2[,i] !=0))
  }
  return(df2)
}

### Taken from elasticnet/PMA packages by Witten et al.###
##########################################################

soft <- function (a, para){
  b <- sort(abs(a))
  b <- abs(a) - para
  b <- (b + abs(b))/2
  b <- sign(a) * b
  b
}
