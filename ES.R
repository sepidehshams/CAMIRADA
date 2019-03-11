ES <- function(L,X){
  n=length(L)
  n1=length(X)
  es<- rep(0, n1) 
  m=0
  es[]=which(L%in%X)
  es=sort(es, decreasing = FALSE)
  x=es-(1:n1)
  e=((1:n1)*((n-n1)/n1))-((x)*(n1/(n-n1)))
  z<-setdiff(L, X)
  y<-which(L%in%z)
  e<-c(e,(-1*y*(n1/(n-n1))))
  return(max(e))
}



