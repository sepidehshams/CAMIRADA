adjacency <- function(PPI){
  nPPI=nrow(PPI)
  c1<-unique(PPI[,1])
  c2<-unique(PPI[,2])
  c_elements<-unique(c(c1,c2))
  m<-length(c_elements)
  A_matrix=matrix(0,m,m)
  for(i in 1:nPPI){
    x=which(PPI[i,1]==c_elements[])
    if(length(x>0))
    {
      y=which(PPI[i,2]==c_elements[])
      if(length(y>0))
      {
        A_matrix[x,y]=1
        A_matrix[y,x]=1
      }
    }
  }
  return(list( A_matrix,c_elements))
}

