
### modified cholesky decomposition
modify.cholesky<-function(mat,delta,beta)
{
  nn=nrow(mat)
  C=matrix(0,nrow=nn,ncol=nn)
  L=C
  d=rep(0,nn)
#  ncol=ncol(mat)
  for (j in 1:(nn-1))
  {
    
    diag(C)=diag(mat)
    maxi=which.max(diag(mat))
    if (maxi!=j)
    {
      mat[,c(j,maxi)]=mat[,c(maxi,j)]
      mat[c(j,maxi),]=mat[c(maxi,j),]
    }
    if (j==1&nn<=1)
    {
      C[(j+1):nn,j]=mat[(j+1):nn,j]
    }
    else
    {
      C[(j+1):nn,j]=mat[(j+1):nn,j]-C[(j+1):nn,1:(j-1)]%*%matrix(L[j,1:(j-1)],ncol=1)
    }
    theta=max(abs(C[(j+1):nn,j]))
    d[j]=max(theta,(theta/beta)^2,delta)
    L[(j+1):nn,j]=C[(j+1):nn,j]/d[j]
    diag(C)[(j+1):nn]=diag(C)[(j+1):nn]-L[(j+1):nn,j]*C[(j+1):nn,j]
  }
  d[nn]=max(C[nn,nn],(C[nn,nn]/beta)^2,delta)
  return(list(L=L+diag(1,nn),D=diag(d),C=C))
}

mat=matrix(c(1,1,2,1,1,3,2,3,1),ncol=3,byrow=T)
mat=matrix(c(0.1,1,1,0.1),ncol=2,byrow=T)

(mat_cho=modify.cholesky(mat,delta=0,beta=10000))

mat_cho$L%*%mat_cho$D%*%t(mat_cho$L)


