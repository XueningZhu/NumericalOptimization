
### simplex algorithm
### A is for coefficient matrix ; B
### keep the relative order permanant

pivot<-function(P,B,L,k)
{
  E=diag(1,nrow(P))
  E[,L]=-P[,k]/P[L,k]
  E[L,L]=1/P[L,k]
  
  P_new=E%*%P
  B_new=E%*%B
#  P_new[,c(l,k)]=P_new[,c(k,l)]
  return(list(P=P_new,B=B_new))
}

get.sigma<-function(P,C,ind)
{
  Z=t(P[,-ind])%*%C[ind,]
  Sigma=C[-ind,]-Z
  return(Sigma)
}

Solution.exist<-function(Sigma,P,ind)
{
  sig_pos_ind=Sigma>0
  if (all(!sig_pos_ind)) return(TRUE)
  if (sum(sig_pos_ind)==1) cond=all(P[,-ind][,sig_pos_ind]<=0)
  else
  {
    cond=apply(P[,-ind][,sig_pos_ind],2,function(v)
      return(all(v<=0)))
    cond=any(cond)
    
  }
  return(!cond)
  return(TRUE)
}

Simplex<-function(A,B,C)
{
  n_base=nrow(B)
  P=cbind(diag(1,n_base),A)
  C=rbind(matrix(c(rep(0,n_base)),ncol=1),C)
  ind=1:n_base    ### ind is the base-ind, initialized with 1:n_base
  flag=T
  n=1
  while(flag&n<=10)
  {
    n=n+1
    Sigma=get.sigma(P,C,ind)
    cat("B:",B,"\n")
    cat("ind:",ind,"\n")
    cat("C:",C[ind],"\n")
    
    cat("Sigma:",Sigma,"\n")
    cat("\n")
    if (all(Sigma<=0)) 
    {
      opt_z=as.vector(t(C[ind,])%*%B)
      opt_x=rep(0,ncol(P))
      opt_x[ind]=as.vector(B)
      return(list(Opt.allX=opt_x,Opt.X=opt_x[-(1:n_base)],Opt.Z=opt_z))
    }
    if(!Solution.exist(Sigma,P,ind))
    {
      cat("No Solution!","\n")
      return(NULL)
    }
    ### decide the in/out vector; update the ind
    in_var=setdiff(1:ncol(P),ind)[which.max(Sigma)[1]]
    
    tmp=B/P[,in_var]
    ind1=which(tmp>0)
    tmp1=tmp[ind1,]
    out_ind=ind1[which.min(tmp1)[1]]
    
    cat("in_var:",in_var,"out_var:",ind[out_ind],"out_ind:",out_ind,"\n")
    #  out_var=ind[out_ind]
    PB=pivot(P,B,L=out_ind,k=in_var)
    ind[out_ind]=in_var
    
    P=PB$P
    B=PB$B
#    cat("P:",P,"\n")
    
  }
  
}

A=matrix(c(1,4,0,2,0,4),ncol=2)
B=matrix(c(8,16,12),ncol=1)
C=matrix(c(2,3),ncol=1)

opt_res=Simplex(A,B,C)


















