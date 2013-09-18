### »Æ½ð·Ö¸î·¨

f<-function(x)
{
  return((x[1]+x[2]^2)^2)
}
df<-function(x)
{
  return(c(2*(x[1]+x[2]^2),4*(x[1]+x[2]^2)*x[2]))
}

alpha.interval<-function(x,y,d,t,h)  ## t>1 h>0
{
  alpha0=0
  alpha1=alpha0+t*h
  while( (f(x+alpha1*d)-f(x+alpha0*d))<0)
  {
    alpha0=alpha1
    alpha1=alpha1+t*h
  }
  cat("alpha.interval: alpha0",alpha0,"alpha1:",alpha1,"\n")
  return(c(alpha0,alpha1))
}

golden.opt<-function(x,y,d,alp_int,eps)
{
  a=alp_int[1]
  b=alp_int[2]
  
  alpha_l=a+(1-0.618)*(b-a)
  alpha_r=a+0.618*(b-a)
  
  y_l=f(x+alpha_l*d)
  y_r=f(x+alpha_r*d)
  
  while(b-a>eps)
  {
#    cat("golden.opt: alpha_l:",alpha_l,"alpha_r:",alpha_r,"\n")
    if(y_l<=y_r)
    {
      b=alpha_r
      alpha_r=alpha_l
      alpha_l=a+(1-0.618)*(b-a)
      y_r=y_l
      y_l=f(x+alpha_l*d)
    }
    else
    {
      a=alpha_l
      alpha_l=alpha_r
      alpha_r=a+0.618*(b-a)
      y_l=y_r
      y_r=f(x+alpha_r*d)
    }
  }
  return((a+b)/2)
}

opt.gradient<-function(f,df,x0=c(3,3),eps_all=10^(-3))
{
  y0=f(x0)
  d0=-df(x0)
  h=sqrt(sum(x0^2))/1000
  alp_int=alpha.interval(x0,y0,t=1.5,h=h,d0)   ### get alpha interval
  eps=alp_int[1]*10^(-4)
  alpha=golden.opt(x0,y0,d0,alp_int,eps)
  x1=x0+alpha*d0
  y1=f(x1)
  while(abs((y1-y0)/y0)>eps_all)
  {
    cat("opt.gradient: y:",y1,"\n")
    x0=x1
    y0=f(x0)
    d0=-df(x0)
    h=sqrt(sum(x0^2))/1000
    alp_int=alpha.interval(x0,y0,t=1.5,h=h,d0)   ### get alpha interval
    eps=10^(-4)
    alpha=golden.opt(x0,y0,d0,alp_int,eps)
    x1=x0+alpha*d0
    y1=f(x1)
  }
  return(list(x=x1,y=y1))
}


min_f=opt.gradient(f,df,x0=c(3,3),eps_all=10^(-4))










