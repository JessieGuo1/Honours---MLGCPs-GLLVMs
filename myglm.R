#we want to update beta to
# score equation is X'(y-mu)
# second deriv is -(X'WX)
# so NEwton Raphson says
# beta_new=beta_old + (X'WX)^-1 X'(y-mu)

# ridge penalty:
# replace second deriv with (X'WX+eps*I)

myglm=function(y,X,beta0,offset=0,nIter=50,ridgeParam=0,wts=rep(1,length(y)),eps=1.e-8)
{
  beta=beta_old=beta0
  p=ncol(X)
  if(length(ridgeParam)!=p & length(ridgeParam)!=1)
  {
    message("Ridge parameter has the wrong dimensions")
    break    
  }
  if(length(ridgeParam)==p)
     penVec=diag(ridgeParam)
  if(length(ridgeParam)==1)
    penVec=ridgeParam*diag(p)
  ll_old=1
  for(iIter in 1:nIter)
  {
    eta = X%*%beta + offset
    mu=poisson()$linkinv(eta)
    ll_new=sum(dpois(y,mu,log=TRUE)*wts)
    deltaLL = (ll_new-ll_old)/ll_old
    print(c(beta,ll_new,deltaLL))
    if(abs(deltaLL)<eps){
      message("convergence")
      break
    } 
    Info = crossprod(X*sqrt(wts*mu))+penVec # does elementwise multiplication, column-by-column
    Xmu = crossprod(X,wts*(y-mu))
    beta = beta_old + solve(Info) %*% t(Xmu)
    beta_old=beta
    ll_old=ll_new
  }
}

n=10
X=matrix(1:n)
y=rpois(n,lambda=X)
myglm(y,X,beta0=0)
