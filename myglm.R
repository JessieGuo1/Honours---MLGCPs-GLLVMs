#we want to update beta to
# score equation is X'(y-mu)
# second deriv is -(X'WX)
# so NEwton Raphson says
# beta_new=beta_old + (X'WX)^-1 X'(y-mu)

# ridge penalty:
# replace second deriv with (X'WX+eps*I)

library(Matrix)
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
    penVec=ridgeParam
  if(length(ridgeParam)==1)
    penVec=rep(ridgeParam,p)
  penMat = Diagonal(p,penVec)
  ll_old=1
  linkinv=poisson()$linkinv
  for(iIter in 1:nIter)
  {
    eta = crossprod(t(X),beta) + offset
    mu=as.vector(linkinv(eta))
    penTerm = crossprod(beta,sqrt(penVec))/2
    ll_new=as.numeric(sum(dpois(y,mu,log=TRUE)*wts)-penTerm)
    deltaLL = (ll_new-ll_old)/ll_old
    print(c(as.numeric(beta),as.numeric(ll_new),as.numeric(deltaLL)))
    if(abs(deltaLL)<eps){
      message("convergence")
      break
    } 
    W=wts*mu
    Info = crossprod(X,X*W)+penMat # does elementwise multiplication, column-by-column
    score = crossprod(X,wts*(y-mu)) - penVec*beta
    cf=Cholesky(Info)
    beta = beta_old + solve(cf,score)
    beta_old=beta
    ll_old=ll_new
  }
}

n=100
p = 0.2 #sparsity of X
X=cbind(1,matrix(1:n)/n*rbinom(n,1,p))
Xs=Matrix(X)
y=rpois(n,lambda=exp(X[,2]))
myglm(y,X,beta0=c(0,0))
glm(y~0+X,family=poisson()) #check
system.time(myglm(y,X,beta0=c(0,0))) #this is faster for v big n (>10^5) and sparse X
system.time(glm(y~0+X,family=poisson()))
