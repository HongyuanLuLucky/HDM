##Generate Data
## For this data I will include beta_0
library(MASS)
library(glmnet)
###########create the training data############
set.seed(99)
ntr=80
p=120
sig=diag(p)
for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
x=mvrnorm(ntr,rep(0,p),sig)
beta=c(0,0,6,5,-5,rep(0,p-5))
p1=1/(1+exp(-x%*%beta))
y=rep(0,ntr)
y[p1>0.5]=1
# getY=function(prob)
# {
#   return(rbinom(1,1,prob))
# }
# y=sapply(p1,getY)


#############create the test data################
nte=200
p=120
sig=diag(p)
for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
xe=mvrnorm(nte,rep(0,p),sig)
beta=c(0,0,6,5,-5,rep(0,p-5))
p1e=1/(1+exp(-xe%*%beta))
# ye=sapply(p1e,getY)
ye=rep(0,nte)
ye[p1e>0.5]=1
####################

beta_ini=rep(0.01,p)

#y=rep(0,ntr)
#y[p1>0.5]=1

#Soft Thresholding
S=function(a,b)
{if(b<abs(a) & a>0){return(a-b)}
  else if(b<abs(a) & a<0){return(a+b)}
  else{return(0)}
}


#The function to run lasso and return the beta
LGRlasso=function(lambda,trx,try,beta_ini)
{p=length(beta_ini)
 ntr=dim(trx)[1]
 maxiter=500
 beta=beta_ini
 
 w=rep(0.25,ntr)
 #ct=0
 c=1
while(c<1000){
   c=c+1
   #print(ct)
   #ct=ct+1
   phat=1/(1+exp(-trx%*%beta))
   phat[phat==0]=0+10^-7
   phat[phat==1]=1-10^-7
   #print("a")
   #print(phat)
   z=trx%*%beta+(try-phat)/((phat)*(1-phat))
   #print("b")
   #print(z)
   beta_new=beta
   for(k in 1:p)
   {
     fk=1/ntr*t((z-trx[,-k]%*%beta[-k]))%*%(w*trx[,k])
     beta_new[k]=S(fk,lambda)/(sum(w*(trx[,k])^2)/ntr)
     #print(c(k,":",beta_new[k]))
     
     
   }
   #print(c("mse:",mean((beta-beta_new)^2)))
   #print("fun")
   #print("are you sure1")
   if(mean((beta-beta_new)^2)<10^-9){break}
   #print("are you sure2")
   #print(c("ak",beta_new))
   #print(c("ak1",beta[1:3]))
   beta=beta_new
   #print(c("ak2",beta[1:3]))

 }  
 return(beta)
}



#####The function to run elastic net and return beta
LGRelsnet=function(lambda,alpha,trx,try,beta_ini)
{p=length(beta_ini)
ntr=dim(trx)[1]
maxiter=500
beta=beta_ini

w=rep(0.25,ntr)
#ct=0
c=1
while(c<1000){
  c=c+1
  #print(ct)
  #ct=ct+1
  phat=1/(1+exp(-trx%*%beta))
  phat[phat==0]=0+10^-7
  phat[phat==1]=1-10^-7
  #print("a")
  #print(phat)
  z=trx%*%beta+(try-phat)/((phat)*(1-phat))
  #print("b")
  #print(z)
  beta_new=beta
  for(k in 1:p)
  {
    fk=1/ntr*t((z-trx[,-k]%*%beta[-k]))%*%(w*trx[,k])
    beta_new[k]=S(fk,alpha*lambda)/(sum(w*(trx[,k])^2)/ntr+lambda*(1-alpha))
    #print(c(k,":",beta_new[k]))
    
    
  }
  #print(c("mse:",mean((beta-beta_new)^2)))
  #print("fun")
  #print("are you sure1")
  if(mean((beta-beta_new)^2)<10^-7){break}
  #print("are you sure2")
  #print(c("ak",beta_new))
  #print(c("ak1",beta[1:3]))
  beta=beta_new
  #print(c("ak2",beta[1:3]))
  
}  
return(beta)
}


#The function to run mcp
LGReMCP=function(lambda,alpha,trx,try,beta_ini)
{p=length(beta_ini)
ntr=dim(trx)[1]
maxiter=500
beta=beta_ini

#w=rep(0.25,ntr)
#x=scale(x)/0.5
c=1
while(c<1000){
  c=c+1
  #print(ct)
  #ct=ct+1
  phat=1/(1+exp(-trx%*%beta))
  phat[phat==0]=0+10^-7
  phat[phat==1]=1-10^-7
  
  w=phat*(1-phat)
  #print("a")
  #print(phat)
  z=trx%*%beta+(try-phat)/((phat)*(1-phat))
  #print("b")
  #print(z)
  beta_new=beta
  
  for(k in 1:p)
  { 
    fk=1/ntr*t((z-trx[,-k]%*%beta[-k]))%*%(w*trx[,k])
    if(abs(fk)<mean(alpha*lambda*w*trx[,k]^2))
    {beta_new[k]=S(fk,lambda)/(sum(w*(trx[,k])^2)/ntr-1/alpha)}else{
     beta_new[k]=fk/(mean(w*trx[,k]^2))  
    }
    #print(c(k,":",beta_new[k]))
    
    
  }
  #print(c("mse:",mean((beta-beta_new)^2)))
  #print("fun")
  #print("are you sure1")
  if(mean((beta-beta_new)^2)<10^-7){break}
  #print("are you sure2")
  #print(c("ak",beta_new))
  #print(c("ak1",beta[1:3]))
  beta=beta_new
  #print(c("ak2",beta[1:3]))
  
}  
return(beta)
}

###Your welcome to play the follwoing functions
r1=LGRlasso(0.15,x,y,beta_ini)
r2=LGRelsnet(0.15,1,x,y,beta_ini)
r3=LGRelsnet(0.146,0.2,x,y,beta_ini)
r4=LGReMCP(0.75,2,x,y,beta_ini)


#the function to run lasso for a sequence of lambda
runlasso=function(lambdas,x,y,xe,ye,beta_ini,beta)
{ tp=c()
  fp=c()
  mc=c() # the number of misclassification. 
  dv=c()
  beta_hat=Null
  std_tp=sum(beta!=0)
  for(lambda in lambdas){
    beta_hat=LGRlasso(lambda,x,y,beta_ini)
    
    phat=1/(1+exp(-xe%*%beta_hat))#
    yhat=rep(0,length(ye))#
    yhat[phat>0.5]=1#
    
    dv[ye==1]=sqrt(2*(log(1+exp(xe%*%beta_hat))-xe%*%beta_hat))[ye==1]
    
    dv[ye==0]=-2*sqrt(2*log(1+exp(xe%*%beta_hat)))[ye==0]
    d=sum(dv^2)
    
    
    mc=c(mc,sum(yhat!=ye))#
    tp=c(tp,sum(beta!=0 & beta_hat!=0))
    fp=c(fp,sum(beta==0 & beta_hat!=0))
    
    # a nice print out to show a table of lambda,dv,misclassification, tp,fp
    print(c(lambda,d,sum(yhat!=ye),sum(beta!=0 & beta_hat!=0),sum(beta==0 & beta_hat!=0)))
  }
  # print("-------")
  # print(tp)
  # print(std_tp)
  # print(fp)
  # print(min(fp))
  # print(length(tp))
  # print(length(fp))
  opt_lambda_index=which(mc==min(mc))#
  opt_lambda=lambdas[opt_lambda_index]
  opt_dv=dv[opt_lambda_index]
  opt_tp=tp[opt_lambda_index]
  opt_fp=fp[opt_lambda_index]
  opt_mc=mc[opt_lambda_index]
  
  opt_dv_index=which(opt_dv==min(opt_dv))[1]
  opt_lambda=opt_lambda[opt_dv_index]
  opt_tp=opt_tp[opt_dv_index]
  opt_fp=opt_fp[opt_dv_index]
  opt_mc=opt_mc[opt_dv_index]
  return(c(opt_lambda,opt_mc,opt_tp,opt_fp))
}

##Your are welcome to play around here 
lambdas=seq(0.9,0.1,-0.01)
re1=runlasso(lambdas,x,y,xe,ye,beta_ini,beta)
r1=LGRlasso(0.11,x,y,beta_ini)


#the function to ran elstic net for a sequence of lambda
runelslasso=function(lambdas,x,y,xe,ye,alpha,beta_ini,beta)
{ tp=c()
  fp=c()
  mc=c() #
  dv=c()
beta_hat=Null
std_tp=sum(beta!=0)
for(lambda in lambdas){
  beta_hat=LGRelsnet(lambda,alpha,x,y,beta_ini)
  
  phat=1/(1+exp(-xe%*%beta_hat))#
  yhat=rep(0,length(ye))#
  yhat[phat>0.5]=1#
  
  dv[ye==1]=sqrt(2*(log(1+exp(xe%*%beta_hat))-xe%*%beta_hat))[ye==1]
  
  dv[ye==0]=-2*sqrt(2*log(1+exp(xe%*%beta_hat)))[ye==0]
  d=sum(dv^2)
  
  mc=c(mc,sum(yhat!=ye))#
  tp=c(tp,sum(beta!=0 & beta_hat!=0))
  fp=c(fp,sum(beta==0 & beta_hat!=0))
  
  # a nice print out to show a table of lambda,dv,misclassification, tp,fp
  print(c(lambda,d,sum(yhat!=ye),sum(beta!=0 & beta_hat!=0),sum(beta==0 & beta_hat!=0)))
}
# print("-------")
# print(tp)
# print(std_tp)
# print(fp)
# print(min(fp))
# print(length(tp))
# print(length(fp))
opt_lambda_index=which(mc==min(mc))#
opt_lambda=lambdas[opt_lambda_index]
opt_dv=dv[opt_lambda_index]
opt_tp=tp[opt_lambda_index]
opt_fp=fp[opt_lambda_index]
opt_mc=mc[opt_lambda_index]

opt_dv_index=which(opt_dv==min(opt_dv))[1]
opt_lambda=opt_lambda[opt_dv_index]
opt_tp=opt_tp[opt_dv_index]
opt_fp=opt_fp[opt_dv_index]
opt_mc=opt_mc[opt_dv_index]
return(c(opt_lambda,opt_mc,opt_tp,opt_fp))

}

#Your are welcome to play the following functions here 
lambdas=seq(0.9,0.1,-0.01)
re2=runelslasso(lambdas,x,y,xe,ye,0.6,beta_ini,beta)
r2=LGRelsnet(0.14,0.5,x,y,beta_ini)

#the function to run MCP for a sequence of lambda
runMCP=function(lambdas,alpha,x,y,xe,ye,beta_ini,beta)
{ tp=c()
fp=c()
mc=c()
dv=c()
beta_hat=Null
std_tp=sum(beta!=0)
for(lambda in lambdas){
  beta_hat=LGReMCP(lambda,alpha,x,y,beta_ini)
  
  phat=1/(1+exp(-xe%*%beta_hat))#
  yhat=rep(0,length(ye))#
  yhat[phat>0.5]=1#
  
  dv[ye==1]=sqrt(2*(log(1+exp(xe%*%beta_hat))-xe%*%beta_hat))[ye==1]
  
  dv[ye==0]=-2*sqrt(2*log(1+exp(xe%*%beta_hat)))[ye==0]
  d=sum(dv^2)
 
  mc=c(mc,sum(yhat!=ye))#
  
  tp=c(tp,sum(beta!=0 & beta_hat!=0))
  fp=c(fp,sum(beta==0 & beta_hat!=0))
  # a nice print out to show a table of lambda,dv,misclassification, tp,fp
  print(c(lambda,d,sum(yhat!=ye),sum(beta!=0 & beta_hat!=0),sum(beta==0 & beta_hat!=0)))
  if(sum(yhat!=ye)>min(mc)){break}
}
# print("-------")
# print(tp)
# print(std_tp)
# print(fp)
# print(min(fp))
# print(length(tp))
# print(length(fp))
opt_lambda_index=which(mc==min(mc))[1]#
opt_lambda=lambdas[opt_lambda_index]
opt_dv=dv[opt_lambda_index]
opt_tp=tp[opt_lambda_index]
opt_fp=fp[opt_lambda_index]
opt_mc=mc[opt_lambda_index]

opt_dv_index=which(opt_dv==min(opt_dv))[1]
opt_lambda=opt_lambda[opt_dv_index]
opt_tp=opt_tp[opt_dv_index]
opt_fp=opt_fp[opt_dv_index]
opt_mc=opt_mc[opt_dv_index]
return(c(opt_lambda,opt_mc,opt_tp,opt_fp))
}

#Your are welcome to play the following functions here
lambdas=seq(0.9,0.1,-0.01)
re3=runMCP(lambdas,3,x,y,xe,ye,beta_ini,beta)
r3=LGReMCP(0.42,3,x,y,beta_ini)

#########Run the lambda vs deviance and lambda vs misclassification
#####lasso Graph
# lambdas=seq(0.9,0.045,-0.001)
# tp=c()
# fp=c()
# mc=c() # the number of misclassification. 
# dv=c()
# dvv=rep(0,length(ye))
# beta_hat=Null
# std_tp=sum(beta!=0)
# for(lambda in lambdas){
#   beta_hat=LGRlasso(lambda,x,y,beta_ini)
#   
#   phat=1/(1+exp(-xe%*%beta_hat))#
#   yhat=rep(0,length(ye))#
#   yhat[phat>0.5]=1#
#   
#   dvv[ye==1]=sqrt(2*(log(1+exp(xe%*%beta_hat))-xe%*%beta_hat))[ye==1]
#   
#   dvv[ye==0]=-2*sqrt(2*log(1+exp(xe%*%beta_hat)))[ye==0]
#   dv=c(dv,sum(dvv^2))
#   
# 
#   
#   mc=c(mc,sum(yhat!=ye))#
#   tp=c(tp,sum(beta!=0 & beta_hat!=0))
#   fp=c(fp,sum(beta==0 & beta_hat!=0))
#   print(c(lambda,sum(dvv^2),sum(yhat!=ye),sum(beta!=0 & beta_hat!=0),sum(beta==0 & beta_hat!=0)))
# }
# 
# par(mfrow=c(1,2))
# plot(log(lambdas)[500:length(lambdas)],dv[500:length(lambdas)],col="red",main="Logistic Regression and Lasso",ylab="deviance",xlab=expression(log(lambda)))
# plot(log(lambdas)[500:length(lambdas)],mc[500:length(lambdas)],col="red",main="Logistic Regression and Lasoo",ylab="misclassfication",xlab=expression(log(lambda)),ylim=c(0,40))

#####elastic Net Graph testing size 200
set.seed(39)
ntr=80
p=120
sig=diag(p)
for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
x=mvrnorm(ntr,rep(0,p),sig)
beta=c(0,0,6,5,-5,rep(0,p-5))
p1=1/(1+exp(-x%*%beta))
y=rep(0,ntr)
y[p1>0.5]=1
nte=200
p=120
sig=diag(p)
for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
xe=mvrnorm(nte,rep(0,p),sig)
beta=c(0,0,6,5,-5,rep(0,p-5))
p1e=1/(1+exp(-xe%*%beta))
ye=rep(0,nte)
ye[p1e>0.5]=1
lambdas=seq(0.9,0.01,-0.001)
tp=c()
fp=c()
mc=c() #
dv=c()

beta_hat=Null
std_tp=sum(beta!=0)
for(lambda in lambdas){
  beta_hat=LGRelsnet(lambda,0.9,x,y,beta_ini)
  
  phat=1/(1+exp(-xe%*%beta_hat))#
  yhat=rep(0,length(ye))#
  yhat[phat>0.5]=1#

  dvv=rep(0,length(ye))
  dvv[ye==1]=sqrt(2*(log(1+exp(xe%*%beta_hat))-xe%*%beta_hat))[ye==1]
  
  dvv[ye==0]=-2*sqrt(2*log(1+exp(xe%*%beta_hat)))[ye==0]
  dv=c(dv,sum(dvv^2))
  
  mc=c(mc,sum(yhat!=ye))#
  tp=c(tp,sum(beta!=0 & beta_hat!=0))
  fp=c(fp,sum(beta==0 & beta_hat!=0))
  print(c(lambda,sum(dvv^2),sum(yhat!=ye),sum(beta!=0 & beta_hat!=0),sum(beta==0 & beta_hat!=0)))
}
par(mfrow=c(1,2))
plot(log(lambdas)[500:length(lambdas)],dv[500:length(lambdas)],col="red",main="Logistic Regression and elastic net",ylab="deviance",xlab=expression(log(lambda)))
plot(log(lambdas)[500:length(lambdas)],mc[500:length(lambdas)],col="red",main="Logistic Regression and elastic net",ylab="misclassication",xlab=expression(log(lambda)),ylim=c(0,40))
#####elastic Net Graph testing size 30
set.seed(39)
ntr=80
p=120
sig=diag(p)
for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
x=mvrnorm(ntr,rep(0,p),sig)
beta=c(0,0,6,5,-5,rep(0,p-5))
p1=1/(1+exp(-x%*%beta))
y=rep(0,ntr)
y[p1>0.5]=1
nte=30
p=120
sig=diag(p)
for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
xe=mvrnorm(nte,rep(0,p),sig)
beta=c(0,0,6,5,-5,rep(0,p-5))
p1e=1/(1+exp(-xe%*%beta))
ye=rep(0,nte)
ye[p1e>0.5]=1
lambdas=seq(0.9,0.01,-0.001)
tp=c()
fp=c()
mc=c() #
dv=c()

beta_hat=Null
std_tp=sum(beta!=0)
for(lambda in lambdas){
  beta_hat=LGRelsnet(lambda,0.9,x,y,beta_ini)
  
  phat=1/(1+exp(-xe%*%beta_hat))#
  yhat=rep(0,length(ye))#
  yhat[phat>0.5]=1#
  
  dvv=rep(0,length(ye))
  dvv[ye==1]=sqrt(2*(log(1+exp(xe%*%beta_hat))-xe%*%beta_hat))[ye==1]
  
  dvv[ye==0]=-2*sqrt(2*log(1+exp(xe%*%beta_hat)))[ye==0]
  dv=c(dv,sum(dvv^2))
  
  mc=c(mc,sum(yhat!=ye))#
  tp=c(tp,sum(beta!=0 & beta_hat!=0))
  fp=c(fp,sum(beta==0 & beta_hat!=0))
  print(c(lambda,sum(dvv^2),sum(yhat!=ye),sum(beta!=0 & beta_hat!=0),sum(beta==0 & beta_hat!=0)))
}
par(mfrow=c(1,2))
plot(log(lambdas)[500:length(lambdas)],dv[500:length(lambdas)],col="red",main="Logistic Regression and elastic net",ylab="deviance",xlab=expression(log(lambda)))
plot(log(lambdas)[500:length(lambdas)],mc[500:length(lambdas)],col="red",main="Logistic Regression and elastic net",ylab="misclassication",xlab=expression(log(lambda)))

##########MCP Graph
# lambdas=seq(0.9,0.01,-0.001)
# tp=c()
# fp=c()
# mc=c()
# dv=c()
# 
# beta_hat=Null
# std_tp=sum(beta!=0)
# for(lambda in lambdas){
#   beta_hat=LGReMCP(lambda,3,x,y,beta_ini)
#   
#   phat=1/(1+exp(-xe%*%beta_hat))#
#   yhat=rep(0,length(ye))#
#   yhat[phat>0.5]=1#
#   dvv=rep(0,length(ye))
#   dvv[ye==1]=sqrt(2*(log(1+exp(xe%*%beta_hat))-xe%*%beta_hat))[ye==1]
#   
#   dvv[ye==0]=-2*sqrt(2*log(1+exp(xe%*%beta_hat)))[ye==0]
#   dv=c(dv,sum(dvv^2))
#   
#   
#   mc=c(mc,sum(yhat!=ye))#
#   tp=c(tp,sum(beta!=0 & beta_hat!=0))
#   fp=c(fp,sum(beta==0 & beta_hat!=0))
#   print(c(lambda,sum(dvv^2),sum(yhat!=ye),sum(beta!=0 & beta_hat!=0),sum(beta==0 & beta_hat!=0)))
# }
# par(mfrow=c(1,2))
# plot(log(lambdas),dv,col="red",main="Logistic Regression and MCP",ylab="deviance",xlab=expression(log(lambda)))
# plot(log(lambdas),mc,col="red",main="Logistic Regression and MCP",ylab="misclassication",xlab=expression(log(lambda)))

########Run(50 Replicates) for testing size of 200
beta_ini=rep(0.01,p)
lasso.lambda=c()
lasso.mc=c()
lasso.tp=c()
lasso.fp=c()

elsnet.lambda=c()
elsnet.mc=c()
elsnet.tp=c()
elsnet.fp=c()


mcp.lambda=c()
mcp.mc=c()
mcp.tp=c()
mcp.fp=c()
for (rep in 1:50)
{#########generate the data
  ntr=80
  p=120
  sig=diag(p)
  for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
  x=mvrnorm(ntr,rep(0,p),sig)
  beta=c(0,0,6,5,-5,rep(0,p-5))
  p1=1/(1+exp(-x%*%beta))
  y=rep(0,ntr)
  y[p1>0.5]=1
  
  nte=200
  p=120
  sig=diag(p)
  for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
  xe=mvrnorm(nte,rep(0,p),sig)
  beta=c(0,0,6,5,-5,rep(0,p-5))
  p1e=1/(1+exp(-xe%*%beta))
  # ye=sapply(p1e,getY)
  ye=rep(0,nte)
  ye[p1e>0.5]=1
  
  print(rep)
  lambdas=seq(0.9,0.1,-0.01)
  re1=runlasso(lambdas,x,y,xe,ye,beta_ini,beta)
  re2=runelslasso(lambdas,x,y,xe,ye,0.5,beta_ini,beta)
  re3=runMCP(lambdas,3,x,y,xe,ye,beta_ini,beta)
  
  
  
  lasso.lambda=c(lasso.lambda,re1[1])
  lasso.mc=c(lasso.mc,re1[2])
  lasso.tp=c(lasso.tp,re1[3])
  lasso.fp=c(lasso.fp,re1[4])
  
  elsnet.lambda=c(elsnet.lambda,re2[1])
  elsnet.mc=c(elsnet.mc,re2[2])
  elsnet.tp=c(elsnet.tp,re2[3])
  elsnet.fp=c(elsnet.fp,re2[4])
  
  mcp.lambda=c(mcp.lambda,re3[1])
  mcp.mc=c(mcp.mc,re3[2])
  mcp.tp=c(mcp.tp,re3[3])
  mcp.fp=c(mcp.fp,re3[4])
  
}

########Run(50 Replicates) for testing size of 30
beta_ini=rep(0.01,p)
lasso.lambda=c()
lasso.mc=c()
lasso.tp=c()
lasso.fp=c()

elsnet.lambda=c()
elsnet.mc=c()
elsnet.tp=c()
elsnet.fp=c()


mcp.lambda=c()
mcp.mc=c()
mcp.tp=c()
mcp.fp=c()
for (rep in 1:50)
{#########generate the data
  ntr=80
  p=120
  sig=diag(p)
  for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
  x=mvrnorm(ntr,rep(0,p),sig)
  beta=c(0,0,6,5,-5,rep(0,p-5))
  p1=1/(1+exp(-x%*%beta))
  y=rep(0,ntr)
  y[p1>0.5]=1
  
  nte=30
  p=120
  sig=diag(p)
  for(i in 1:p){for(j in 1:p){sig[i,j]=0.5^abs(i-j)}} 
  xe=mvrnorm(nte,rep(0,p),sig)
  beta=c(0,0,6,5,-5,rep(0,p-5))
  p1e=1/(1+exp(-xe%*%beta))
  # ye=sapply(p1e,getY)
  ye=rep(0,nte)
  ye[p1e>0.5]=1
  
  print(rep)
  lambdas=seq(0.9,0.1,-0.01)
  re1=runlasso(lambdas,x,y,xe,ye,beta_ini,beta)
  re2=runelslasso(lambdas,x,y,xe,ye,0.5,beta_ini,beta)
  re3=runMCP(lambdas,3,x,y,xe,ye,beta_ini,beta)
  
  
  
  lasso.lambda=c(lasso.lambda,re1[1])
  lasso.mc=c(lasso.mc,re1[2])
  lasso.tp=c(lasso.tp,re1[3])
  lasso.fp=c(lasso.fp,re1[4])
  
  elsnet.lambda=c(elsnet.lambda,re2[1])
  elsnet.mc=c(elsnet.mc,re2[2])
  elsnet.tp=c(elsnet.tp,re2[3])
  elsnet.fp=c(elsnet.fp,re2[4])
  
  mcp.lambda=c(mcp.lambda,re3[1])
  mcp.mc=c(mcp.mc,re3[2])
  mcp.tp=c(mcp.tp,re3[3])
  mcp.fp=c(mcp.fp,re3[4])
  
} 
