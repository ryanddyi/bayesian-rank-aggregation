# BARC
# simulation
# partial list 

setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")
source("./BARCW.R")

jobID=as.numeric(commandArgs(trailingOnly=T)[1])
#jobID=201

folderID = (jobID-1) %/% 100 + 1
folder.char = paste('./Output/sim-6/batch', folderID, sep = '')
dir.create(folder.char, showWarnings = F, recursive = T)

param=expand.grid(seedID=c(1:100),m=10,g=c(1,2,4,8,10,16),std=5)
set.seed(param[jobID,1])

Beta0 = c(3,2,1) #true value
Sigma0 = param[jobID,4]^2 #sd of patient error term

MCMC = 3000

############### data read ###############
M = param[jobID,2] #number of rankers
N = 80 #number of entities
p = length(Beta0) #number of covariate
nComp = 3 # number of components
nGroup = param$g[jobID] # number of gruops
nEach = N/nGroup # number of entities in each group


X<-array(dim=c(N,p)) #covariate matrix

# Generate X
rho=0.5
CovMat=diag(p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    CovMat[i,j]=rho^(abs(i-j))
    CovMat[j,i]=rho^(abs(i-j))
  }
}

alpha0 = rep(0,N)
mu0 = rep(0,N)

#
Y0<-array(dim=c(M,N)) #true value of Y*, generated as following


for(j in 1:N){
  X[j,] = mvrnorm(1,rep(0,p),CovMat)
  alpha0[j] = sum(X[j,]^2)
  mu0[j] = X[j,]%*%Beta0 + alpha0[j]
  Y0[,j] = rnorm(M, mu0[j], sqrt(Sigma0))
}

# subgroups of entities for all ranker (could be different for different rankers)
subgroups = list()
for(i in 1:nGroup){
  subgroups[[i]] = nEach*(i-1)+1:nEach
}

ranked_entities = list() # decreasingly ranked entities
for(i in 1:M){
  ranked_entities[[i]] = list()
}

Y = array(dim=c(M,N,MCMC)) #latent value of Y*

# find out ranked entities
for (j in 1:nGroup){
  for(i in 1:M){
    rank_subgroup = nEach+1-rank(Y0[i,subgroups[[j]]])
    Y[i,subgroups[[j]],1] = rank(Y0[i,subgroups[[j]]]) # set initial value of Y*
    ranked_entities[[i]][[j]] = order(rank_subgroup) + (j-1)*nEach
  }
}

KDist = c()

#####################

# hyper-param
lambda=10
s0=1

burn = 500

# define variables
alpha = array(dim=c(N,MCMC)) 
Beta = array(dim=c(p,MCMC)) 
weight = array(dim=c(M,MCMC))

# introduce expansion variable
expansion = rep(1,MCMC)

# initial value for alpha and beta
Beta[,1] = 0
alpha[,1] = 0
expansion[1] = 1
weight[,1] = 1

#########BARC(1)

for(n in 1:(MCMC-1)){
  if(n%%100==0){
    print(n)
    print(expansion[n])
  }
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], Beta[,n], alpha[,n], weight[,n], expansion[n])
  }
  
  temp = Calculate_S_B(Y[,,n+1], alpha[,n], weight[,n], X)
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  Beta[,n+1] = mvrnorm(1,temp$Beta_hat/expansion[n+1],temp$V_hat)
  
  alpha[,n+1] = Gibbs_step_alpha(Y[,,n+1], Beta[,n+1], alpha[,n], weight[,n], lambda)
  
  weight[,n+1]=1
}


mu = X%*%Beta[,burn:MCMC]+alpha[,burn:MCMC]
for(n in 1:dim(mu)[2]){
  mu[,n] = mu[,n]-mean(mu[,n])
}

Y_agg = apply(mu,1,mean)

KDist=cbind(KDist,cor(Y_agg,mu0,method="kendall"))

#########BARC(0.5)
s0=0.5

for(n in 1:(MCMC-1)){
  if(n%%100==0){
    print(n)
    print(expansion[n])
  }
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], Beta[,n], alpha[,n], weight[,n], expansion[n])
  }
  
  temp = Calculate_S_B(Y[,,n+1], alpha[,n], weight[,n], X)
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  Beta[,n+1] = mvrnorm(1,temp$Beta_hat/expansion[n+1],temp$V_hat)
  
  alpha[,n+1] = Gibbs_step_alpha(Y[,,n+1], Beta[,n+1], alpha[,n], weight[,n], lambda)
  
  weight[,n+1]=1
}


mu = X%*%Beta[,burn:MCMC]+alpha[,burn:MCMC]
for(n in 1:dim(mu)[2]){
  mu[,n] = mu[,n]-mean(mu[,n])
}

Y_agg = apply(mu,1,mean)

KDist=cbind(KDist,cor(Y_agg,mu0,method="kendall"))



######## BAR #######
s0=1
for(n in 1:(MCMC-1)){
  if(n%%100==0){
    print(n)
    print(expansion[n])
  }
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], Beta[,n], alpha[,n], weight[,n], expansion[n])
  }
  
  temp = Calculate_S_B(Y[,,n+1], alpha[,n], weight[,n], X)
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  Beta[,n+1]=0#mvrnorm(1,temp$Beta_hat/expansion[n+1],temp$V_hat)
  
  alpha[,n+1] = Gibbs_step_alpha(Y[,,n+1], Beta[,n+1], alpha[,n], weight[,n], lambda)
  
  weight[,n+1]=1
}


mu = X%*%Beta[,burn:MCMC]+alpha[,burn:MCMC]
for(n in 1:dim(mu)[2]){
  mu[,n] = mu[,n]-mean(mu[,n])
}

Y_agg = apply(mu,1,mean)

KDist=cbind(KDist,cor(Y_agg,mu0,method="kendall"))

write.table(KDist,paste(folder.char, "/KD",jobID,".txt",sep=""), col.names=F,row.names=F)
