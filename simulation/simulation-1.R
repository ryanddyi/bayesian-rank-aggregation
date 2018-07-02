# BARC(W)
# simulation

#setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")
source("./BARCW.R")

jobID=as.numeric(commandArgs(trailingOnly=T)[1])
#jobID=2

folderID = (jobID-1) %/% 500 + 1
folder.char = paste('./Output/sim-1/batch', folderID, sep = '')
dir.create(folder.char, showWarnings = F, recursive = T)

library(MCMCpack)
library(TopKLists)
library(RankAggreg)
library(StatRank)

param=expand.grid(seedID=c(1:500),m=10,n=50,std=c(5,10,20,40))
#hyperparam=expand.grid(lambda=c(5,10),s2=c(0.5,1))
set.seed(param[jobID,1])

Beta0 = c(3,2,1,0.5) #true value
Sigma0 = param[jobID,4]^2 #sd of patient error term

MCMC = 3000

############### data read ###############
M = param[jobID,2] #number of rankers
N = param[jobID,3]  #number of entities
p = length(Beta0) #number of covariate
nComp = 3 # number of components
nGroup = 1 # number of gruops
nEach = N/nGroup # number of entities in each group


X<-array(dim=c(N,p)) #covariate matrix


# Generate X
rho=0.2
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


#========================Naive methods=====================

#BC
order = N+1-Y[,,1]
KDist=cor(-colSums(order),mu0,method="kendall")

#========================MC methods========================

#input1=as.list(data.frame(t(showorder[,1,])))

input1 = list()
for(i in 1:M){
  input1[[i]] = ranked_entities[[i]][[1]]
}

outMC=MC(input1)
MC_agg<-array(dim=c(N,3)) #aggregated Y
for(i in 1:N){
  MC_agg[outMC$MC1.TopK[i],1]=i
  MC_agg[outMC$MC2.TopK[i],2]=i
  MC_agg[outMC$MC3.TopK[i],3]=i
}

#MC1
KDist=cbind(KDist,cor(-MC_agg[,1],mu0,method="kendall"))
#MC2
KDist=cbind(KDist,cor(-MC_agg[,2],mu0,method="kendall"))
#MC3
KDist=cbind(KDist,cor(-MC_agg[,3],mu0,method="kendall"))

#========================BARCW=======================
# hyper-param
lambda=5
s0=1

burn =500

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

# Gibbs iterations

for(n in 1:(MCMC-1)){
  #if(n%%100==0){
  #  print(n)
  #  print(expansion[n])
  #}
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], Beta[,n], alpha[,n], weight[,n], expansion[n])
  }
  
  temp = Calculate_S_B(Y[,,n+1], alpha[,n], weight[,n])
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  Beta[,n+1]=mvrnorm(1,temp$Beta_hat/expansion[n+1],temp$V_hat)

  alpha[,n+1] = Gibbs_step_alpha(Y[,,n+1], Beta[,n+1], alpha[,n], weight[,n], lambda)
  
  weight[,n+1]=1
}


mu = X%*%Beta[,burn:MCMC]+alpha[,burn:MCMC]
for(n in 1:dim(mu)[2]){
  mu[,n] = mu[,n]-mean(mu[,n])
}

Y_agg = apply(mu,1,mean)

KDist=cbind(KDist,cor(Y_agg,mu0,method="kendall"))

####################### PL-base method #######################
input.pl = array(dim=c(M,N))
for(i in 1:M){
  input.pl[i,] = ranked_entities[[i]][[1]]
}

res.pl = Estimation.PL.MLE(input.pl)

KDist=cbind(KDist,cor(res.pl$Mean^-1,mu0,method="kendall"))

write.table(KDist,paste(folder.char, "/KD",jobID,".txt",sep=""), col.names=F,row.names=F)
