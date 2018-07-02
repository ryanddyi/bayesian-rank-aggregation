# BARCM-Collapsed Gibbs
# simulation

#setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")
source("./BARCM.R")

jobID=as.numeric(commandArgs(trailingOnly=T)[1])
#jobID=1

folderID = (jobID-1) %/% 100 + 1
folder.char = paste('./Output/sim-7/batch', folderID, sep = '')
dir.create(folder.char, showWarnings = F, recursive = T)

param=expand.grid(seedID=c(1:100))
set.seed(param[jobID,1])

############### data read ###############
M = 69 #number of rankers
N = 108 #number of entities
p = 11 #number of covariate
nComp = 1 # number of components
nGroup = 9 # number of gruops
nEach = N/nGroup # number of entities in each group

MCMC = 1500 #number of samples

Beta0 = array(dim=c(p,nComp))
for(i in 1:nComp){
  Beta0[,i] = rnorm(p)
}

Sigma0 = 1 #sd of patient error term

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

idCat0 = rep(1,M)
alpha0 = array(dim=c(N,nComp)) 

#
Y0<-array(dim=c(M,N)) #true value of Y*, generated as following

for(i in 1:nComp){
  alpha0[,i] = rnorm(N)*2
}

for(j in 1:N){
  X[j,] = mvrnorm(1,rep(0,p),CovMat)
  Y0[,j] = rnorm(M, X[j,]%*%Beta0[,idCat0]+alpha0[j,idCat0], Sigma0)
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

#========================BARCM=======================
source("./BARCM.R")
# hyper-param
s_alpha = 1
s_beta = 100

# define variables
idCat = array(dim=c(M,MCMC)) #category of each ranker
alpha = array(0,dim=c(N,M,MCMC)) 
Beta = array(0,dim=c(p,M,MCMC)) 

# introduce expansion variable
expansion = rep(1,MCMC)

# initial value
idCat[,1]= 1
expansion[1] = 1

# V matrix in paper (cbind(I,X))
IX = cbind(diag(N),X)
Lambda_inv = diag(c(rep(1/s_alpha^2,N), rep(1/s_beta^2, p)))

# matrix for future use
# H[[i]] is V (L^-1 + i*V'V)^-1 V
H = list()
for(i in 1:M){
  H[[i]] = IX%*%solve(i*t(IX)%*%IX+Lambda_inv)%*%t(IX)
}

# Gibbs iterations

accuracy_all = rep(0,4)
nCat_all = rep(0,4)
for(log_dp_count in 1:4){
for(n in 1:(MCMC-1)){
  if(n%%50==0){
    print(n)
    print(expansion[n])
    print(idCat[,n])
  }
  
  idCat[,n+1] = Gibbs_step_idCat_DP(Y[,,n], idCat[,n], H, M^(log_dp_count/2-1.5))
  
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y_DP(i, Y[,,n], Beta[,i,n], alpha[,i,n])
  }
  
  temp = Calculate_S_DP(Y[,,n+1], idCat[,n+1], IX, Lambda_inv)
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  for(category in unique(idCat[,n+1])){
    eta = mvrnorm(1,temp$eta_hat[[category]]/expansion[n+1],temp$Sigma_hat[[category]])
    ranker_cat = which(idCat[,n+1]==category)
    alpha[,ranker_cat,n+1] = head(eta, N)
    Beta[,ranker_cat,n+1] = tail(eta, p)
  }
}

burn = 500

mu = array(dim=c(N,M,MCMC-burn+1)) #aggregated Y
for(id in 1:M){
  mu[,id,] = X%*%Beta[,id,burn:MCMC]+alpha[,id,burn:MCMC]
  for(i in 1:dim(mu)[3]){
    mu[,id,i] = mu[,id,i] - mean(mu[,id,i])
  }
}

Y_agg = array(dim = c(N,M))
for(id in 1:M){
  Y_agg[,id] = apply(mu[,id,],1,mean)
}

idCat_est = idCat[,MCMC]
for(i in 1:M){
  idCat_est[i] = as.numeric(names(which.max(table(idCat[i,burn:MCMC]))))
}
nCat = rep(0,MCMC)
for(i in 1:MCMC){
  nCat[i] = length(unique(idCat[,i]))
}

rand_count = 0
for(i in 1:(M-1)){
  for(j in (i+1):M){
    if ((idCat0[i]==idCat0[j])==(idCat_est[i]==idCat_est[j])){
      rand_count = rand_count+1
    }
  }
}

accuracy = rand_count/(M*(M-1)/2)

accuracy_all[log_dp_count] = accuracy
nCat_all[log_dp_count] = mean(nCat[burn:MCMC])
}


saveRDS(accuracy_all, paste0(folder.char, "/acc",jobID,".rds"))
saveRDS(nCat_all, paste0(folder.char, "/nCat",jobID,".rds"))


#========================BARCW=======================
source("./BARCW.R")
# hyper-param
s_alpha = 1
s_beta = 100

# define variables
alpha = array(dim=c(N,MCMC)) 
Beta = array(dim=c(p,MCMC)) 
weight = array(dim=c(M,MCMC))

# introduce expansion variable
expansion = rep(1,MCMC)

# V matrix in paper (cbind(I,X))
IX = cbind(diag(N),X)
Lambda_inv = diag(c(rep(1/s_alpha^2,N), rep(1/s_beta^2, p)))

# initial value for alpha and beta
Beta[,1] = 0
alpha[,1] = 0
expansion[1] = 1
weight[,1] = 1

# Gibbs iterations

for(n in 1:(MCMC-1)){
  if(n%%50==0){
    print(n)
    print(expansion[n])
    print(weight[,n])
  }
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], Beta[,n], alpha[,n], weight[,n], expansion[n])
  }
  
  temp = Calculate_S(Y[,,n+1], weight[,n], IX, Lambda_inv)
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  eta = mvrnorm(1,temp$eta_hat/expansion[n+1],temp$Sigma_hat)
  alpha[,n+1] = head(eta, N)
  Beta[,n+1] = tail(eta, p)
  
  weight[,n+1] = Gibbs_step_weight(Y[,,n+1], Beta[,n+1], alpha[,n+1], weight[,n], X)
}

weight_est = rowMeans(weight[,burn:MCMC])

saveRDS(weight_est, paste0(folder.char, "/weight",jobID,".rds"))

