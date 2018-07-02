# BARCM-Collapsed Gibbs
# simulation

setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")
source("./BARCM.R")

jobID=as.numeric(commandArgs(trailingOnly=T)[1])

folderID = (jobID-1) %/% 100 + 1
folder.char = paste('./Output/sim-0/batch', folderID, sep = '')
dir.create(folder.char, showWarnings = F, recursive = T)

param=expand.grid(seedID=c(1:100),std=c(1,2,5,10,20))
set.seed(param[jobID,1])

############### data read ###############
M = 69 #number of rankers
N = 108 #number of entities
p = 11 #number of covariate
nComp = 3 # number of components
nGroup = 9 # number of gruops
nEach = N/nGroup # number of entities in each group

MCMC = 3000 #number of samples
burn = 1000

Beta0 = array(dim=c(p,nComp))
for(i in 1:nComp){
  Beta0[,i] = rnorm(p)
}

Sigma0 = param[jobID,2]^2 #var of patient error term

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

idCat0 = sample(nComp, M, replace = T, prob = c(0.5,0.3,0.2))
#idCat0 = c(rep(1,40), rep(2,20), rep(3,9))
alpha0 = array(dim=c(N,nComp)) 

#
Y0<-array(dim=c(M,N)) #true value of Y*, generated as following

for(i in 1:nComp){
  alpha0[,i] = rnorm(N)*sqrt(2)
}

for(j in 1:N){
  X[j,] = mvrnorm(1,rep(0,p),CovMat)
  Y0[,j] = rnorm(M, X[j,]%*%Beta0[,idCat0]+alpha0[j,idCat0], sqrt(Sigma0))
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

accuracy_all = rep(0,5)
  
for (cat_num in 3:5){
# hyper-param
lambda=5
s0=1
nComp=cat_num

# define variables
idCat = array(dim=c(M,MCMC)) #category of each ranker
alpha = array(dim=c(N,nComp,MCMC)) 
Beta = array(dim=c(p,nComp,MCMC)) 
pi_comps = array(dim=c(nComp,MCMC))

# introduce expansion variable
expansion = rep(1,MCMC)

# initial value for alpha and beta
Beta[,,1] = 0#solve(t(X)%*%X)%*%t(X)%*%Y[1,,1]
alpha[,,1] = 0
idCat[,1]= sample(nComp, M, replace = T)
expansion[1] = 1
pi_comps[,1] = 1/nComp

# matrix for future use
H = X%*%solve(t(X)%*%X)%*%t(X)

# Gibbs iterations

for(n in 1:(MCMC-1)){
  if(n%%50==0){
    print(n)
    print(expansion[n])
    print(idCat[,n])
    print(pi_comps[,n])
  }
  
  pi_comps[,n+1] = Gibbs_step_pi(pi_comps[,n], idCat[,n])
  
  idCat[,n+1] = Gibbs_step_idCat_collapsed(Y[,,n], idCat[,n], pi_comps[,n+1], H)
  #idCat[,n+1] = Gibbs_step_idCat(Y[,,n], idCat[,n], Beta[,,n], alpha[,,n], pi_comps[,n+1])
  
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], idCat[,n+1], Beta[,,n], alpha[,,n], expansion[n])
  }
  
  temp = Calculate_S_B(Y[,,n+1], idCat[,n+1], alpha[,,n])
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  for(category in 1:dim(Beta)[2]){
    Beta[,category,n+1]=mvrnorm(1,temp$Beta_hat[[category]]/expansion[n+1],temp$V_hat[[category]])
  }
  
  alpha[,,n+1] = Gibbs_step_alpha(Y[,,n+1], idCat[,n+1], Beta[,,n+1], alpha[,,n], lambda)
}


rowMeans(pi_comps[,burn:MCMC])

Y_agg = array(dim=c(N,nComp)) #aggregated Y
for(id in 1:nComp){
  Y_agg[,id] = apply(X%*%Beta[,id,burn:MCMC]+alpha[,id,burn:MCMC],1,mean)
  Y_agg[,id] = Y_agg[,id] - mean(Y_agg[,id])
}

idCat_est = idCat[,MCMC]
for(i in 1:M){
  idCat_est[i] = as.numeric(names(which.max(table(idCat[i,burn:MCMC]))))
}

mu0 = cbind(X%*%Beta0[,1]+alpha0[,1], X%*%Beta0[,2]+alpha0[,2], X%*%Beta0[,3]+alpha0[,3])
cor_mat_k = cor(mu0, Y_agg, method = 'kendall')
cor_mat = cor(mu0, Y_agg)

matching = rep(0,dim(cor_mat)[1])
for(i in 1:dim(cor_mat)[1]){
  matching[i] = which.max(cor_mat[i,])
}
KDist = c(cor_mat_k[1, matching[1]], cor_mat_k[2, matching[2]], cor_mat_k[3, matching[3]])

# accuracy
sum = 0
for(i in 1:M){
  if (matching[idCat0[i]]==idCat_est[i]) sum = sum + 1
}
accuracy = sum/M
accuracy_all[cat_num] = accuracy
}

saveRDS(accuracy_all, paste0(folder.char, "/acc",jobID,".rds"))
