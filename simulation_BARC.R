# BARC(W)
# simulation

setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")
source("./BARCW.R")

############### data read ###############
M = 10 #number of rankers
N = 50 #number of entities
p = 3 #number of covariate
nComp = 3 # number of components
nGroup = 1 # number of gruops
nEach = N/nGroup # number of entities in each group

MCMC = 1000 #number of samples

Beta0 = cbind(c(3,2,1), c(1,-1,1), c(1,1,-1)) #true value
Sigma0 = 5 #sd of patient error term


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

idCat0 = sample(nComp, M, replace = T, prob = c(1,0,0))
alpha0 = array(dim=c(N,nComp)) 

#
Y0<-array(dim=c(M,N)) #true value of Y*, generated as following

for(i in 1:nComp){
  alpha0[,i] = 2*rnorm(N)
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

#========================BARCW=======================
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
  if(n%%100==0){
    print(n)
    print(expansion[n])
  }
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], Beta[,n], alpha[,n], weight[,n])
  }
  
  temp = Calculate_S(Y[,,n+1], weight[,n], IX, Lambda_inv)
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  eta = mvrnorm(1,temp$eta_hat/expansion[n+1],temp$Sigma_hat)
  alpha[,n+1] = head(eta, N)
  Beta[,n+1] = tail(eta, p)

  weight[,n+1]=1
}

burn = 200

mu = X%*%Beta[,burn:MCMC]+alpha[,burn:MCMC]
for(n in 1:dim(mu)[2]){
  mu[,n] = mu[,n]-mean(mu[,n])
}

Y_agg = apply(mu,1,mean)

rowMeans(Beta[,burn:MCMC])

order(Y_agg)

order(X%*%Beta0[,1]+alpha0[,1])

cor(Y_agg, X%*%Beta0[,1]+alpha0[,1])

acf(Beta[1,burn:MCMC],100, main='')

