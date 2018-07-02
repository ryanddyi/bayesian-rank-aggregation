#setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")

#nfl data

covariates<-as.matrix(read.table("data-read/QB_cov.txt"))
EBranking<-read.table("data-read/QB_EBF.txt",sep="\t")
ranking<-as.matrix(EBranking[,-c(1,9)])

#beta independent+equal known variance.

library("MCMCpack")

MCMC=5000 #number of samples
M=dim(ranking)[2] #number of rankers
N=dim(ranking)[1] #number of entities in each group
G=1 #number of groups
K=dim(covariates)[2] #number of covariate

order<-array(dim=c(M,G,N)) #largest one rank first
showorder<-array(dim=c(M,G,N)) #indicator of the i-th item
X<-array(dim=c(G,N,K)) #covariate matrix
Y<-array(dim=c(M,G,N,MCMC)) #latent value of Y*
#gamma<-array(dim=c(G,N,MCMC)) #the part that doctors agree
alpha<-array(dim=c(G,N,MCMC)) #the part that doctors agree
Beta<-array(dim=c(K,MCMC)) 
weight<-array(dim=c(M,MCMC))

lambda=5
s0=1

#read covariates in different groups

for(i in 1:dim(covariates)[2]){
  covariates[,i] = (covariates[,i] - mean(covariates[,i]))/sd(covariates[,i])
}

X[1,,]=covariates

Y[,1,,1]=25-t(ranking)

#find out the order and set initial value of Y*
for(i in 1:M){
	for(g in 1:G){
		for(j in 1:N){
			count=0
			for(k in 1:N){
				if(Y[i,g,k,1]>=Y[i,g,j,1]) count=count+1
			}
		order[i,g,j]=count
		}
	}
}

#find out indicator of different order
for(i in 1:M){
	for(j in 1:G){
		for(k in 1:N){
			for(m in 1:N){
				if(order[i,j,m]==k){
					showorder[i,j,k]=m
				}
			}
		}
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

for(n in 1:(MCMC-1)){
  if(n%%50==0){
    print(n)
    print(expansion[n])
    print(idCat[,n])
  }
  
  idCat[,n+1] = Gibbs_step_idCat_DP(Y[,,n], idCat[,n], H, 1)
  
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

mu_cat = array(dim = c(N,max(idCat_est)))
for (cat in unique(idCat_est)){
  if (sum(idCat_est==cat)==1){
    mu_cat[,cat] = Y_agg[,idCat_est==cat]
  } else {
    mu_cat[,cat] = rowMeans(Y_agg[,idCat_est==cat])
  }
}


pi_est = table(idCat[,burn:MCMC])/sum(table((idCat[,burn:MCMC])))
#pi_order = order(pi_est, decreasing=T)
pi_order = as.numeric(names(table(idCat_est))[order(table(idCat_est),decreasing=T)])
#pi_order = order(weight_cat,decreasing=T)

cat_count = c()
re_order = c()
for(i in 1:length(pi_order)){
  re_order = c(re_order, which(idCat_est==pi_order[i]))
  cat_count = c(cat_count,length(re_order))
}

Y_est = Y[,,MCMC]
for(i in 1:M){
  Y_est[i,] = apply(Y[i,,burn:MCMC],1,mean)
}

require(reshape2)
require(ggplot2)


cor_mat = cor(t(Y_est[re_order,]), method='kendall')
#cor_mat = cor(Y_agg[,re_order], method='kendall')

melted_covmat <- melt(cor_mat)
p <- ggplot(data = melted_covmat, aes(x=Var1, y=Var2, fill=value)) + theme_bw() + geom_tile() + 
  scale_fill_gradient(low="white", high="red") + xlab('') + ylab('') + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  #theme(legend.position="none") + 
  geom_vline(xintercept = 0.5+cat_count[1:(M-1)], size = 0.1) + 
  geom_hline(yintercept = 0.5+cat_count[1:(M-1)], size = 0.1)
p

pdf('corrplot_z.pdf', height = 5, width = 6)
print(p)
dev.off()


cor_mat = cor(Y_agg[,re_order], method='kendall')

melted_covmat <- melt(cor_mat)
p <- ggplot(data = melted_covmat, aes(x=Var1, y=Var2, fill=value)) + theme_bw() + geom_tile() + 
  scale_fill_gradient(low="white", high="red") + xlab('') + ylab('') + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0)) + 
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  #theme(legend.position="none") + 
  geom_vline(xintercept = 0.5+cat_count[1:(M-1)], size = 0.1) + 
  geom_hline(yintercept = 0.5+cat_count[1:(M-1)], size = 0.1)
p

pdf('corrplot_mu.pdf', height = 5, width = 6)
print(p)
dev.off()
