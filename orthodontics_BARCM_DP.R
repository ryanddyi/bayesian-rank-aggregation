# orthodontics data - rank aggregation mixture model

setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")
source("./BARCM.R")
covariates = as.matrix(read.csv("data-read/covariates.csv"))

##### combine covariates
#covariates0 = covariates
#covariates1 = array(dim = c(dim(covariates0)[1],5))
#covariates1[,1] = rowSums(covariates0[,1:6])
#covariates1[,2] = rowSums(covariates0[,7:8])
#covariates1[,3:5] = rowSums(covariates0[,9:11])
#colnames(covariates1) = c('dmm','bom', 'ojm','obm','clm')
#covariates = covariates1
##########################

MCMC = 50000 #number of samples
burn = 1000

############### data read ###############
M = 69 #number of rankers
N = 108 #number of entities
p=dim(covariates)[2] #number of covariate
nGroup = 9 # number of gruops
nEach = 12 # number of entities in each group

X = array(dim=c(N,p)) #covariate matrix
Y = array(dim=c(M,N,MCMC)) #latent value of Y*

# read covariates
X = covariates[1:N,]
for(i in 1:dim(X)[2]){
  X[,i] = (X[,i] - mean(X[,i]))/sd(X[,i])
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

# find out ranked entities
for (j in 1:nGroup){
  rank_subgroup = t(as.matrix(read.table(paste0("data-read/g",j,".txt"),sep="\t",header=T)[,-1]))
  Y[,subgroups[[j]],1] = nEach+1-rank_subgroup # set initial value of Y*
  for(i in 1:M){
    ranked_entities[[i]][[j]] = order(rank_subgroup[i,]) + (j-1)*nEach
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


##################

Y_agg_all = rowMeans(Y_agg)
Y_agg_other = rowMeans(Y_agg[,re_order[35:69]])

entityname=c(paste0("A",1:12),paste0("B",1:12),paste0("C",1:12),paste0("D",1:12), paste0("E",1:12),
             paste0("F",1:12),paste0("G",1:12),paste0("H",1:12),paste0("I",1:12))
#result=cbind(Y_agg_w, Y_agg_all, mu_cat[,pi_order[1:3]])
result=cbind(Y_agg_all, mu_cat[,pi_order[1:3]], Y_agg_other)
#(1-cor(Y_agg_w, Y_agg_all, method = 'kendall'))/2 # kendall tau distance between BARCM and BARCW.

top_five = array(dim=c(5,5))
for(i in 1:dim(top_five)[2]){
  top_five[,i] = entityname[order(result[,i],decreasing=T)[1:5]]
}
top_five

bot_five = array(dim=c(5,5))
for(i in 1:dim(bot_five)[2]){
  bot_five[,i] = entityname[order(result[,i])[1:5]]
}
bot_five

xtable(bot_five)

########################

# clustering plots
library("cluster")
library("factoextra")
library("magrittr")

pdf('distplot_mu.pdf', height = 8, width = 9)
res.dist <- get_dist(t(Y_agg), stand = TRUE, method = "kendall")/2
#res.dist <- get_dist(Y_est, stand = TRUE, method = "kendall")
#fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_dist(res.dist, gradient = list(low = "#FC4E07", high = "white"))
dev.off()

# diagnositc plots
library('coda')
library('mcmcplots')

for(j in 1:3){
  chain = as.mcmc(mu[1,j,])
  pdf(paste0('orth_densityplot',j,'.pdf'), height = 5, width = 5)
  densplot(chain)
  dev.off()
  pdf(paste0('orth_acfplot',j,'.pdf'), height = 5, width = 5)
  acf(chain, main="", 2000)
  dev.off()
  pdf(paste0('orth_traceplot',j,'.pdf'), height = 5, width = 5)
  traceplot(chain)
  dev.off()
}


#effectiveSize(mu[1,3,1:10000])

#eff_size = array(dim = dim(mu[,,1]))
#for (i in 1:108){
#  for(j in 1:69){
#    eff_size[i,j] = effectiveSize(mu[i,j,1:10000])
#  }
#}

