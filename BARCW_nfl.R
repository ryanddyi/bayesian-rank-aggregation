setwd('/Users/ydd/Documents/BARC/code')
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

#========================BARCW=======================
#combine X into a gn*k vector
X_comb=X[1,,]
if(G>1){
  for(i in 2:G){
    X_comb=rbind(X_comb,X[i,,])
  }
}

#set initial value for beta
Y_comb=Y[1,1,,1]
if(G>1){
  for(j in 2:G){
    Y_comb=c(Y_comb,Y[1,j,,1])
  }
}

Beta[,1]=solve(t(X_comb)%*%X_comb)%*%t(X_comb)%*%Y_comb
alpha[,,1]=0
expansion<-rep(1,MCMC)
weight[,1]=1

#==========hyperparam seting ==========

  ## all possible weight 
  allweight=c(0.5,1,2)
  WeightProb<-function(weight,sumsq,num){
    prob=rep(0,length(allweight))
    logprob=(num/2)*log(weight)-weight*sumsq/2
    if((sort(logprob)[3]-sort(logprob)[2])>100){
      prob[order(logprob)]=c(0,0,1)
    } else if(sort(logprob)[2]-sort(logprob[1])>100){
      prob[order(logprob)[1]]=0
      prob[order(logprob)[-1]]=exp(sort(logprob)[-1]-sort(logprob)[2])/sum(exp(sort(logprob)[-1]-sort(logprob)[2]))
    } else {
      prob[order(logprob)]=exp(sort(logprob)-sort(logprob)[2])/sum(exp(sort(logprob)-sort(logprob)[2]))
    }
    prob
  }
  
  for(n in 1:(MCMC-1)){
    if(n%%200==0){
      print(n)
      print(expansion[n])
      print(weight[,n])
    }
    for(i in 1:M){
      for(j in 1:G){
        Y[i,j,,(n+1)]=Y[i,j,,n]
        Y[i,j,showorder[i,j,1],n+1]=expansion[n]*alpha[j,showorder[i,j,1],n]+betarightnorm(X[j,showorder[i,j,1], ], expansion[n]*Beta[ ,n], expansion[n]/weight[i,n], Y[i,j,showorder[i,j,2],n+1]-expansion[n]*alpha[j,showorder[i,j,1],n])
        for(m in 2:(N-1)){
          Y[i,j,showorder[i,j,m],n+1]=expansion[n]*alpha[j,showorder[i,j,m],n]+betamidnorm(X[j,showorder[i,j,m], ], expansion[n]*Beta[ ,n], expansion[n]/weight[i,n], Y[i,j,showorder[i,j,m+1],n+1]-expansion[n]*alpha[j,showorder[i,j,m],n], Y[i,j,showorder[i,j,m-1],n+1]-expansion[n]*alpha[j,showorder[i,j,m],n])
        }
        Y[i,j,showorder[i,j,N],n+1]=expansion[n]*alpha[j,showorder[i,j,N],n]+betaleftnorm(X[j,showorder[i,j,N], ], expansion[n]*Beta[ ,n], expansion[n]/weight[i,n], Y[i,j,showorder[i,j,(N-1)],n+1]-expansion[n]*alpha[j,showorder[i,j,N],n])
      }
    }
    Y[,,,n+1]=Y[,,,n+1]/expansion[n]
    
    A=0
    for(j in 1:G){
      for(i in 1:M){
        A=A+t(Y[i,j,,n+1]-alpha[j,,n])%*%(Y[i,j,,n+1]-alpha[j,,n])*weight[i,n]
      }
    }
    B=rep(0,K)
    for(j in 1:G){
      for(i in 1:M){
        B=B+t(X[j,,])%*%(Y[i,j,,n+1]-alpha[j,,n])*weight[i,n]
      }
    }
    
    C=array(0,dim=c(K,K))
    for(j in 1:G){
      C=C+t(X[j,,])%*%X[j,,]
    }
    C=C*sum(weight[,n])
    H=solve(C)
    expansion[n+1]=sqrt((A-t(B)%*%H%*%B)/rchisq(1,G*N*M))
    
    Beta[,n+1]=mvrnorm(1,H%*%B/expansion[n+1],H)
    for(j in 1:G){
      for(k in 1:N){
        alpha[j,k,n+1]=midnorm(sum((Y[,j,k,n+1]-X[j,k,]%*%Beta[,n+1])*weight[,n])/(sum(weight[,n])+1/s0),sqrt(1/(sum(weight[,n])+1/s0)),-lambda,lambda)
      }
    }
    
    for(i in 1:M){
      temp=0
      for(j in 1:G){
        for(k in 1:N){
          temp=temp+(Y[i,j,k,n+1]-X[j,k,]%*%Beta[,n+1]-alpha[j,k,n+1])^2
        }
      }
      #weight[i,n+1]=midchisq(G*N+1,temp+1,0.4,2)
      weight[i,n+1]=sample(allweight,size=1,prob=WeightProb(allweight,temp,G*N))
      #print(WeightProb(allweight,temp,G*N))
    }
  }
  
  Y_agg<-array(dim=c(G,N)) #aggregated Y
  for(i in 1:G){
    for(j in 1:N){
      Y_agg[i,j]=X[i,j,]%*%apply(Beta[,500:MCMC],1,mean)+mean(alpha[i,j,500:MCMC])
    }
  }

Y_u<-array(dim=c(G,N,(MCMC-500+1)))
for(i in 1:G){
  for(j in 1:N){
    Y_u[i,j,]=X[i,j,]%*%Beta[,500:MCMC]+alpha[i,j,500:MCMC]
  }
}

rank_u=Y_u

for(k in 1:(MCMC-500+1)){
  rank_u[1,,k]=N+1-rank(Y_u[1,,k])
}

apply(Y_u[1,,],1,function(x){quantile(x,probs=0.025)})
apply(Y_u[1,,],1,function(x){quantile(x,probs=0.975)})

apply(rank_u[1,,],1,max)
apply(rank_u[1,,],1,function(x){quantile(x,probs=0.025)})
apply(rank_u[1,,],1,function(x){quantile(x,probs=0.975)})

rank_agg=array(dim=c(24,2))
#rank_agg[,1]=as.vector(EBranking[,1])
rank_agg[,1]=Y_agg[1,]
rank_agg[,2]=N+1-rank(Y_agg)

#===========draw table/plot of rank with mean and CI.===============
rank_table=array(dim=c(24,4))
rank_table[,1]=as.vector(EBranking[,1])
rank_table[,2]=rank_agg[,2]
rank_table[,3]=apply(rank_u[1,,],1,function(x){quantile(x,probs=0.025)})
rank_table[,4]=apply(rank_u[1,,],1,function(x){quantile(x,probs=0.975)})

rank_table=rank_table[order(as.numeric(rank_table[,2])),]
rank_df<-data.frame(
  qb = factor(rank_table[,1],rank_table[order(as.numeric(rank_table[,2]),decreasing=T),1]),
  group= factor(rank_table[,2],rank_table[,2]),
  left = rank_table[,3],
  right = rank_table[,4]
)
setwd('/Users/ydd/Documents/BARC/tex/')
pdf('rank-CI.pdf', height = 6, width = 9)
p <- ggplot(rank_df, aes(group,qb))
p<-p + geom_point()+geom_errorbarh(aes(xmax = right, xmin = left))+
  labs(x = 'Rank', y = 'Quarterback') +
  theme_bw()
print(p)
dev.off()

barcw=t(weight[,500:MCMC])
boxplot(barcw)

#========================BARC=======================
MCMC=3000

for(n in 1:(MCMC-1)){
  if(n%%200==0){
    print(n)
    print(expansion[n])
    print(weight[,n])
  }
  for(i in 1:M){
    for(j in 1:G){
      Y[i,j,,(n+1)]=Y[i,j,,n]
      Y[i,j,showorder[i,j,1],n+1]=expansion[n]*alpha[j,showorder[i,j,1],n]+betarightnorm(X[j,showorder[i,j,1], ], expansion[n]*Beta[ ,n], expansion[n]/weight[i,n], Y[i,j,showorder[i,j,2],n+1]-expansion[n]*alpha[j,showorder[i,j,1],n])
      for(m in 2:(N-1)){
        Y[i,j,showorder[i,j,m],n+1]=expansion[n]*alpha[j,showorder[i,j,m],n]+betamidnorm(X[j,showorder[i,j,m], ], expansion[n]*Beta[ ,n], expansion[n]/weight[i,n], Y[i,j,showorder[i,j,m+1],n+1]-expansion[n]*alpha[j,showorder[i,j,m],n], Y[i,j,showorder[i,j,m-1],n+1]-expansion[n]*alpha[j,showorder[i,j,m],n])
      }
      Y[i,j,showorder[i,j,N],n+1]=expansion[n]*alpha[j,showorder[i,j,N],n]+betaleftnorm(X[j,showorder[i,j,N], ], expansion[n]*Beta[ ,n], expansion[n]/weight[i,n], Y[i,j,showorder[i,j,(N-1)],n+1]-expansion[n]*alpha[j,showorder[i,j,N],n])
    }
  }
  Y[,,,n+1]=Y[,,,n+1]/expansion[n]
  
  A=0
  for(j in 1:G){
    for(i in 1:M){
      A=A+t(Y[i,j,,n+1]-alpha[j,,n])%*%(Y[i,j,,n+1]-alpha[j,,n])*weight[i,n]
    }
  }
  B=rep(0,K)
  for(j in 1:G){
    for(i in 1:M){
      B=B+t(X[j,,])%*%(Y[i,j,,n+1]-alpha[j,,n])*weight[i,n]
    }
  }
  
  C=array(0,dim=c(K,K))
  for(j in 1:G){
    C=C+t(X[j,,])%*%X[j,,]
  }
  C=C*sum(weight[,n])
  H=solve(C)
  expansion[n+1]=sqrt((A-t(B)%*%H%*%B)/rchisq(1,G*N*M))
  
  Beta[,n+1]=mvrnorm(1,H%*%B/expansion[n+1],H)
  for(j in 1:G){
    for(k in 1:N){
      alpha[j,k,n+1]=midnorm(sum((Y[,j,k,n+1]-X[j,k,]%*%Beta[,n+1])*weight[,n])/(sum(weight[,n])+1/s0),sqrt(1/(sum(weight[,n])+1/s0)),-lambda,lambda)
    }
  }
  
  for(i in 1:M){
    weight[i,n+1]=1
  }
}

Y_agg<-array(dim=c(G,N)) #aggregated Y
for(i in 1:G){
  for(j in 1:N){
    Y_agg[i,j]=X[i,j,]%*%apply(Beta[,500:MCMC],1,mean)+mean(alpha[i,j,500:MCMC])
  }
}

rank_agg=cbind(rank_agg,Y_agg[1,])
rank_agg=cbind(rank_agg,N+1-rank(Y_agg))
#rank_agg=cbind(rank_agg,order(Y_agg[1,],decreasing=T))

Y_u<-array(dim=c(G,N,(MCMC-500+1)))
for(i in 1:G){
  for(j in 1:N){
    Y_u[i,j,]=X[i,j,]%*%Beta[,500:MCMC]+alpha[i,j,500:MCMC]
  }
}

#----------trace plot-------------
burn=1000
Y_u<-array(dim=c(G,N,(MCMC-burn+1)))
for(i in 1:G){
  for(j in 1:N){
    Y_u[i,j,]=X[i,j,]%*%Beta[,burn:MCMC]+alpha[i,j,burn:MCMC]
  }
}
sample_all=Y_u[1,,]
sample_demeaned=sample_all
for(k in 1:24){
  sample_demeaned[k,]=sample_all[k,]-apply(sample_all,2,mean)
}
x=sample_demeaned[1,]
effectiveSize(x)
setwd('/Users/ydd/Documents/BARC/aoas-template/')
pdf('traceplot1.pdf', height = 5, width = 5)
traceplot(as.mcmc(x), ylab='mu_1')
dev.off()
pdf('densityplot1.pdf', height = 5, width = 5)
plot(density(x),ylab = 'density',main='')
dev.off()
pdf('acfplot1.pdf', height = 5, width = 5)
acf(x,100,main='')
dev.off()

#---------------------------------
rank_u=Y_u

for(k in 1:(MCMC-500+1)){
  rank_u[1,,k]=N+1-rank(Y_u[1,,k])
}

apply(Y_u[1,,],1,function(x){quantile(x,probs=0.025)})
apply(Y_u[1,,],1,function(x){quantile(x,probs=0.975)})

apply(rank_u[1,,],1,max)
apply(rank_u[1,,],1,function(x){quantile(x,probs=0.025)})
apply(rank_u[1,,],1,function(x){quantile(x,probs=0.975)})

rank_agg=array(dim=c(24,2))
#rank_agg[,1]=as.vector(EBranking[,1])
rank_agg[,1]=Y_agg[1,]
rank_agg[,2]=N+1-rank(Y_agg)

#===========draw table/plot of rank with mean and CI.===============
rank_table=array(dim=c(24,4))
rank_table[,1]=as.vector(EBranking[,1])
rank_table[,2]=rank_agg[,2]
rank_table[,3]=apply(rank_u[1,,],1,function(x){quantile(x,probs=0.025)})
rank_table[,4]=apply(rank_u[1,,],1,function(x){quantile(x,probs=0.975)})

rank_table=rank_table[order(as.numeric(rank_table[,2])),]
rank_df<-data.frame(
  qb = factor(rank_table[,1],rank_table[order(as.numeric(rank_table[,2]),decreasing=T),1]),
  group= factor(rank_table[,2],rank_table[,2]),
  left = rank_table[,3],
  right = rank_table[,4]
)
setwd('/Users/ydd/Documents/BARC/tex/')
pdf('rank-CI-BARC2.pdf', height = 6, width = 9)
p <- ggplot(rank_df, aes(group,qb))
p<-p + geom_point()+geom_errorbarh(aes(xmax = right, xmin = left))+
  labs(x = 'Rank', y = 'Quarterback') +
  theme_bw()
print(p)
dev.off()

barcw=t(weight[,500:MCMC])
boxplot(barcw)


#=======AriM========

rank_agg=cbind(rank_agg,colSums(order[,1,])/M)
rank_agg=cbind(rank_agg,rank(colSums(order[,1,])/M))
#rank_agg=cbind(rank_agg,order(colSums(order[,1,])/M))

#======MC=====

library(TopKLists)
library(RankAggreg)

input1=as.list(data.frame(t(showorder[,1,])))
outMC=MC(input1)
MC_agg<-array(dim=c(N,3)) #aggregated Y
for(i in 1:N){
  MC_agg[outMC$MC1.TopK[i],1]=i
  MC_agg[outMC$MC2.TopK[i],2]=i
  MC_agg[outMC$MC3.TopK[i],3]=i
}

rank_agg=cbind(rank_agg,outMC$MC1.Prob[MC_agg[,1]])
rank_agg=cbind(rank_agg,MC_agg[,1])
rank_agg=cbind(rank_agg,outMC$MC2.Prob[MC_agg[,2]])
rank_agg=cbind(rank_agg,MC_agg[,2])
rank_agg=cbind(rank_agg,outMC$MC3.Prob[MC_agg[,3]])
rank_agg=cbind(rank_agg,MC_agg[,3])

#agg=RankAggreg(showorder[,1,], N, method="CE",distance="Spearman", verbose = F)$top.list

rank_agg=cbind(as.vector(EBranking[,1]),rank_agg)
result=rank_agg[order(as.numeric(rank_agg[,2]),decreasing = T),]
for(i in 1:6){
  result[,2*i]=round(as.numeric(result[,2*i]),3)
}
xtable(result[,-3])

#=========accuracy=========
accuracy=as.matrix(read.csv(file="/Users/ydd/Documents/BARC/QBdata/accuracy.csv"))
plot(accuracy[3,-1])

ranker.validate=rbind(apply(barcw,2,mean),accuracy[,-1])
ranker.sort=ranker.validate[,order(as.numeric(ranker.validate[1,]),decreasing=T)]

plot(ranker.sort[3,],ranker.sort[1,],ylab="Posterior mean of weight w",xlab="Accuracy in Season 2014")
weight<-as.numeric(ranker.sort[1,])
accuracy<-as.numeric(ranker.sort[3,])
lm.w<-lm(weight~accuracy)
plot(accuracy,weight,ylab="Posterior mean of weight w",xlab="Accuracy in Season 2014")
abline(lm.w,col='blue')
         
reg<-lm(ranker.sort[1,]~ranker.sort[3,])
abline(reg)
plot(ranker.sort[4,],ranker.sort[1,],ylab="Posterior mean of weight w",xlab="Accuracy of week 13")

weight<-as.numeric(ranker.validate[1,])
accuracy<-as.numeric(ranker.validate[4,])
lm.w<-lm(accuracy~weight)
anova(lm.w)


#as.numeric(ranker.validate[4,])%*%as.numeric(ranker.validate[1,])/sum(as.numeric(ranker.validate[1,]))
#mean(as.numeric(ranker.validate[4,]))

plot(ranker.sort[1,])
par(new=T)
plot(ranker.sort[3,],col="red")
par(new=T)
plot(ranker.sort[4,],col="green")

library(ggplot2)
library(RColorBrewer)
library(reshape)
setwd('/Users/ydd/Documents/BARC/tex/')
pdf('qb-weight.pdf', height = 4, width = 8)
pdmelt <- melt(barcw)
p<-ggplot(pdmelt, aes(x=X2, y=value,group = X2))+
  ggtitle("BARCW result - QB ranking") +
  geom_boxplot(fill="darkseagreen4")+
  stat_summary(fun.y=mean, colour="red", geom="point", 
               shape=3, size=3,show_guide = FALSE)+
  #geom_jitter(alpha=0.5,position = position_jitter(width = .2))
  labs(x = 'Ranker', y = 'Weight') +
  theme_bw()
print(p)
dev.off()

barplot(Y_agg[1,]-apply(alpha[1,,500:MCMC],1,mean))
barplot(apply(alpha[1,,500:MCMC],1,mean))

#plot coefficient

library("coefplot")
my.ModelCI = data.frame(Value = apply(Beta[,500:MCMC],1,mean), 
                        Coefficient = c("G","Pct","Att", "Avg", "Yds","TD","Int","RAtt","RAvg","RYds","R1st"),
                        HighInner=apply(Beta[,500:MCMC],1,f<-function(x){quantile(x,0.975)}),
                        LowInner=apply(Beta[,500:MCMC],1,f<-function(x){quantile(x,0.025)}),
                        HighOuter=apply(Beta[,500:MCMC],1,f<-function(x){quantile(x,0.9)}),
                        LowOuter=apply(Beta[,500:MCMC],1,f<-function(x){quantile(x,0.1)}),
                        Model="BARC"
                        )

custom.coefplot.default = function(custom.modelCI = my.ModelCI, title = "Coefficient Plot", xlab = "Value", 
          ylab = "Coefficient", innerCI = 1, outerCI = 2, lwdInner = 1, 
          lwdOuter = 0, pointSize = 3, color = "blue", shape = 16, 
          cex = 0.8, textAngle = 0, numberAngle = 0, zeroColor = "grey", 
          zeroLWD = 1, zeroType = 2, facet = FALSE, scales = "free", 
          sort = c("natural", "magnitude", "alphabetical"), decreasing = FALSE, 
          numeric = FALSE, fillColor = "grey", alpha = 1/2, horizontal = FALSE, 
          factors = NULL, only = NULL, shorten = TRUE, intercept = TRUE, 
          interceptName = "(Intercept)", coefficients = NULL, predictors = NULL, 
          strict = FALSE, newNames = NULL, plot = TRUE, ...) 
{
  theDots <- list(...)
  sort <- match.arg(sort)
  modelCI <- custom.modelCI
  if (!plot) {
    return(modelCI)
  }
  p <- coefplot:::buildPlotting.default(modelCI = modelCI, title = title, 
                             xlab = xlab, ylab = ylab, lwdInner = lwdInner, lwdOuter = lwdOuter, 
                             pointSize = pointSize, color = color, cex = cex, textAngle = textAngle, 
                             numberAngle = numberAngle, zeroColor = zeroColor, zeroLWD = zeroLWD, 
                             outerCI = outerCI, innerCI = innerCI, multi = FALSE, 
                             zeroType = zeroType, numeric = numeric, fillColor = fillColor, 
                             alpha = alpha, horizontal = horizontal, facet = facet, 
                             scales = scales)
  return(p)
}


custom.ci = my.ModelCI
custom.ci$Coefficient = factor(custom.ci$Coefficient, my.ModelCI$Coefficient[order(my.ModelCI$Value, decreasing = F)])
custom.ci = custom.ci[order(custom.ci$Value, decreasing = F), ]
custom.ci
custom.coefplot.default(custom.ci)

####################### PL-base method #######################

library(StatRank)

input.pl = showorder[,1,]

res.pl = Estimation.PL.MLE(input.pl)

p_vec = res.pl$Mean^-1

p_vec = p_vec/sum(p_vec)

print(cbind(p_vec,res.pl$order))
