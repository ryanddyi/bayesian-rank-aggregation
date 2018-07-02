# orthodontics data - rank aggregation mixture model

setwd('/Users/ydd/Documents/BARC/code')
source("./truncated-normal.R")
source("./BARCM.R")
covariates = as.matrix(read.csv("data-read/covariates.csv"))

MCMC=5000
burn=1000

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

#========================BARCW=======================
source("./BARCW.R")
# hyper-param
lambda=5
s0=1

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
  if(n%%50==0){
    print(n)
    print(expansion[n])
    print(weight[,n])
  }
  for(i in 1:M){
    Y[i,,(n+1)] = Gibbs_step_Y(i, Y[,,n], Beta[,n], alpha[,n], weight[,n], expansion[n])
  }
  
  temp = Calculate_S_B(Y[,,n+1], alpha[,n], weight[,n], X)
  
  expansion[n+1] = sqrt(temp$S/rchisq(1,N*M))
  
  Beta[,n+1]=mvrnorm(1,temp$Beta_hat/expansion[n+1],temp$V_hat)
  
  alpha[,n+1] = Gibbs_step_alpha(Y[,,n+1], Beta[,n+1], alpha[,n], weight[,n], lambda)
  
  weight[,n+1] = Gibbs_step_weight(Y[,,n+1], Beta[,n+1], alpha[,n+1], weight[,n], X)
}

weight_est = rowMeans(weight[,burn:MCMC])

pi_order = order(pi_est, decreasing=T)

idCat4 = idCat_est_dp

for(i in 1:3){
  idCat4[idCat_est_dp==pi_order[i]] = i
}
idCat4[!idCat_est_dp%in%pi_order[1:3]] = 'others'

df_weight = cbind(idCat4,weight_est)
df_weight = as.data.frame(df_weight)

df_weight$weight_est = as.numeric(as.character(df_weight$weight_est))

p<-ggplot(df_weight, aes(x=idCat4, y=weight_est,group = idCat4))+
  theme_bw() +
  ggtitle("BARCW & BARCM -- Rankers' weights by cluster") +
  geom_boxplot()+
  labs(x = 'Clusters', y = 'Weight')
print(p)

pdf('weights-cat.pdf', height = 4, width = 6)
print(p)
dev.off()

################################

mu = X%*%Beta[,burn:MCMC]+alpha[,burn:MCMC]
for(n in 1:dim(mu)[2]){
  mu[,n] = mu[,n]-mean(mu[,n])
}

Y_agg_w = apply(mu,1,mean)

###############################
library("coefplot")

burnIn=500
my.ModelCI = data.frame(Value = apply(Beta[,burnIn:MCMC],1,mean), 
                        Coefficient = c("Upper right segment","Upper anterior segment","Upper left segment", "Lower right segment", "Lower anterior segment","Lower left segment","Right buccal occlusion","Left buccal occlusion","Overjet","Overbit","Centerline"),
                        HighInner=apply(Beta[,burnIn:MCMC],1,f<-function(x){quantile(x,0.975)}),
                        LowInner=apply(Beta[,burnIn:MCMC],1,f<-function(x){quantile(x,0.025)}),
                        HighOuter=apply(Beta[,burnIn:MCMC],1,f<-function(x){quantile(x,0.9)}),
                        LowOuter=apply(Beta[,burnIn:MCMC],1,f<-function(x){quantile(x,0.1)}),
                        Model="BARCW"
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


