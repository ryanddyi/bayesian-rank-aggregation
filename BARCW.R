
Gibbs_step_Y = function(ranker_id, Y_last, Beta_last, alpha_last, weights_last, exp_param=1){
  
  Y_new = Y_last[ranker_id,]
  
  Beta_cat = Beta_last
  alpha_cat = alpha_last
  weight = weights_last[ranker_id]
  
  for(j in 1:length(ranked_entities[[ranker_id]])){
    ranked = ranked_entities[[ranker_id]][[j]]
    
    intercept = exp_param*alpha_cat[ranked[1]]
    Y_new[ranked[1]] =  intercept + 
      betarightnorm(X[ranked[1],], exp_param*Beta_cat, exp_param/sqrt(weight), Y_new[ranked[2]]-intercept)
    
    for(m in 2:(nEach-1)){
      intercept = exp_param*alpha_cat[ranked[m]]
      Y_new[ranked[m]] = intercept + 
        betamidnorm(X[ranked[m],], exp_param*Beta_cat, exp_param/sqrt(weight), 
                    Y_new[ranked[m+1]]-intercept, Y_new[ranked[m-1]]-intercept)
    }
    
    intercept = exp_param*alpha_cat[ranked[nEach]]
    Y_new[ranked[nEach]] = intercept + 
      betaleftnorm(X[ranked[nEach],], exp_param*Beta_cat, exp_param/sqrt(weight), Y_new[ranked[nEach-1]]-intercept)
  }
  Y_new=Y_new/exp_param
  Y_new
}

Calculate_S_B = function(Y_last, alpha_last, weights_last, X){
  S = 0
  A = 0
  B = rep(0,p)
  
  for(i in 1:dim(Y_last)[1]){
    A = A+t(Y_last[i,]-alpha_last)%*%(Y_last[i,]-alpha_last)*weights_last[i]
    B = B+t(X)%*%(Y_last[i,]-alpha_last)*weights_last[i]
  }
  
  C = t(X)%*%X*sum(weights_last)
  H = solve(C)
  
  S = S+A-t(B)%*%H%*%B
  Beta_hat = H%*%B
  V_hat = H

  list(S=S, Beta_hat=Beta_hat, V_hat=V_hat)
}

Calculate_S = function(Y_last, weights_last, IX, Lambda_inv){
  
  S = 0
  A = 0
  B = rep(0,dim(IX)[2])
  
  for(i in 1:dim(Y_last)[1]){
    A = A+t(Y_last[i,])%*%Y_last[i,]*weights_last[i]
    B = B+t(IX)%*%Y_last[i,]*weights_last[i]
  }
  
  C = t(IX)%*%IX*sum(weights_last)
  H = solve(C+Lambda_inv)
  
  S = S+A-t(B)%*%H%*%B
  eta_hat = H%*%B
  Sigma_hat = H
  
  list(S=S, eta_hat=eta_hat, Sigma_hat=Sigma_hat) 
}

Gibbs_step_alpha = function(Y_last, Beta_last, alpha_last, weights_last, lambda){
  alpha_new = alpha_last
  N = dim(Y_last)[2]
 
  for(k in 1:N){
    alpha_new[k] = midnorm(sum((Y_last[,k]-X[k,]%*%Beta_last)*weights_last)/(sum(weights_last)+1/s0),
              sqrt(1/(sum(weights_last) + 1/s0)),-lambda,lambda)
  }
  alpha_new
}

WeightProb = function(weight,sumsq,num){
  prob=rep(0,length(weight))
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

Gibbs_step_weight = function(Y_last, Beta_last, alpha_last, weights_last, X){
  weights_new = weights_last
  allweight = c(0.5,1,2)
  for(i in 1:M){
    temp = sum((Y_last[i,]-X%*%Beta_last-alpha_last)^2)
    weights_new[i] = sample(allweight,size=1,prob=WeightProb(allweight,temp,N))
  }
  weights_new
}
