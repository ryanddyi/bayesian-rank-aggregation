Gibbs_step_Y_DP = function(ranker_id, Y_last, Beta_last, alpha_last){
  
  Y_new = Y_last[ranker_id,]
  
  for(j in 1:length(ranked_entities[[ranker_id]])){
    ranked = ranked_entities[[ranker_id]][[j]]
    
    intercept = alpha_last[ranked[1]]
    Y_new[ranked[1]] =  intercept + 
      betarightnorm(X[ranked[1],], Beta_last, 1, Y_new[ranked[2]]-intercept)
    
    for(m in 2:(nEach-1)){
      intercept = alpha_last[ranked[m]]
      Y_new[ranked[m]] = intercept + 
        betamidnorm(X[ranked[m],], Beta_last, 1, Y_new[ranked[m+1]]-intercept, Y_new[ranked[m-1]]-intercept)
    }
    
    intercept = alpha_last[ranked[nEach]]
    Y_new[ranked[nEach]] = intercept + 
      betaleftnorm(X[ranked[nEach],], Beta_last, 1, Y_new[ranked[nEach-1]]-intercept)
  }
  Y_new
}

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

Gibbs_step_Y_CL = function(ranker_id, Y_last, idCat_last, Beta_last, alpha_last, exp_param=1){
  
  Y_new = Y_last[ranker_id,]
  
  # category of ranker i at step n
  category = idCat_last[ranker_id]
  Beta_cat = Beta_last[,category]
  alpha_cat = alpha_last[,category]
  
  for(j in 1:length(ranked_entities[[ranker_id]])){
    ranked = ranked_entities[[ranker_id]][[j]]
    
    intercept = exp_param*alpha_cat[ranked[1]]
    Y_new[ranked[1]] =  intercept + 
      betarightnorm(X[ranked[1],], exp_param*Beta_cat, exp_param, Y_new[ranked[2]]-intercept)
    
    for(m in 2:(nEach-1)){
      intercept = exp_param*alpha_cat[ranked[m]]
      Y_new[ranked[m]] = intercept + 
        betamidnorm(X[ranked[m],], exp_param*Beta_cat, exp_param, 
                    Y_new[ranked[m+1]]-intercept, Y_new[ranked[m-1]]-intercept)
    }
    
    intercept = exp_param*alpha_cat[ranked[nEach]]
    Y_new[ranked[nEach]] = intercept + 
      betaleftnorm(X[ranked[nEach],], exp_param*Beta_cat, exp_param, Y_new[ranked[nEach-1]]-intercept)
  }
  Y_new=Y_new/exp_param
  Y_new
}

Calculate_S_B = function(Y_last, idCat_last, alpha_last){
  S = 0
  Beta_hat = list() 
  V_hat = list()
  for(category in 1:dim(alpha_last)[2]){
    alpha_cat = alpha_last[,category]
    
    # all rankers in certain category
    ranker_cat = which(idCat_last==category)
    
    A = 0
    B = rep(0,p)
    
    for(i in ranker_cat){
      A = A+t(Y_last[i,]-alpha_cat)%*%(Y_last[i,]-alpha_cat)
      B = B+t(X)%*%(Y_last[i,]-alpha_cat)
    }
    
    C = t(X)%*%X*M
    H = solve(C)
    
    S = S+A-t(B)%*%H%*%B
    Beta_hat[[category]] = H%*%B
    V_hat[[category]] = H
  }
  list(S=S, Beta_hat=Beta_hat, V_hat=V_hat)
}

Calculate_S_B_DP = function(Y_last, idCat_last, alpha_last){
  S = 0
  Beta_hat = list() 
  V_hat = list()
  for(category in unique(idCat_last)){
    alpha_cat = alpha_last[,category]
    
    # all rankers in certain category
    ranker_cat = which(idCat_last==category)
    
    A = 0
    B = rep(0,p)
    
    for(i in ranker_cat){
      A = A+t(Y_last[i,]-alpha_cat)%*%(Y_last[i,]-alpha_cat)
      B = B+t(X)%*%(Y_last[i,]-alpha_cat)
    }
    
    C = t(X)%*%X*M
    H = solve(C)
    
    S = S+A-t(B)%*%H%*%B
    Beta_hat[[category]] = H%*%B
    V_hat[[category]] = H
  }
  list(S=S, Beta_hat=Beta_hat, V_hat=V_hat)
}

Calculate_S_DP = function(Y_last, idCat_last, IX, Lambda_inv){
  
  S = 0
  eta_hat = list() 
  Sigma_hat = list()

  for(category in unique(idCat_last)){

    # all rankers in certain category
    ranker_cat = which(idCat_last==category)
    
    A = 0
    B = rep(0,dim(IX)[2])
    
    for(i in ranker_cat){
      A = A+t(Y_last[i,])%*%Y_last[i,]
      B = B+t(IX)%*%Y_last[i,]
    }
    
    C = t(IX)%*%IX*length(ranker_cat)
    H = solve(C+Lambda_inv)
    
    S = S+A-t(B)%*%H%*%B
    eta_hat[[category]] = H%*%B
    Sigma_hat[[category]] = H
  }
  
  list(S=S, eta_hat=eta_hat, Sigma_hat=Sigma_hat) 
}

Gibbs_step_idCat_collapsed = function(Y_last, idCat_last, pi_last, H){

  idCat_new = idCat_last
  nComp = length(pi_last)
  M = length(idCat_last)
  for(i in 1:M){
    log_prob_temp = rep(0,nComp)
    
    # set ranker i to no category
    idCat_new[i] = 0
    
    for(category in 1:nComp){
      # all rankers in certain category
      ranker_cat = which(idCat_new==category)
      cat_size = length(ranker_cat)

      Y0 = rep(0,dim(Y_last)[2])
      
      for(ranker in ranker_cat){
        Y0 = Y0 + Y_last[ranker,]
      }
      
      Y1 = Y0 + Y_last[i,]
            
      S_diff = 0
      if(cat_size>0){
        S_diff = S_diff + t(Y0)%*%Y0/(cat_size+1) + t(Y0)%*%H%*%Y0/(cat_size+1)/cat_size 
      }
      S_diff = S_diff - t(Y1)%*%Y1/(cat_size+2) - t(Y1)%*%H%*%Y1/(cat_size+2)/(cat_size+1)
    
      #print(S_diff)
      
      log_prob_temp[category] = -S_diff/2
    }
    #print(pi_last*exp(log_prob_temp-max(log_prob_temp)))
    idCat_new[i] = sample(nComp,1,prob = pi_last*exp(log_prob_temp-max(log_prob_temp)))
    #print(idCat_new)
  }
  idCat_new
}

Gibbs_step_idCat_DP = function(Y_last, idCat_last, H, dp_count=1){
  
  idCat_new = idCat_last
  
  M = length(idCat_last)
  for(i in 1:M){
    log_prob_temp = rep(0,M)
    count_temp = rep(0,M)
    
    # set ranker i to no category
    idCat_new[i] = 0
    
    for(category in unique(idCat_new[-i])){
      # all rankers in certain category
      ranker_cat = which(idCat_new==category)
      cat_size = length(ranker_cat)
      
      Y0 = rep(0,dim(Y_last)[2])
      
      for(ranker in ranker_cat){
        Y0 = Y0 + Y_last[ranker,]
      }
      
      Y1 = Y0 + Y_last[i,]
      
      if(cat_size>0) {
        log_prob_temp[category] = 0.5*(t(Y1)%*%H[[cat_size+1]]%*%Y1-t(Y0)%*%H[[cat_size]]%*%Y0)
      } else {
        log_prob_temp[category] = 0.5*t(Y1)%*%H[[cat_size+1]]%*%Y1
      }
      count_temp[category] = cat_size
    }
    empty_cat = min(which(count_temp==0))
    count_temp[empty_cat] = dp_count
    
    Y1 = Y_last[i,]
    cat_size = 0
    log_prob_temp[empty_cat] = 0.5*t(Y1)%*%H[[cat_size+1]]%*%Y1
    
    #print(empty_cat)
    #print(log_prob_temp)
    idCat_new[i] = sample(M, 1, prob = count_temp*exp(log_prob_temp-max(log_prob_temp)))
  }
  idCat_new
}

Gibbs_step_pi = function(pi_last, idCat_last, ps_count=2){
  nComp = length(pi_last)
  count_comps = rep(0,nComp)
  for(id in 1:nComp){
    count_comps[id] = sum(idCat_last==id)
  }
  rdirichlet(1, ps_count/nComp + count_comps)
}
