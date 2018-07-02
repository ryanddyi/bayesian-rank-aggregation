#functions we use
library("MCMCpack")

#sample from right truncated standard normal
rightsdnorm<-function(c){
	lambda=(c+sqrt(c^2+4))/2
	accept=0
	while(accept==0){	
		x=rexp(1,lambda)
		g=dnorm(x+c,log=T)-((lambda^2-2*lambda*c)/2-1/2*log(2*pi)-lambda*x)
		if(log(runif(1,0,1))<=g) accept=1
	}
	x+c
}

#sample from right-truncated normal
rightnorm<-function(mu, sigma, cut){
	if(cut>=mu){
		h=rightsdnorm((cut-mu)/sigma)
		x=h*sigma+mu
	}
	if(cut<mu){
		accept=0
		while(accept==0){
			x=rnorm(1,mu,sigma)
			if(x>=cut){accept=1}
		}
	}
	x
}

#sample from left-truncated normal
leftnorm<-function(mu, sigma, cut){
	h=rightnorm(mu, sigma, 2*mu-cut)
	x=2*mu-h
}

#sample from mid-truncated normal
midnorm<-function(mu, sigma, lcut, rcut){
	if(lcut>rcut) print("invalid input")
	if(mu>=lcut&mu<=rcut){h=dnorm(mu,mu,sigma)}
	if(mu<lcut){h=dnorm(lcut,mu,sigma)}
	if(mu>rcut){h=dnorm(rcut,mu,sigma)}
	accept=0
	while(accept==0){
		x=runif(1,lcut,rcut)
		g=dnorm(x,mu,sigma,log=T)-log(h)
		if(log(runif(1,0,1))<=g){accept=1}
	}
	x
}

betarightnorm<-function(x, beta, delta, cut){
	d=rightnorm(x%*%beta, delta, cut)
	d	
}

betaleftnorm<-function(x, beta, delta, cut){
	d=leftnorm(x%*%beta, delta, cut)
	d	
}

betamidnorm<-function(x, beta, delta, lcut, rcut){
	d=midnorm(x%*%beta, delta, lcut, rcut)
	d	
}

midchisq<-function(df,scale,lcut,rcut){
  if(df<=2){print("not applicable")}
  mode=(df-2)/scale
  if(lcut<rcut){
    if(mode>=lcut&mode<=rcut){logh=dchisq(mode*scale,df,log=T)}
    if(mode<lcut){logh=dchisq(lcut*scale,df,log=T)}
    if(mode>rcut){logh=dchisq(rcut*scale,df,log=T)}
    accept=0
    while(accept==0){
      x=runif(1,lcut,rcut)
      g=dchisq(x*scale,df,log=T)-logh
      if(log(runif(1,0,1))<=g){accept=1}
    }
  }
  if(lcut==rcut){x=lcut}
  x
}
