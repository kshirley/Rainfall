########################################################################################################################

# Functions for MCMC from multisite probit (frequency) model:

# new way to draw z:
draw.z.ind <- function(bounds,w,alpha,tau,site.mat,arc){
  T <- dim(w)[1]
  L.sum <- dim(bounds)[1]
  z.all <- matrix(NA,T,L.sum)
  for (l in 1:L.sum){
  	z.new <- rep(0,T)
  	mean <- w[,site.mat[l,1]]+alpha[site.mat[l,1]]*arc[l]
    z.score <- mean/tau[site.mat[l,1]]
    w1 <- which(z.score < -3.2 & bounds[l,,1]==0) # low z.score and truncated on the left by zero:
    w2 <- which(z.score > 3.2 & bounds[l,,2]==0) # high z.score and truncated on the right by zero:
    w3 <- 1:T
    if (length(w1)+length(w2)>0) w3 <- (1:T)[-c(w1,w2)]
    if (length(w3)>0) z.new[w3] <- rtnorm(n=length(w3),mean=mean[w3],sd=tau[site.mat[l,1]],lower=bounds[l,w3,1],upper=bounds[l,w3,2])
    if (sum(is.na(z.new))>0) print(paste("Iteration",g,"NA in z w3, site:",l))
    if (sum(z.new==Inf)>0) print(paste("Iteration",g,"Inf in z w3, site:",l))
    if (sum(z.new==-Inf)>0) print(paste("Iteration",g,"-Inf in z w3, site:",l))
    if (length(w1)>0) z.new[w1] <- rtnorm.exp(n=length(w1),mean=mean[w1],sd=tau[site.mat[l,1]],z.lower=-z.score[w1],right=0,iter.max=100)
    if (sum(is.na(z.new))>0) print(paste("Iteration",g,"NA in z w1, site:",l))
    if (sum(z.new==Inf)>0) print(paste("Iteration",g,"Inf in z w1, site:",l))
    if (sum(z.new==-Inf)>0) print(paste("Iteration",g,"-Inf in z w1, site:",l))
    if (length(w2)>0) z.new[w2] <- rtnorm.exp(n=length(w2),mean=mean[w2],sd=tau[site.mat[l,1]],z.lower=z.score[w2],right=1,iter.max=100)
    if (sum(is.na(z.new))>0) print(paste("Iteration",g,"NA in z w2, site:",l))
    if (sum(z.new==Inf)>0) print(paste("Iteration",g,"Inf in z w2, site:",l))
    if (sum(z.new==-Inf)>0) print(paste("Iteration",g,"-Inf in z w2, site:",l))
    z.all[,l] <- z.new
  }
  return(z.all)
}

# draw alpha:
alpha.draw <- function(S,z,site.mat,w,arc,T,tau,mu.alpha,tau.alpha,L,ws){
  new.alpha <- numeric(S)
  b0 <- mu.alpha; sig0 <- matrix(tau.alpha^2,1,1)
  for (s in 1:S){
  	ym <- matrix(z[,ws[[s]]] - w[,s],ncol=1)
    xm <- matrix(x.arc[,ws[[s]]],ncol=1)
    v.post <- solve(t(xm)%*%xm+sig0)  	
    b.hat <- v.post %*% (t(xm)%*%ym+sig0%*%b0)
    new.alpha[s] <- rnorm(1,b.hat,tau[s]*sqrt(v.post))
  }
  return(new.alpha)  	
}

# draw tau.alpha with a metropolis step:
ta.draw <- function(alpha,mu.alpha,tau.alpha,jump.ta,g,adapt){
  ll.data <- sum(log(dnorm(alpha,mu.alpha,tau.alpha)))
  lp.old <- ll.data + log(dst(tau.alpha,nu=3,mu=0,sigma=0.5))
  tau.alpha.star <- abs(tau.alpha + rnorm(1,0,jump.ta))
  ll.data.star <- sum(log(dnorm(alpha,mu.alpha,tau.alpha.star)))
  lp.star <- ll.data.star + log(dst(tau.alpha.star,nu=3,mu=0,sigma=0.5))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.ta <- ifelse(r > 0.44,jump.ta*A1,jump.ta*B1)
  #if (runif(1) < 0.44) tau.alpha <- tau.alpha.star
  if (runif(1) < r) tau.alpha <- tau.alpha.star
  return(list(tau.alpha=tau.alpha,jump.ta=jump.ta))
}

# draw tau.alpha with a metropolis step:
tau.draw <- function(tau,z,w,alpha,arc,site.mat,T,jump.tau,g,adapt,L,ws,x.arc){
  new.tau <- tau
  for (s in 1:S){
  	ns <- L[s]
    ll.data <- sum(dnorm(z[,ws[[s]]],mean=rep(w[,s],L[s])+alpha[s]*x.arc[,ws[[s]]],sd=tau[s],log=TRUE))
    lp.old <- ll.data + log(dst(tau[s],nu=3,mu=0,sigma=1))
    tau.star <- abs(tau[s] + rnorm(1,0,jump.tau[s]))
    ll.data.star <- sum(dnorm(z[,ws[[s]]],mean=rep(w[,s],L[s])+alpha[s]*x.arc[,ws[[s]]],sd=tau.star,log=TRUE))
    lp.star <- ll.data.star + log(dst(tau.star,nu=3,mu=0,sigma=1))
    r <- exp(lp.star - lp.old)
    if (g < adapt) jump.tau[s] <- ifelse(r > 0.44,jump.tau[s]*A1,jump.tau[s]*B1)
    #if (runif(1) < 0.44) new.tau[s] <- tau.star
	if (runif(1) < r) new.tau[s] <- tau.star
  }
  return(list(tau=new.tau,jump.tau=jump.tau))
}

# draw w:
#w.draw <- function(T,S,z,lambda,d.mat,L,tau,x.mat,beta,ws,alpha.expand){
#  z.bar.mat <- t(apply(rbind(rep(0,T),apply(z-alpha.expand,1,cumsum)[c(2,4,6,8,10,15),]),2,diff)/L)
#  R.inv <- solve(R.cov(lambda,d.mat))
#  SigN <- solve(R.inv + diag(L/tau^2))
#  muN <- SigN %*% (R.inv%*%t(x.mat%*%beta) + diag(L/tau^2)%*%t(z.bar.mat))
#  new.w <- mvrnorm(T,rep(0,S),SigN) + t(muN)
#  return(new.w)
#}

# draw lambda with a metropolis step:
lambda.draw <- function(x.mat,beta,lambda,d.mat,w,S,jump.lambda){
  w.mu <- x.mat %*% beta
  R <- R.cov(lambda,d.mat)
  ll.data <- sum(dmvnorm(w-w.mu,rep(0,S),R,log=TRUE))
  lp.old <- ll.data + log(dgamma(lambda,shape=50,scale=0.03))
  lambda.star <- abs(lambda + rnorm(1,0,jump.lambda))
  R.star <- R.cov(lambda.star,d.mat)
  ll.data.star <- sum(dmvnorm(w-w.mu,rep(0,S),R.star,log=TRUE)) 
  lp.star <- ll.data.star + log(dgamma(lambda.star,shape=50,scale=0.03))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.lambda <- ifelse(r > 0.44,jump.lambda*A1,jump.lambda*B1)
  #if (runif(1) < 0.44) lambda <- lambda.star
  if (runif(1) < r) lambda <- lambda.star
  return(list(lambda=lambda,jump.lambda=jump.lambda))
}

# draw beta given z, R:
draw.beta <- function(R,vec,xtx.ones,z,x.mat,sigma,mu){
  S <- dim(R)[1]
  P <- dim(beta)[1]
  R.inv <- solve(R)
  B <- matrix(R.inv[vec]*as.numeric(xtx.ones),S*P,S*P)
  W <- matrix(0,P,S)   
  for (j in 1:P){
    for (s in 1:S){
      temp <- 0
      for (l in 1:S) temp <- temp + R.inv[s,l]*sum(z[l,]*x.mat[,j])
      W[j,s] <- W[j,s] + temp
    }
  }
  W <- matrix(W,ncol=1)
  sig.prior <- diag(1/rep(sigma^2,S))
  B <- B + sig.prior
  B.inv <- solve(B)
  beta.prior <- rep(mu,S)
  beta.hat <- B.inv %*% (W + sig.prior %*% beta.prior)    
  beta <- matrix(mvrnorm(1,beta.hat,B.inv),P,S)
  return(list(beta=beta,B=B))
}


# script for metropolis-hastings draw of mean and sds of beta:
mu.sigma.draw <- function(beta,mu,sigma,S,P,jump.sigma,g,adapt){
  ll.data <- apply(matrix(dnorm(beta,mu,sigma,log=TRUE),ncol=S),1,sum)
  lp.old <- ll.data + log(dst(sigma,nu=3,mu=0,sigma=1))
  sigma.star <- abs(sigma+rnorm(P,0,jump.sigma))
  ll.data.star <- apply(matrix(dnorm(beta,mu,sigma.star,log=TRUE),ncol=S),1,sum)    
  lp.star <- ll.data.star + log(dst(sigma.star,nu=3,mu=0,sigma=1))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.sigma <- ifelse(r > 0.44,jump.sigma*A1,jump.sigma*B1)
  promote <- runif(P) < r
  sigma[promote] <- sigma.star[promote]
  mu <- rnorm(P,apply(beta,1,mean),sigma/sqrt(S))
  return(list(mu=mu,sigma=sigma,jump.sigma=jump.sigma))
}


########################################################################################################################

# sample from a truncated normal on one side using an exponential rejection sampling scheme:
rtnorm.exp <- function(n=1,mean=0,sd=1,z.lower=-Inf,right=1,iter.max=10){
  a.star <- (z.lower + sqrt(z.lower^2+4))/2
  ret.val <- rep(NA,n)
  w <- 1:n
  for (i in 1:iter.max){
    e <- rexp(n,rate=a.star) + z.lower
    rho <- exp(-0.5*(e-a.star)^2)
    u <- runif(n) < rho
    ret.val[w][which(u[w])] <- e[w][which(u[w])]
    n.new <- sum(!u[w])
    if (n.new==0){    
      break
    } else {
      #w <- which(!u[w])
      w <- which(is.na(ret.val))
    }
  }
  return(ret.val*(-1)^right*sd+mean)
}

# create exponential covariance matrix:
R.cov <- function(lambda,d){
  R <- diag(dim(d)[1])
  R[upper.tri(R)] <- exp(-lambda*d[upper.tri(d)])
  R[lower.tri(R)] <- exp(-lambda*d[lower.tri(d)])
  return(R)
}


########################################################################################################################

# Functions for multisite hierarchical tobit model:

# new way to draw z:
draw.W <- function(R,Z,beta.arc,Sigma,gamma.mat,X.arc,W,up,w.draw){
  S <- length(R)
  T <- dim(R[[1]])[2]
  J <- unlist(lapply(R,dim))[seq(1,by=2,length=S)]
  for (s in 1:S){
    mu.mat <- t(matrix(Z[,s],T,J[s]) + matrix(rep(X.arc[[s]]*beta.arc[s],each=T),ncol=J[s]))  # J[s] x T matrix
    for (j in 1:J[s]){ # draw j given -j
      Sigma.temp <- Sigma[[s]][j,-j] %*% solve(Sigma[[s]][-j,-j])
      cond.mu <- mu.mat[j,] + Sigma.temp %*% (W[[s]][-j,] - mu.mat[-j,])
      cond.Sigma <- (Sigma[[s]][j,j] - Sigma.temp %*% Sigma[[s]][-j,j])/gamma.mat[,s]
      z.score <- cond.mu/as.numeric(sqrt(cond.Sigma))
      sel <- w.draw[[s]][j,]
      w2 <- which(z.score[sel] > 3.2 & up[[s]][j,sel]==0) # high z.score and truncated on the right by zero:
      w3 <- 1:sum(sel)
      if (length(w2)>0) w3 <- w3[-w2]
      if (length(w3)>0) W[[s]][j,sel][w3] <- rtnorm(n=length(w3),mean=cond.mu[sel][w3],sd=sqrt(cond.Sigma[sel][w3]),lower=-Inf,upper=up[[s]][j,sel][w3])
      if (sum(is.na(W[[s]][j,]))>0) print(paste("Iteration",g,"NA in W w3"))
      if (sum(W[[s]][j,]==Inf)>0) print(paste("Iteration",g,"Inf in W w3"))
      if (sum(W[[s]][j,]==-Inf)>0) print(paste("Iteration",g,"-Inf in W w3"))
      if (length(w2)>0) W[[s]][j,sel][w2] <- rtnorm.exp(n=length(w2),mean=cond.mu[sel][w2],sd=sqrt(cond.Sigma[sel][w2]),z.lower=z.score[sel][w2],right=1,iter.max=100)
      if (sum(is.na(W[[s]][j,]))>0) print(paste("Iteration",g,"NA in W w2"))
      if (sum(W[[s]][j,]==Inf)>0) print(paste("Iteration",g,"Inf in W w2"))
      if (sum(W[[s]][j,]==-Inf)>0) print(paste("Iteration",g,"-Inf in W w2"))
    }
  }
  return(W)
}


# new way to draw z:
draw.W.start <- function(R,Z,beta.arc,Sigma,gamma.mat,X.arc,up,w.draw,W){
  S <- length(R)
  T <- dim(R[[1]])[2]
  J <- unlist(lapply(R,dim))[seq(1,by=2,length=S)]
  for (s in 1:S){
    mu.mat <- t(matrix(Z[,s],T,J[s]) + t(matrix(X.arc[[s]]*beta.arc[s],J[s],T)))  # J[s] x T matrix
    for (j in 1:J[s]){ # draw j given -j
      W[[s]][j,] <- R[[s]][j,]
      cond.mu <- mu.mat[j,]
      cond.Sigma <- rep(Sigma[[s]][j,j],T)
      z.score <- cond.mu/as.numeric(sqrt(cond.Sigma))
      sel <- w.draw[[s]][j,]
      w2 <- which(z.score[sel] > 3.2 & up[[s]][j,sel]==0) # high z.score and truncated on the right by zero:
      w3 <- 1:sum(sel)
      if (length(w2)>0) w3 <- w3[-w2]
      if (length(w3)>0) W[[s]][j,sel][w3] <- rtnorm(n=length(w3),mean=cond.mu[sel][w3],sd=sqrt(cond.Sigma[sel][w3]),lower=-Inf,upper=up[[s]][j,sel][w3])
      if (sum(is.na(W[[s]][j,]))>0) print(paste("Iteration",g,"NA in W w3"))
      if (sum(W[[s]][j,]==Inf)>0) print(paste("Iteration",g,"Inf in W w3"))
      if (sum(W[[s]][j,]==-Inf)>0) print(paste("Iteration",g,"-Inf in W w3"))
      if (length(w2)>0) W[[s]][j,sel][w2] <- rtnorm.exp(n=length(w2),mean=cond.mu[sel][w2],sd=sqrt(cond.Sigma[sel][w2]),z.lower=z.score[sel][w2],right=1,iter.max=100)
      if (sum(is.na(W[[s]][j,]))>0) print(paste("Iteration",g,"NA in W w2"))
      if (sum(W[[s]][j,]==Inf)>0) print(paste("Iteration",g,"Inf in W w2"))
      if (sum(W[[s]][j,]==-Inf)>0) print(paste("Iteration",g,"-Inf in W w2"))
    }
  }
  return(W)
}




# rwishart from package bayesm:
rwishart <- function (nu, V){
	m = nrow(V)
	df = (nu + nu - m + 1) - (nu - m + 1):nu
	if (m > 1) {
		T = diag(sqrt(rchisq(c(rep(1, m)), df)))
		T[lower.tri(T)] = rnorm((m * (m + 1)/2 - m))
	}
	else {
		T = sqrt(rchisq(1, df))
	}
	U = chol(V)
	C = t(T) %*% U
	CI = backsolve(C, diag(m))
	return(list(W = crossprod(C), IW = crossprod(t(CI)), C = C, 
							CI = CI))
}


# draw alpha with a metropolis step:
alpha.draw <- function(gamma.mat,a,b,jump.alpha){
  ll.data <- sum(dgamma(gamma.mat,shape=alpha/2,scale=2/alpha,log=TRUE))
  lp.old <- ll.data + dgamma(alpha,shape=a,scale=b,log=TRUE)
  alpha.star <- abs(alpha + rnorm(1,0,jump.alpha))
  ll.data.star <- sum(dgamma(gamma.mat,shape=alpha.star/2,scale=2/alpha.star,log=TRUE))
  lp.star <- ll.data.star + dgamma(alpha.star,shape=a,scale=b,log=TRUE)
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.alpha <- ifelse(r > 0.44,jump.alpha*A1,jump.alpha*B1)
  if (runif(1) < r) alpha <- alpha.star
  return(list(alpha=alpha,jump.alpha=jump.alpha))
}

# draw alpha with a metropolis step:
alpha.draw.unif <- function(gamma.mat,a,b,jump.alpha){
  ll.data <- sum(dgamma(gamma.mat,shape=alpha/2,scale=2/alpha,log=TRUE))
  lp.old <- ll.data + dunif(alpha,a,b,log=TRUE)
  alpha.star <- alpha + rnorm(1,0,jump.alpha)
  if (alpha.star > a & alpha.star < b){
    ll.data.star <- sum(dgamma(gamma.mat,shape=alpha.star/2,scale=2/alpha.star,log=TRUE))
    lp.star <- ll.data.star + dunif(alpha,a,b,log=TRUE)
    r <- exp(lp.star - lp.old)
  } else {
    r <- 0
  }
  if (g < adapt) jump.alpha <- ifelse(r > 0.44,jump.alpha*A1,jump.alpha*B1)
  if (runif(1) < r) alpha <- alpha.star
  return(list(alpha=alpha,jump.alpha=jump.alpha))
}


# draw beta given z, R:
draw.beta.tobit <- function(tau,lambda,vec,xtx.ones,Z,X,sigma,mu){
  Covar <- tau^2*R.cov(lambda,d.mat)
  S <- dim(Covar)[1]
  P <- dim(beta)[1]
  Covar.inv <- solve(Covar)
  B <- matrix(Covar.inv[vec]*as.numeric(xtx.ones),S*P,S*P)
  W.temp <- matrix(0,P,S)   
  for (j in 1:P){
    for (s in 1:S){
      temp <- 0
      for (l in 1:S) temp <- temp + Covar.inv[s,l]*sum(Z[,l]*X[,j])
      W.temp[j,s] <- W.temp[j,s] + temp
    }
  }
  W.temp <- matrix(W.temp,ncol=1)
  sig.prior <- diag(1/rep(sigma^2,S))
  B <- B + sig.prior
  B.inv <- solve(B)
  beta.prior <- rep(mu,S)
  beta.hat <- B.inv %*% (W.temp + sig.prior %*% beta.prior)    
  beta <- matrix(mvrnorm(1,beta.hat,B.inv),P,S)
  return(list(beta=beta,B=B))
}

tau.draw.tobit <- function(Z,beta,lambda,a,b){
  b.post <- b + 0.5*sum(((Z-X%*%beta) %*% solve(R.cov(lambda,d.mat)))*(Z-X%*%beta))
  a.post <- a + T*S*0.5
  return(sqrt(1/rgamma(1,shape=a.post,rate=b.post)))
}


# draw lambda with a metropolis step:
lambda.draw.tobit <- function(X,beta,tau,lambda,d.mat,Z,S,jump.lambda){
  mu.Z <- X %*% beta
  Covar <- tau^2*R.cov(lambda,d.mat)
  ll.data <- sum(dmvnorm(Z-mu.Z,rep(0,S),Covar,log=TRUE))
  lp.old <- ll.data + log(dgamma(lambda,shape=50,scale=0.03))
  lambda.star <- abs(lambda + rnorm(1,0,jump.lambda))
  Covar.star <- tau^2*R.cov(lambda.star,d.mat)
  ll.data.star <- sum(dmvnorm(Z-mu.Z,rep(0,S),Covar.star,log=TRUE)) 
  lp.star <- ll.data.star + log(dgamma(lambda.star,shape=50,scale=0.03))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.lambda <- ifelse(r > 0.44,jump.lambda*A1,jump.lambda*B1)
  if (runif(1) < r) lambda <- lambda.star
  return(list(lambda=lambda,jump.lambda=jump.lambda))
}

# draw gamma:
draw.gamma <- function(S,Z,X.arc,beta.arc,T,J,Sigma,W,alpha,gamma.temp){
  for (s in 1:S){      
    mu.W <- t(Z[,s] + matrix(rep(X.arc[[s]]*beta.arc[s],each=T),ncol=J[s]))
    S.inv <- solve(Sigma[[s]])
    lhs <- t(W[[s]] - mu.W)
    out <- .C("draw_gamma",params=as.integer(c(J[s],T)),Sinv_in=as.double(S.inv),lhs_in=as.double(lhs),rate_in=as.double(rep(0,T)),alpha_in=as.double(alpha))
    rate.vec <- out$rate_in
    gamma.temp[,s] <- rgamma(T,shape=(J[s]+alpha)/2,rate=rate.vec)
  }
  return(gamma.temp)
}


# Draw Z | W, gamma, beta.arc, Sigma, beta, lambda, tau:
draw.Z <- function(S,Sigma,T,J,X.arc,beta.arc,W,lambda,d.mat,tau,X,beta,gamma.mat){
  # Compute mean vector for each draw, from the data:
  sigma.z <- numeric(S)
  for (s in 1:S) sigma.z[s] <- sum(solve(Sigma[[s]]))
  mu.z <- matrix(0,S,T)
  for (s in 1:S) mu.z[s,] <- (matrix(-1,1,J[s]) %*% solve(Sigma[[s]]) %*% (matrix(X.arc[[s]]*beta.arc[s],J[s],T) - W[[s]]))/sigma.z[s]
  # Compute prior covariance:
  prior.inv <- solve(R.cov(lambda,d.mat))/tau^2
  # Compute prior mean of multivariate z (from beta):
  Pxb <- prior.inv %*% t(X%*%beta)
  # Compute covariance matrix for each time point (an (S x S x T) array?)
  Sig <- t(gamma.mat)*sigma.z # Here, Sig is the inverse of Sigma^z_t, i.e. the data precision matrix
  # Multiply Sig by mu.z to get the left addition element in the right hand side of the mean of z_t:
  Sig.mu <- Sig*mu.z
  z.draw <- .C("zdraw",params=as.integer(c(S,T)),Smu_in=as.double(Sig.mu),Pxb_in=as.double(Pxb),Sig_in=as.double(Sig),Prior_in=as.double(prior.inv))
  Z.out <- t(matrix(z.draw$Smu_in,S,T))
  return(Z.out)
}


# Draw beta.arc | W, gamma, Z, Sigma, mu.arc, tau.arc:
draw.beta.arc <- function(S,Sigma,gamma.mat,X.arc,Z,T,J,W,tau.arc,mu.arc){
  beta.arc.out <- numeric(S)
  mu.ba <- numeric(S); sigma.ba <- numeric(S)
  for (s in 1:S){
    solve.S <- solve(Sigma[[s]])
    denom <- sum(gamma.mat[,s])*t(X.arc[[s]]) %*% solve.S %*% X.arc[[s]]
    numer <- -t(X.arc[[s]]) %*% solve.S %*% apply(gamma.mat[,s]*(matrix(Z[,s],T,J[s])-t(W[[s]])),2,sum)
    mu.w <- numer/denom
    ssq.w <- 1/denom
    d2 <- 1/ssq.w + 1/tau.arc^2
    mu.ba[s] <- (mu.w/ssq.w + mu.arc/tau.arc^2)/d2
    sigma.ba[s] <- sqrt(1/d2)
    beta.arc.out[s] <- rnorm(1,mu.ba[s],sigma.ba[s])
  }
  return(beta.arc.out)
}  


# Sample Sigma | W, mu (Z, X.arc, beta.arc), gamma:
draw.Sigma <- function(S,Z,X.arc,beta.arc,T,J,W,gamma.mat,v.0,Lambda.0.inv,Sigma.null){
  for (s in 1:S){
    mu.W <- t(Z[,s] + matrix(rep(X.arc[[s]]*beta.arc[s],each=T),ncol=J[s]))      
    SS <- (t(t(W[[s]])*gamma.mat[,s]) - t(t(mu.W)*gamma.mat[,s])) %*% t(W[[s]]-mu.W)
  	Sigma.null[[s]] <- rwishart(nu = v.0[s] + T, V = solve(Lambda.0.inv[[s]] + SS))$IW
  }
  return(Sigma.null)
}



########################################################################################################################

# Old MCMC functions and scripts for correlated frequency model (single series per site):

# script for metropolis-hastings draw of mean and sds of beta:
rwmh.mean.sd <- quote({
  ll.data <- apply(matrix(dnorm(beta,b,sig,log=TRUE),ncol=S),1,sum)
  #lp.old <- ll.data + dunif(sig,0,3,log=TRUE)
  lp.old <- ll.data + log(dst(sig,nu=3,mu=0,sigma=1))
  sig.star <- abs(sig+rnorm(P,0,jump))
  ll.data.star <- apply(matrix(dnorm(beta,b,sig.star,log=TRUE),ncol=S),1,sum)    
  #lp.star <- ll.data.star + dunif(sig.star,0,3,log=TRUE)
  lp.star <- ll.data.star + log(dst(sig.star,nu=3,mu=0,sigma=1))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump <- ifelse(r > 0.44,jump*A1,jump*B1)
  promote <- runif(P) < r
  sig[promote] <- sig.star[promote]
  jump.gibbs[k,g,] <- jump
  sig.gibbs[k,g,] <- sig
  b <- rnorm(P,apply(beta,1,mean),sig/sqrt(S))
  b.gibbs[k,g,] <- b    
})

# script for draw of mean and sd using normal-inv-chisq prior:
nics.mean.sd <- quote({
  beta.bar <- apply(beta,1,mean)
  beta.ss <- apply((beta-beta.bar)^2,1,sum)
  draw <- nics(S,beta.bar,beta.ss,k0,mu0,v0,sig.sq0)
  b <- draw$mu
  sig <- sqrt(draw$sigma.squared)
  b.gibbs[k,g,] <- b
  sig.gibbs[k,g,] <- sig
})

# function for posterior of normal-inverse-chi-square model:
nics <- function(n,y.bar,s.sq,k0,mu0,v0,sig.sq0){
  mu.n <- k0*mu0/(k0+n) + n*y.bar/(k0+n)
  k.n <- k0 + n
  v.n <- v0 + n
  sig.sq.n <- (v0*sig.sq0 + (n-1)*s.sq + k0*n*(y.bar-mu0)^2/(k0+n))/v.n
  sigma.squared <- rics(length(y.bar),v.n,sig.sq.n)
  mu <- rnorm(length(y.bar),mu.n,sqrt(sigma.squared/k.n))
  return(list(mu=mu,sigma.squared=sigma.squared))
}

# get function to evaluate determinant of sigma replacing one value:
get.det <- function(R,r,i,j){
  new.R <- R
  new.R[i,j] <- r; new.R[j,i] <- r
  return(det(new.R))
}

# function to solve for interval of positive definite correlations:
get.int <- function(R,i,j){
  a <- (get.det(R,1,i,j) + get.det(R,-1,i,j) - 2*get.det(R,0,i,j))/2
  b <- (get.det(R,1,i,j) - get.det(R,-1,i,j))/2
  c <- get.det(R,0,i,j)
  return(sort((-b +c(-1,1)*sqrt(b^2-4*a*c))/(2*a)))
}

# impute the zs given y, beta, and R:
draw.z <- function(R,x.mat,beta,z.old,bounds){
  N <- dim(z.old)[2]
  S <- dim(R)[1]
  mu <- x.mat %*% beta
  for (s in 1:S){
    R.temp <- R[s,-s]%*%solve(R[-s,-s])
    cond.mu <- matrix(mu[,s],nrow=1) + R.temp%*%(z.old[-s,]-t(mu[,-s]))
    cond.R <- R[s,s] - R.temp%*%R[-s,s]
    z.old[s,] <- rtnorm(N,mean=cond.mu,sd=sqrt(cond.R),lower=bounds[s,,1],upper=bounds[s,,2])
  }
  return(z.old)
}

# new way to draw z:
draw.z.new <- function(R,x.mat,beta,z.old,bounds){
  N <- dim(z.old)[2]
  S <- dim(R)[1]
  mu <- x.mat %*% beta
  for (s in 1:S){
    R.temp <- R[s,-s]%*%solve(R[-s,-s])
    cond.mu <- matrix(mu[,s],nrow=1) + R.temp%*%(z.old[-s,]-t(mu[,-s]))
    cond.R <- R[s,s] - R.temp%*%R[-s,s]
    z.score <- cond.mu/as.numeric(sqrt(cond.R))
    w1 <- which(z.score < -3.2 & bounds[s,,1]==0) # low z.score and truncated on the left by zero:
    w2 <- which(z.score > 3.2 & bounds[s,,2]==0) # high z.score and truncated on the right by zero:
    w3 <- 1:N
    if (length(w1)+length(w2)>0) w3 <- (1:N)[-c(w1,w2)]
    if (length(w3)>0) z.old[s,w3] <- rtnorm(n=length(w3),mean=cond.mu[1,w3],sd=sqrt(cond.R),lower=bounds[s,w3,1],upper=bounds[s,w3,2])
    if (sum(is.na(z.old))>0) print(paste("Iteration",g,"NA in z w3"))
    if (sum(z.old==Inf)>0) print(paste("Iteration",g,"Inf in z w3"))
    if (sum(z.old==-Inf)>0) print(paste("Iteration",g,"-Inf in z w3"))
    if (length(w1)>0) z.old[s,w1] <- rtnorm.exp(n=length(w1),mean=cond.mu[1,w1],sd=sqrt(cond.R),z.lower=-z.score[w1],right=0,iter.max=100)
    if (sum(is.na(z.old))>0) print(paste("Iteration",g,"NA in z w1"))
    if (sum(z.old==Inf)>0) print(paste("Iteration",g,"Inf in z w1"))
    if (sum(z.old==-Inf)>0) print(paste("Iteration",g,"-Inf in z w1"))
    if (length(w2)>0) z.old[s,w2] <- rtnorm.exp(n=length(w2),mean=cond.mu[1,w2],sd=sqrt(cond.R),z.lower=z.score[w2],right=1,iter.max=100)
    if (sum(is.na(z.old))>0) print(paste("Iteration",g,"NA in z w2"))
    if (sum(z.old==Inf)>0) print(paste("Iteration",g,"Inf in z w2"))
    if (sum(z.old==-Inf)>0) print(paste("Iteration",g,"-Inf in z w2"))
    #z.old[s,] <- rtnorm(N,mean=cond.mu,sd=sqrt(cond.R),lower=bounds[s,,1],upper=bounds[s,,2])
  }
  return(z.old)
}

# draw R given z and beta:
draw.R.da <- function(S,R,z,x.mat,beta,N){    
  alpha <- rgamma(S,shape=(S+1)/2,rate=1)
  D <- diag(sqrt(diag(solve(R))/(2*alpha)))
  e.star <- D %*% (z - t(x.mat%*%beta))
  Sigma <- solve(my.rwish(N+S+1,solve(e.star %*% t(e.star))))
  R <- cov2cor(Sigma)
  return(R)
}

# function to draw from inverse cdf:
icdf.draw <- function(x,y,d){
  n.grid <- length(x)
  pc <- sum(diff(x)[1]*mean(y))*n.grid
  cdf <- c(0,cumsum(y/pc)*diff(x)[1])
  u <- runif(1)
  upper <- min(which(cdf>u))
  a <- (u-cdf[upper-1])/(cdf[upper]-cdf[upper-1])
  return(a*d[upper] + (1-a)*d[upper-1])
}

# function to draw the correlation parameters from their conditional distributions:
rij.draw <- function(x.mat,beta,l.grid=40,S,R,z){
  mu <- x.mat %*% beta
  l.grid <- 40
  for (i in 1:(S-1)){
    for (j in (i+1):S){
      int <- get.int(R,i,j)
      draw.seq <- seq(int[1],int[2],length=l.grid+1)
      r.seq <- (seq(int[1],int[2],length=l.grid+1) + 0.5*(int[2]-int[1])/(l.grid))[1:l.grid]
      p.grid <- numeric(l.grid)
      for (l in 1:l.grid){
        R.temp <- R; R.temp[i,j] <- r.seq[l]; R.temp[j,i] <- r.seq[l]
        p.grid[l] <- sum(dmvnorm(t(z)-mu,rep(0,S),R.temp,log=TRUE))
      }
      new.r <- icdf.draw(r.seq,exp(p.grid-max(p.grid)),draw.seq)
      R[i,j] <- new.r; R[j,i] <- new.r
    }
  }
  return(R)
}

# function to draw correlations using Metropolis-Hastings
rij.mh <- function(S,R,z,x.mat,beta,jump,adapt=adapt,g=g){
  mu <- x.mat %*% beta
  for (i in 1:(S-1)){
    for (j in (i+1):S){
      ll.data <- sum(dmvnorm(t(z)-mu,rep(0,S),R,log=TRUE))
      lp.old <- ll.data
      r.star <- R[i,j] + rnorm(1,0,jump[i,j])
      int <- get.int(R,i,j)
      if (r.star > int[1] & r.star < int[2]){
        R.temp <- R
        R.temp[i,j] <- r.star; R.temp[j,i] <- r.star
        ll.star <- sum(dmvnorm(t(z)-mu,rep(0,S),R.temp,log=TRUE))
        lp.star <- ll.star
        r <- exp(lp.star - lp.old)
        if (g < adapt) jump[i,j] <- ifelse(r > 0.44,jump[i,j]*A1,jump[i,j]*B1)
        if (runif(1) < r){
          R[i,j] <- r.star
          R[j,i] <- r.star
        }
      } else {
         jump[i,j] <- jump[i,j]*B1
      }
    }
  }
  return(list(R=R,jump=jump))
}


# another way to draw tau:
tau.draw.ig <- function(tau,z,w,alpha,arc,site.mat,T,jump.tau,g,adapt,L,ws,x.arc){
  new.tau <- tau
  for (s in 1:S){
  	yv <- z[,ws[[s]]]-rep(w[,s],L[s])+alpha[s]*x.arc[,ws[[s]]]
    #new.tau[s] <- sqrt(1/rgamma(1,shape=4+0.5*L[s]*T,scale=0.5+0.5*sum((yv-mean(yv))^2)))
    new.tau[s] <- sqrt(1/rgamma(1,shape=4+0.5*L[s]*T,scale=0.5+0.5*sum((yv-0)^2)))
  }
  return(new.tau)
}



























# EOF
