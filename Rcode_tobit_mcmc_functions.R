#2345678901234567890123456789012345678901234567890123456789012345678901234567890

# Simulate Function W.new 
# function to simulate rainfall given parameters from gibbs output:
sim.W <- function(alpha, beta, lambda, tau, beta.arc, Sigma, X, X.arc, 
                  na.preserve = TRUE, na.mat){
  Z.sim <- mvrnorm(N, rep(0, S), tau^2*R.cov(lambda, d.mat)) + X %*% beta
  gamma.mat <- matrix(rgamma(N*S, shape = alpha/2, scale = 2/alpha), N, S)  
  W.new <- as.list(rep(NA, S))
  for (s in 1:S){
    mu.W <- t(Z.sim[, s] + t(matrix(X.arc[[s]]*beta.arc[s], J[s], N)))
    temp <- t(mvrnorm(n = N, mu = rep(0, J[s]), 
                      Sigma = Sigma[[s]])/sqrt(gamma.mat[, s]))
    temp <- temp + mu.W
    temp[temp < 0] <- 0
    W.new[[s]] <- temp
    if (na.preserve) W.new[[s]][na.mat[[s]] == 1] <- NA
  }
  return(W.new)
}


# Traceplot Function  
# trace plot function from myRfunctions.R:
tp <- function(gibbs.input, burn = 0, end = dim(gibbs.input)[2], nc = 3, 
               ylim = range(gibbs.input[, (burn + 1):end]), thin = 1, ...){
  if(nc == 1){
    z <- matrix(gibbs.input, ncol=1)
  } else {
    z <- matrix(t(gibbs.input), ncol = nc)
  }  
  G.local <- dim(z)[1]
  thin.seq <- seq(thin, G.local, by = thin)
  lt <- length(thin.seq)
  plot(thin.seq, z[thin.seq, 1], col = 1, type = "l", ylim = ylim, ..., 
       xaxt = "n")
  axis(1, at = seq(0, G.local, length = 5))
  if(nc > 1){
    for (i in 2:nc) lines(thin.seq, z[thin.seq, i], col = i)
  }
}


# define student t density function:
dst <- function(x, nu, mu, sigma) {
  gamma((nu + 1)/2)/(gamma(nu/2)*sqrt(nu*pi*sigma^2))*(1 + 1/nu*((x - mu)/sigma)^2)^(-1*(nu + 1)/2)
}

# rwishart from package bayesm:
rwishart <- function (nu,  V){
  m = nrow(V)
  df = (nu + nu - m + 1) - (nu - m + 1):nu
  if (m  >  1) {
	T <- diag(sqrt(rchisq(c(rep(1, m)), df)))
	T[lower.tri(T)] = rnorm((m * (m + 1)/2 - m))
  } else {
    T <- sqrt(rchisq(1, df))
  }
  U = chol(V)
  C = t(T) %*% U
  CI = backsolve(C, diag(m))
  return(list(W = crossprod(C), IW = crossprod(t(CI)), C = C, CI = CI))
}

# sample from a truncated normal on one side using an exponential rejection 
# sampling scheme:
rtnorm.exp <- function(n = 1, mean = 0, sd = 1, z.lower = -Inf, right = 1, 
                       iter.max = 10) {
  a.star <- (z.lower + sqrt(z.lower^2 + 4))/2
  ret.val <- rep(NA, n)
  w <- 1:n
  for (i in 1:iter.max){
    e <- rexp(n, rate = a.star) + z.lower
    rho <- exp(-0.5*(e - a.star)^2)
    u <- runif(n) < rho
    ret.val[w][which(u[w])] <- e[w][which(u[w])]
    n.new <- sum(!u[w])
    if (n.new == 0){
      break
    } else {
      w <- which(is.na(ret.val))
    }
  }
  return(ret.val*(-1)^right*sd + mean)
}

# create exponential covariance matrix:
R.cov <- function(lambda, d){
  R <- diag(dim(d)[1])
  R[upper.tri(R)] <- exp(-lambda*d[upper.tri(d)])
  R[lower.tri(R)] <- exp(-lambda*d[lower.tri(d)])
  return(R)
}

# Draw W:
draw.W.start <- function(R, Z, beta.arc, Sigma, gamma.mat, X.arc, up, w.draw, 
                         W) {
  S <- length(R)
  T <- dim(R[[1]])[2]
  J <- unlist(lapply(R, dim))[seq(1, by = 2, length = S)]
  for (s in 1:S){
    # J[s] x T matrix
    mu.mat <- t(matrix(Z[, s], T, J[s]) + 
              t(matrix(X.arc[[s]]*beta.arc[s], J[s], T)))
    for (j in 1:J[s]){ # draw j given -j
      W[[s]][j, ] <- R[[s]][j, ]
      cond.mu <- mu.mat[j, ]
      cond.Sigma <- rep(Sigma[[s]][j, j], T)
      z.score <- cond.mu/as.numeric(sqrt(cond.Sigma))
      sel <- w.draw[[s]][j, ]
      # high z.score and truncated on the right by zero:
      w2 <- which(z.score[sel] > 3.2 & up[[s]][j, sel] == 0)
      w3 <- 1:sum(sel)
      if (length(w2) > 0) w3 <- w3[-w2]
      if (length(w3) > 0) {
        W[[s]][j, sel][w3] <- rtnorm(n = length(w3), 
      	                             mean = cond.mu[sel][w3], 
      	                             sd = sqrt(cond.Sigma[sel][w3]), 
      	                             lower = -Inf, 
      	                             upper = up[[s]][j, sel][w3])
      }
      if (sum(is.na(W[[s]][j, ])) > 0) {
      	print(paste("Iteration", g, "NA in W w3"))
      }
      if (sum(W[[s]][j, ] == Inf) > 0) {
      	print(paste("Iteration", g, "Inf in W w3"))
      }
      if (sum(W[[s]][j, ] == -Inf) > 0) {
      	print(paste("Iteration", g, "-Inf in W w3"))
      }
      if (length(w2) > 0) {
        W[[s]][j, sel][w2] <- rtnorm.exp(n = length(w2), 
      	                                 mean = cond.mu[sel][w2], 
      	                                 sd = sqrt(cond.Sigma[sel][w2]), 
      	                                 z.lower = z.score[sel][w2], 
      	                                 right = 1, 
      	                                 iter.max = 100)
      }
      if (sum(is.na(W[[s]][j, ])) > 0) {
      	print(paste("Iteration", g, "NA in W w2"))
      }
      if (sum(W[[s]][j, ] == Inf) > 0) {
      	print(paste("Iteration", g, "Inf in W w2"))
      }
      if (sum(W[[s]][j, ] == -Inf) > 0) {
      	print(paste("Iteration", g, "-Inf in W w2"))
      }
    }
  }
  return(W)
}

# draw gamma:
draw.gamma <- function(S, Z, X.arc, beta.arc, T, J, Sigma, W, alpha, 
                       gamma.temp) {
  for (s in 1:S){      
    mu.W <- t(Z[, s] + matrix(rep(X.arc[[s]]*beta.arc[s], each = T), 
                              ncol = J[s]))
    S.inv <- solve(Sigma[[s]])
    lhs <- t(W[[s]] - mu.W)
    out <- .C("draw_gamma", 
              params = as.integer(c(J[s], T)), 
              Sinv_in = as.double(S.inv), 
              lhs_in = as.double(lhs), 
              rate_in = as.double(rep(0, T)), 
              alpha_in = as.double(alpha))
    rate.vec <- out$rate_in
    gamma.temp[, s] <- rgamma(T, shape = (J[s] + alpha)/2, rate = rate.vec)
  }
  return(gamma.temp)
}

# Draw Z | W,  gamma,  beta.arc,  Sigma,  beta,  lambda,  tau:
draw.Z <- function(S, Sigma, T, J, X.arc, beta.arc, W, lambda, d.mat, tau, X, 
                   beta, gamma.mat) {
  # Compute mean vector for each draw, from the data:
  sigma.z <- numeric(S)
  for (s in 1:S) sigma.z[s] <- sum(solve(Sigma[[s]]))
  mu.z <- matrix(0, S, T)
  for (s in 1:S) {
  	A <- matrix(-1, 1, J[s])
  	B <- solve(Sigma[[s]])
  	C <- matrix(X.arc[[s]]*beta.arc[s], J[s], T) - W[[s]]
    mu.z[s, ] <- (A %*% B %*% C)/sigma.z[s]
  }
  # Compute prior covariance:
  prior.inv <- solve(R.cov(lambda, d.mat))/tau^2
  # Compute prior mean of multivariate z (from beta):
  Pxb <- prior.inv %*% t(X%*%beta)
  # Compute covariance matrix for each time point (an (S x S x T) array?)
  Sig <- t(gamma.mat)*sigma.z
  # Here,  Sig is the inverse of Sigma^z_t,  i.e. the data precision matrix
  # Multiply Sig by mu.z to get the left addition element in the right hand 
  # side of the mean of z_t:
  Sig.mu <- Sig*mu.z
  z.draw <- .C("zdraw", 
               params = as.integer(c(S, T)), 
               Smu_in = as.double(Sig.mu), 
               Pxb_in = as.double(Pxb), 
               Sig_in = as.double(Sig), 
               Prior_in = as.double(prior.inv))
  Z.out <- t(matrix(z.draw$Smu_in, S, T))
  return(Z.out)
}

# Draw beta.arc | W,  gamma,  Z,  Sigma,  mu.arc,  tau.arc:
draw.beta.arc <- function(S, Sigma, gamma.mat, X.arc, Z, T, J, W, tau.arc, 
                          mu.arc) {
  beta.arc.out <- numeric(S)
  mu.ba <- numeric(S)
  sigma.ba <- numeric(S)
  for (s in 1:S){
    solve.S <- solve(Sigma[[s]])
    denom <- sum(gamma.mat[, s])*t(X.arc[[s]]) %*% solve.S %*% X.arc[[s]]
    tmp <- apply(gamma.mat[, s]*(matrix(Z[, s], T, J[s]) - t(W[[s]])), 2, sum)
    numer <- -t(X.arc[[s]]) %*% solve.S %*% tmp
    mu.w <- numer/denom
    ssq.w <- 1/denom
    d2 <- 1/ssq.w + 1/tau.arc^2
    mu.ba[s] <- (mu.w/ssq.w + mu.arc/tau.arc^2)/d2
    sigma.ba[s] <- sqrt(1/d2)
    beta.arc.out[s] <- rnorm(1, mu.ba[s], sigma.ba[s])
  }
  return(beta.arc.out)
}  


# draw tau.alpha with a metropolis step:
ta.draw <- function(alpha, mu.alpha, tau.alpha, jump.ta, g, adapt) {
  ll.data <- sum(log(dnorm(alpha, mu.alpha, tau.alpha)))
  lp.old <- ll.data + log(dst(tau.alpha, nu = 3, mu = 0, sigma = 0.5))
  tau.alpha.star <- abs(tau.alpha + rnorm(1, 0, jump.ta))
  ll.data.star <- sum(log(dnorm(alpha, mu.alpha, tau.alpha.star)))
  lp.star <- ll.data.star + 
             log(dst(tau.alpha.star, nu = 3, mu = 0, sigma = 0.5))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.ta <- ifelse(r > 0.44, jump.ta*A1, jump.ta*B1)
  if (runif(1) < r) tau.alpha <- tau.alpha.star
  return(list(tau.alpha = tau.alpha, jump.ta = jump.ta))
}

# Sample Sigma | W,  mu (Z,  X.arc,  beta.arc),  gamma:
draw.Sigma <- function(S, Z, X.arc, beta.arc, T, J, W, gamma.mat, v.0, 
                       Lambda.0.inv, Sigma.null) {
  for (s in 1:S){
    mu.W <- t(Z[, s] + 
              matrix(rep(X.arc[[s]]*beta.arc[s], each = T), ncol = J[s]))
    SS <- (t(t(W[[s]])*gamma.mat[, s]) - t(t(mu.W)*gamma.mat[, s])) %*% 
          t(W[[s]] - mu.W)
  	Sigma.null[[s]] <- rwishart(nu = v.0[s] + T, 
  	                            V = solve(Lambda.0.inv[[s]] + SS))$IW
  }
  return(Sigma.null)
}

# draw beta given z, R:
draw.beta.tobit <- function(tau, lambda, vec, xtx.ones, Z, X, sigma, mu) {
  Covar <- tau^2*R.cov(lambda, d.mat)
  S <- dim(Covar)[1]
  P <- dim(beta)[1]
  Covar.inv <- solve(Covar)
  B <- matrix(Covar.inv[vec]*as.numeric(xtx.ones), S*P, S*P)
  W.temp <- matrix(0, P, S)   
  for (j in 1:P){
    for (s in 1:S){
      temp <- 0
      for (l in 1:S) temp <- temp + Covar.inv[s, l]*sum(Z[, l]*X[, j])
      W.temp[j, s] <- W.temp[j, s] + temp
    }
  }
  W.temp <- matrix(W.temp, ncol = 1)
  sig.prior <- diag(1/rep(sigma^2, S))
  B <- B + sig.prior
  B.inv <- solve(B)
  beta.prior <- rep(mu, S)
  beta.hat <- B.inv %*% (W.temp + sig.prior %*% beta.prior)    
  beta <- matrix(mvrnorm(1, beta.hat, B.inv), P, S)
  return(list(beta = beta, B = B))
}

tau.draw.tobit <- function(Z, beta, lambda, a, b, T) {
  V <- solve(R.cov(lambda, d.mat))
  b.post <- b + 0.5*sum(((Z - X%*%beta) %*% V)*(Z - X%*%beta))
  a.post <- a + T*S*0.5
  return(sqrt(1/rgamma(1, shape = a.post, rate = b.post)))
}


# draw lambda with a metropolis step:
lambda.draw.tobit <- function(X, beta, tau, lambda, d.mat, Z, S, jump.lambda) {
  mu.Z <- X %*% beta
  Covar <- tau^2*R.cov(lambda, d.mat)
  ll.data <- sum(dmvnorm(Z-mu.Z, rep(0, S), Covar, log = TRUE))
  lp.old <- ll.data + log(dgamma(lambda, shape = 50, scale = 0.03))
  lambda.star <- abs(lambda + rnorm(1, 0, jump.lambda))
  Covar.star <- tau^2*R.cov(lambda.star, d.mat)
  ll.data.star <- sum(dmvnorm(Z-mu.Z, rep(0, S), Covar.star, log = TRUE)) 
  lp.star <- ll.data.star + log(dgamma(lambda.star, shape = 50, scale = 0.03))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.lambda <- ifelse(r > 0.44, jump.lambda*A1, jump.lambda*B1)
  if (runif(1) < r) lambda <- lambda.star
  return(list(lambda = lambda, jump.lambda = jump.lambda))
}


# script for metropolis-hastings draw of mean and sds of beta:
mu.sigma.draw <- function(beta, mu, sigma, S, P, jump.sigma, g, adapt) {
  ll.data <- apply(matrix(dnorm(beta, mu, sigma, log = TRUE), ncol = S), 1, sum)
  lp.old <- ll.data + log(dst(sigma, nu = 3, mu = 0, sigma = 1))
  sigma.star <- abs(sigma + rnorm(P, 0, jump.sigma))
  ll.data.star <- apply(matrix(dnorm(beta, mu, sigma.star, log = TRUE), 
                               ncol = S), 1, sum)
  lp.star <- ll.data.star + log(dst(sigma.star, nu = 3, mu = 0, sigma = 1))
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.sigma <- ifelse(r > 0.44, jump.sigma*A1, jump.sigma*B1)
  promote <- runif(P) < r
  sigma[promote] <- sigma.star[promote]
  mu <- rnorm(P, apply(beta, 1, mean), sigma/sqrt(S))
  return(list(mu = mu, sigma = sigma, jump.sigma = jump.sigma))
}

# Draw W, the multivariate student-t random variable:
draw.W <- function(R, Z, beta.arc, Sigma, gamma.mat, X.arc, W, up, w.draw) {
  S <- length(R)
  T <- dim(R[[1]])[2]
  J <- unlist(lapply(R, dim))[seq(1, by = 2, length = S)]
  for (s in 1:S){
    # J[s] x T matrix
    mu.mat <- t(matrix(Z[, s], T, J[s]) + 
              matrix(rep(X.arc[[s]]*beta.arc[s], each = T), ncol = J[s]))
    for (j in 1:J[s]){ # draw j given -j
      Sigma.temp <- Sigma[[s]][j, -j] %*% solve(Sigma[[s]][-j, -j])
      cond.mu <- mu.mat[j, ] + Sigma.temp %*% (W[[s]][-j, ] - mu.mat[-j, ])
      cond.numerator <- (Sigma[[s]][j, j] - Sigma.temp %*% Sigma[[s]][-j, j])
      cond.Sigma <- cond.numerator/gamma.mat[, s]
      z.score <- cond.mu/as.numeric(sqrt(cond.Sigma))
      sel <- w.draw[[s]][j, ]
      # high z.score and truncated on the right by zero:
      w2 <- which(z.score[sel] > 3.2 & up[[s]][j, sel] == 0)
      w3 <- 1:sum(sel)
      if (length(w2) > 0) w3 <- w3[-w2]
      if (length(w3) > 0) {
      	W[[s]][j, sel][w3] <- rtnorm(n = length(w3), 
      	                             mean = cond.mu[sel][w3], 
      	                             sd = sqrt(cond.Sigma[sel][w3]), 
      	                             lower = -Inf, 
      	                             upper = up[[s]][j, sel][w3])
      }
      if (sum(is.na(W[[s]][j, ])) > 0) {
      	print(paste("Iteration", g, "NA in W w3"))
      }
      if (sum(W[[s]][j, ] == Inf) > 0) {
      	print(paste("Iteration", g, "Inf in W w3"))
      }
      if (sum(W[[s]][j, ] == -Inf) > 0) {
      	print(paste("Iteration", g, "-Inf in W w3"))
      }
      if (length(w2) > 0) {
      	W[[s]][j, sel][w2] <- rtnorm.exp(n = length(w2), 
      	                                 mean = cond.mu[sel][w2], 
      	                                 sd = sqrt(cond.Sigma[sel][w2]), 
      	                                 z.lower = z.score[sel][w2], 
      	                                 right = 1, 
      	                                 iter.max = 100)
      }
      if (sum(is.na(W[[s]][j, ])) > 0) {
      	print(paste("Iteration", g, "NA in W w2"))
      }
      if (sum(W[[s]][j, ] == Inf) > 0) {
      	print(paste("Iteration", g, "Inf in W w2"))
      }
      if (sum(W[[s]][j, ] == -Inf) > 0) {
      	print(paste("Iteration", g, "-Inf in W w2"))
      }
    }
  }
  return(W)
}



##### Not currently used #####

# draw alpha with a metropolis step:
alpha.draw <- function(gamma.mat, a, b, jump.alpha){
  ll.data <- sum(dgamma(gamma.mat, shape = alpha/2, scale = 2/alpha, 
                        log = TRUE))
  lp.old <- ll.data + dgamma(alpha, shape = a, scale = b, log = TRUE)
  alpha.star <- abs(alpha + rnorm(1, 0, jump.alpha))
  ll.data.star <- sum(dgamma(gamma.mat, shape = alpha.star/2, 
                             scale = 2/alpha.star, log = TRUE))
  lp.star <- ll.data.star + dgamma(alpha.star, shape = a, scale = b, log = TRUE)
  r <- exp(lp.star - lp.old)
  if (g < adapt) jump.alpha <- ifelse(r > 0.44, jump.alpha*A1, jump.alpha*B1)
  if (runif(1) < r) alpha <- alpha.star
  return(list(alpha = alpha, jump.alpha = jump.alpha))
}

# draw alpha with a metropolis step:
alpha.draw.unif <- function(gamma.mat, a, b, jump.alpha){
  ll.data <- sum(dgamma(gamma.mat, shape = alpha/2, scale = 2/alpha, 
                        log = TRUE))
  lp.old <- ll.data + dunif(alpha, a, b, log = TRUE)
  alpha.star <- alpha + rnorm(1, 0, jump.alpha)
  if (alpha.star > a & alpha.star < b) {
    ll.data.star <- sum(dgamma(gamma.mat, shape = alpha.star/2, 
                               scale = 2/alpha.star, log = TRUE))
    lp.star <- ll.data.star + dunif(alpha, a, b, log = TRUE)
    r <- exp(lp.star - lp.old)
  } else {
    r <- 0
  }
  if (g < adapt) jump.alpha <- ifelse(r > 0.44, jump.alpha*A1, jump.alpha*B1)
  if (runif(1)  <  r) alpha <- alpha.star
  return(list(alpha = alpha, jump.alpha = jump.alpha))
}

