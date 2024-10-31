##===============================================
## Bayesian Nonparametric LCA Using NIMBLE
## Chinese Restaurant Process (CRP) 
##===============================================

##===============================================
## Required libraries & self-defined functions
##===============================================

# Install packages if required
if (!require("nimble")) install.packages("nimble") 
if (!require("label.switching")) install.packages("label.switching") 
if (!require("devtools")) install.packages("devtools") 
library(devtools); install_github("sarawade/mcclust.ext") 

# Load libraries 
lab <- c('nimble', 'mcclust.ext', 'label.switching', 
         'poLCA', 'MASS', 'abind', 'mcmcr') 
lapply(lab, require, character.only = TRUE) 

# Function used to simulate binary data with values 0 and 1
lcadat.sim <- function(J, K, N, props, ips) { 
  trueclass <- sample(1:K, prob=props, size=N, replace=TRUE) 
  dat <- matrix(0, nrow=N, ncol=J)
  for(i in 1:N) {
    dat[i,] <- rbinom(n=J, size=1, p=ips[trueclass[i],])
  }
  dat <- as.data.frame(dat) 
  colnames(dat) <- paste0('Y',1:J) 
  res <- list(dat = dat, trueclass = trueclass) 
  return(res) 
}

# Function used to sample from Normal-Gamma distribution
rnormgamma <- function(n, mu0, lambda0, alpha0, beta0) { 
  if (length(n) > 1) 
    n <- length(n) 
  tau <- rgamma(n, alpha0, beta0) 
  mu <- rnorm(n, mu0, sqrt(1/(lambda0*tau))) 
  data.frame(mu = mu, tau = tau) 
} 

# Function used to simulate multivariate Gaussian mixture data
mvndat.sim <- function(J, K, N, props, mus, vars) {
  trueclass <- sample(1:K, prob = props, size = N, replace = TRUE)
  dat <- matrix(0, nrow = N, ncol = J)
  for(i in 1:N) {
    dat[i, ] <- mvrnorm(n = 1, mu = mus[[trueclass[i]]], 
                        Sigma = vars[[trueclass[i]]])
  }
  res <- list(dat = dat, trueclass = trueclass)
  return(res)
}

# Function used to simulate mixture of binary and continous variables
# J: number of binary variables
# L: number of continuous variables
mixdat.sim <- function(J, L, K, N, props, ips, mus, vars) {
  trueclass <- sample(1:K, prob = props, size = N, replace = TRUE)
  # simulate binary variables 
  dat.b <- matrix(0, nrow=N, ncol=J)
  for(i in 1:N) {
    dat.b[i,] <- rbinom(n=J, size=1, p=ips[trueclass[i],])
  }
  # simulate continuous variables 
  dat.c <- matrix(0, nrow = N, ncol = L)
  for(i in 1:N) {
    dat.c[i, ] <- mvrnorm(n = 1, mu = mus[[trueclass[i]]], 
                          Sigma = vars[[trueclass[i]]])
  }
  # combine into one dataset
  dat <- cbind(dat.b, dat.c)
  cnames <- c(paste0('B', 1:J), paste0('C', 1:L))
  colnames(dat) <- cnames
  res <- list(dat = dat, trueclass = trueclass)
  return(res)
}


##===================================
## 1. BNP-LCA with binary variables
##===================================

# Simulation conditions
K <- 3   # number of classes in the population
J <- 6   # number of binary items
N <- 500 # sample size
P <- c(0.5, 0.3, 0.2) # class weights
ip <- matrix(c(rep(0.9,J), rep(c(0.9,0.1),c(J/2,J/2)), rep(0.1,J)), # item response probabilities (IRPs)
             nrow = K, ncol = J, byrow = T)

# Simulate a dataset
set.seed(123)
simdat <- lcadat.sim(J, K, N, P, ip)

# DPM-LCA starts here
# Nimble model
bcode <- nimbleCode({
  for(i in 1:N) {
    for(j in 1:J) {
      y[i, j] ~ dbern(ip[z[i], j])
    }
  }
  # specify beta prior to IRPs
  for(i in 1:N) {
    for(j in 1:J) {
      ip[i, j] ~ dbeta(shape1 = a1, shape2 = a2)
    }
  }
  # CRP construction
  alpha ~ dgamma(shape = 2, rate = 2)
  z[1:N] ~ dCRP(conc = alpha, size = N)
})

# Run the MCMC
nMCMC <- 10000
constants <- list(N = N, J = J, a1 = 1, a2 = 1)
set.seed(123)
inits <- list(ip = matrix(rbeta(N*J, shape1 = 1, shape2 = 1), nrow = N, ncol = J),
              z = sample(1:5, size = constants$N, replace = TRUE))
dat.b <- list(y = as.matrix(simdat$dat))
model <- nimbleModel(bcode, data = dat.b, inits = inits, constants = constants,
                     dimensions = list(ip = c(N, J))) 
cmodel <- compileNimble(model) 
conf <- configureMCMC(model, monitors = c('z','alpha','ip'), print = TRUE) 
mcmc <- buildMCMC(conf) 
cmcmc <- compileNimble(mcmc, project = model)
samples.b <- runMCMC(cmcmc,  niter = nMCMC, nburnin = 0, thin = 1, setSeed = TRUE) # no burning and thinning
zSamples <- samples.b[, grep('z', colnames(samples.b))]
aSamples <- samples.b[, grep('alpha', colnames(samples.b))]
ipSamples <- samples.b[, grep('ip', colnames(samples.b))]
nGroups <- apply(zSamples, 1, function(x) length(unique(x)))

# Perform burning and thinning based on the convergence and auto-correlation plots
nburn <- 1000 
nthin <- 2 
z.post <- zSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 
a.post <- aSamples[nburn+seq(1, (nMCMC-nburn), nthin)] 
ip.post <- ipSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 

# Align the within-cluster parameters (i.e., ip) with the cluster weights 
nMCMC.post <- (nMCMC-nburn)/nthin 
nGroups.post <- nGroups[nburn+seq(1, (nMCMC-nburn), nthin)] 
w.list <- ip.list <- list() 
for(iter in 1:nMCMC.post) { 
  w.list[[iter]] <- as.numeric(table(z.post[iter,]))/N
  ids <- sort(unique(z.post[iter, ]))
  ips <- matrix(nrow = J, ncol = nGroups.post[iter])
  ip.m <- matrix(ip.post[iter,], nrow = J, ncol = N, byrow = T)
  ip.list[[iter]] <- sapply(ids, function(i) ip.m[, i])
}

# Fill the gap
z.nogap <- matrix(0, nrow = nMCMC.post, ncol = N)
for (iter in 1:nMCMC.post) {
  ids <- sort(unique(as.numeric(z.post[iter, ])))
  z.nogap[iter, ] <- match(z.post[iter, ], ids)
}

# Find the "best" partition using VI
psm <- comp.psm(z.nogap) 
x11() 
plotpsm(psm) # heat map of posterior similarity matrix 
out.VI <- minVI(psm, z.nogap, method = ('all'), include.greedy = TRUE) 
K.est <- length(unique(out.VI$cl[1,])) 
table(out.VI$cl[1,])/N 

# Subset the draws with the estimated number of clusters
nClusters <- apply(z.nogap, 1, function(x) length(unique(x)))
indx <- which(nClusters == K.est)
ndraw <- length(indx) # number of draws filtered
z.filter <- z.nogap[indx,]
w.filter <- w.list[indx]
ip.filter <- ip.list[indx]
w.matrix <- matrix(unlist(w.filter), nrow = ndraw, ncol = K.est, byrow = T)
ip.array <- abind(ip.filter, along = 3)
ip.array <- aperm(ip.array, c(3,2,1)) 

# Relabeling across MCMCs
run <- ecr.iterative.1(z = z.filter, K = K.est)
neworder <- run$permutations
w.reorder <- matrix(0, nrow = ndraw, ncol = K.est)
for(i in 1:ndraw) {
  w.reorder[i,] <- w.matrix[i,][neworder[i,]]
}
colnames(w.reorder) <- paste0('Class',1:K.est)
ip.reorder <- permute.mcmc(ip.array, neworder)

# Summarize the posteriors
P.postmean  <-  apply(w.reorder, 2, mean) # posterior means of proportions
samplesSummary(as.mcmc(w.reorder))  
ip.postmean <- apply(ip.reorder$output, c(2,3), mean)
dimnames(ip.postmean) <- list(paste0('Class',1:K.est), paste0('Item',1:J))
pars.postmean <- list(Mix.weights = P.postmean, IRPs = t(ip.postmean)) 
pars.postmean 


##========================================
## 2. BNP-LCA with polytomous variables
##========================================

# Simulation conditions
J <- 5   # number of items
N <- 500 # sample size
P <- c(0.5,0.3,0.2) # class weights
# 3 classes, 5 polytomous variables with different number of categories
ip <- list(matrix(c(0.6,0.1,0.3,     0.6,0.3,0.1,     0.3,0.1,0.6    ),ncol=3,byrow=TRUE), # Y1
           matrix(c(0.1,0.9,         0.9,0.1,         0.1,0.9        ),ncol=2,byrow=TRUE), # Y2
           matrix(c(0.2,0.7,0.1,     0.1,0.2,0.7,     0.2,0.7,0.1    ),ncol=3,byrow=TRUE), # Y3
           matrix(c(0.1,0.1,0.5,0.3, 0.5,0.3,0.1,0.1, 0.3,0.1,0.1,0.5),ncol=4,byrow=TRUE), # Y4
           matrix(c(0.1,0.1,0.8,     0.1,0.8,0.1,     0.8,0.1,0.1    ),ncol=3,byrow=TRUE)) # Y5

# Simulate a dataset
set.seed(123)
simdat <- poLCA.simdata(N = N, probs = ip, P = P)
C <- apply(simdat$dat, 2, max) # number of categories for each item
C.max <- max(C)

# DPM-LCA starts here
# Nimble model
pcode <- nimbleCode({
  for(i in 1:N) {
    for(j in 1:J) {
      y[i, j] ~ dcat(ip[z[i], j, 1:C[j]])
    }
  }
  # specify Dirichlet prior to IRPs
  for(i in 1:N) {
    for(j in 1:J) {
      ip[i, j, 1:C[j]] ~ ddirch(beta[1:C[j]])
    }
  }
  # CRP construction
  alpha ~ dgamma(shape = 2, rate = 2)
  z[1:N] ~ dCRP(conc = alpha, size = N)
})

# Run the MCMC
nMCMC <- 10000
constants <- list(N = N, J = J, C = C, beta = rep(1, C.max))
set.seed(123)
ip <- array(0, dim = c(N, J, C.max)) # initial values for ip
for(i in 1:N) {
  for(j in 1:J) {
    pp <- MCMCpack::rdirichlet(1, alpha = rep(1,C[j]))
    if(C[j] != C.max) {
      ip[i,j,] <- c(pp, rep(0, C.max-C[j]))
    } else {
      ip[i,j,] <- pp
    }
  }
}
inits <- list(ip = ip, z = sample(1:5, size = constants$N, replace = TRUE)) 
dat.p <- list(y = as.matrix(simdat$dat)) 
model <- nimbleModel(pcode, data = dat.p, inits = inits, constants = constants) 
cmodel <- compileNimble(model) 
conf <- configureMCMC(model, monitors = c('z', 'alpha', 'ip'), print = TRUE) 
mcmc <- buildMCMC(conf) 
cmcmc <- compileNimble(mcmc, project = model)
samples.p <- runMCMC(cmcmc,  niter = nMCMC, nburnin = 0, thin = 1, setSeed = TRUE)
zSamples <- samples.p[, grep('z', colnames(samples.p))]
aSamples <- samples.p[, grep('alpha', colnames(samples.p))]
ipSamples <- samples.p[, grep('ip', colnames(samples.p))]
nGroups <- apply(zSamples, 1, function(x) length(unique(x)))

# Perform burning and thinning based on the convergence and auto-correlation plots
nburn <- 1000 
nthin <- 2 
z.post <- zSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 
a.post <- aSamples[nburn+seq(1, (nMCMC-nburn), nthin)] 
ip.post <- ipSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 

# Align the within-cluster parameters (i.e., ip) with the cluster weights 
nMCMC.post <- (nMCMC-nburn)/nthin 
nGroups.post <- nGroups[nburn+seq(1, (nMCMC-nburn), nthin)] 
w.list <- ip.list <- list() 
for(iter in 1:nMCMC.post) { 
  w.list[[iter]] <- as.numeric(table(z.post[iter,]))/N
  ids <- sort(unique(z.post[iter, ]))
  ips <- array(0, dim=c(J, C.max, nGroups.post[iter]))
  ip.array <- array(ip.post[iter,], dim=c(N, J, C.max))
  ip.list[[iter]] <- sapply(ids, function(i) ip.array[i,,])
}

# Fill the gap
z.nogap <- matrix(0, nrow = nMCMC.post, ncol = N)
for (iter in 1:nMCMC.post) {
  ids <- sort(unique(as.numeric(z.post[iter, ])))
  z.nogap[iter, ] <- match(z.post[iter, ], ids)
}

# Find the "best" partition using VI
psm <- comp.psm(z.nogap)
x11() 
plotpsm(psm) # heat map of posterior similarity matrix
out.VI <- minVI(psm, z.nogap, method = ('all'), include.greedy = TRUE)
K.est <- length(unique(out.VI$cl[1,]))
table(out.VI$cl[1,])/N

# Subset the draws with the estimated number of clusters
nClusters <- apply(z.nogap, 1, function(x) length(unique(x)))
indx <- which(nClusters == K.est)
ndraw <- length(indx) # number of draws filtered
z.filter <- z.nogap[indx, ]
w.filter <- w.list[indx]
ip.filter <- ip.list[indx]
w.matrix <- matrix(unlist(w.filter), nrow = ndraw, ncol = K.est, byrow = T)
ip.array <- abind(ip.filter, along = 3)
ip.array <- aperm(ip.array, c(3,2,1))

# Relabeling across MCMCs
run <- ecr.iterative.1(z = z.filter, K = K.est)
neworder <- run$permutations
w.reorder <- matrix(0, nrow = ndraw, ncol = K.est)
for(i in 1:ndraw) {
  w.reorder[i,] <- w.matrix[i,][neworder[i,]]
}

# Summarize the posteriors
P.postmean <- apply(w.reorder, 2, mean) # posterior means of proportions
names(P.postmean) <- paste0('Class',1:K.est)
ip.reorder <- list()
for(r in 1:ndraw) {
  ip.reorder[[r]] <- ip.array[r,neworder[r,],] # reorder within-cluster parameters
}
ip.temp <- abind(ip.reorder, along = 3)        # combine into a 3D array
ip.temp <- apply(ip.temp, c(1,2), mean, na.rm = T) # take mean across the 3-rd dim
ip.postmean <- lapply(seq_len(nrow(ip.temp)),      # reorganize the parameters
                      function(i) matrix(ip.temp[i,], nrow = J, ncol = C.max, byrow = F))
ip.postmean <- abind(ip.postmean, along = 3)
dimnames(ip.postmean) <- list(paste0('Item',1:J), 
                              paste0('Level',1:C.max), 
                              paste0('Class',1:K.est))
pars.postmean <- list(Mix.weights = P.postmean, IRPs = ip.postmean)
pars.postmean # Note. The zeros indicate empty categories of an item.


##========================================
## 3. BNP-LCA with continuous variables
##========================================

# Simulation conditions
K <- 3   # number of classes 
J <- 8   # number of items
N <- 500 # sample size
P <- c(0.5, 0.3, 0.2) # class weights
mu1 <- rep(1, J); mu2 <- rep(0, J); mu3 <- rep(-1, J)
mus <- list(mu1, mu2, mu3) # item means
Sig <- diag(1, nrow = J)
vars <- replicate(K, Sig, simplify = FALSE) # class-specific covariance matrix

# Simulate a dataset 
set.seed(123)
simdat <- mvndat.sim(J, K, N, props=P, mus, vars)

# DPM-LCA starts here
# Nimble model
ccode <- nimbleCode({
  for(i in 1:N) {
    for(j in 1:J) {
      y[i, j] ~ dnorm(mean = mu[z[i], j], tau = tau[z[i], j])
    }
  }
  for(i in 1:N) {
    for(j in 1:J) {
      # set independent priors to item mean and precision
      mu[i, j] ~ dnorm(mean = nu1, tau = nu2)
      tau[i, j] ~ dgamma(a, b) # tau = 1/variance
    }
  }
  # CRP construction
  alpha ~ dgamma(shape = 2, rate = 2) 
  z[1:N] ~ dCRP(conc = alpha, size = N) 
})

# Run the MCMC
nMCMC <- 10000
constants <- list(N = N, J = J, a = 1, b = 1, nu1 = 0, nu2 = 1/1000)
set.seed(1234) 
mutau.int <- rnormgamma(J*N, 0,1, 1,1) # generate initial values for mu and tau
inits <- list(z = sample(1:5, size = N, replace = TRUE), alpha = 1,      
              mu = matrix(mutau.int$mu, nrow = N, ncol = J), 
              tau = matrix(mutau.int$tau, nrow = N, ncol = J))
dat.c <- list(y = simdat$dat)
model <- nimbleModel(ccode, data = dat.c, inits = inits, constants = constants, 
                     dimensions = list(mu = c(N, J), tau = c(N, J)))
cmodel <- compileNimble(model)
conf <- configureMCMC(model, monitors = c('z','alpha','mu','tau'), print = TRUE)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)
samples.c <- runMCMC(cmcmc, niter = nMCMC, nburnin = 0, thin = 1, setSeed = TRUE)
zSamples <- samples.c[, grep('z', colnames(samples.c))]  
aSamples <- samples.c[, grep('alpha', colnames(samples.c))]   
muSamples <- samples.c[, grep('mu', colnames(samples.c))]
tauSamples <- samples.c[, grep('tau', colnames(samples.c))]
nGroups <- apply(zSamples, 1, function(x)  length(unique(x)))

# Perform burning and thinning based on the convergence and auto-correlation plots
nburn <- 100
nthin <- 2
z.post <- zSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 
a.post <- aSamples[nburn+seq(1, (nMCMC-nburn), nthin)] 
mu.post <- muSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 
tau.post <- tauSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 

# Align the within-cluster parameters with the cluster weights 
nMCMC.post <- (nMCMC-nburn)/nthin 
nGroups.post <- nGroups[nburn+seq(1, (nMCMC-nburn), nthin)] 
w.list <- mu.list <- tau.list <- list()
for(iter in 1:nMCMC.post) {
  w.list[[iter]] <- as.numeric(table(z.post[iter,]))/N
  ids <- sort(unique(z.post[iter, ]))
  mus <- taus <- matrix(nrow = J, ncol = nGroups.post[iter])
  mu.m <- matrix(mu.post[iter,], nrow = J, ncol = N, byrow = T)
  tau.m <- matrix(tau.post[iter,], nrow = J, ncol = N, byrow = T)
  mu.list[[iter]] <- sapply(ids, function(i) mu.m[, i])
  tau.list[[iter]] <- sapply(ids, function(i) tau.m[, i])
}

# Fill the gap
z.nogap <- matrix(0, nrow = nMCMC.post, ncol = N)
for (iter in 1:nMCMC.post) {
  ids <- sort(unique(as.numeric(z.post[iter, ])))
  z.nogap[iter, ] <- match(z.post[iter, ], ids)
}

# Find the "best" partition using VI
psm <- comp.psm(z.nogap) 
x11() 
plotpsm(psm) # heat map of posterior similarity matrix 
out.VI <- minVI(psm, z.nogap, method = ('all'), include.greedy = TRUE) 
K.est <- length(unique(out.VI$cl[1,])) # number of clusters determined by VI 
table(out.VI$cl[1,])/N 

# Subset the draws with the estimated number of clusters
nClusters <- apply(z.nogap, 1, function(x) length(unique(x)))
indx <- which(nClusters == K.est)
ndraw <- length(indx) # number of draws filtered
z.filter <- z.nogap[indx,]
w.filter <- w.list[indx]
mu.filter <- mu.list[indx]
tau.filter <- tau.list[indx]
w.matrix <- matrix(unlist(w.filter), nrow = ndraw, ncol = K.est, byrow = T)
mu.array <- abind(mu.filter, along = 3)
mu.array <- aperm(mu.array, c(3,2,1))
tau.array <- abind(tau.filter, along = 3)
tau.array <- aperm(tau.array, c(3,2,1))

# Relabeling across MCMCs
run <- ecr.iterative.1(z = z.filter, K = K.est)
neworder <- run$permutations
w.reorder <- matrix(0, nrow = ndraw, ncol = K.est)
for(i in 1:ndraw) {
  w.reorder[i,] <- w.matrix[i,][neworder[i,]]
}

# Summarize the posteriors
P.postmean <- apply(w.reorder, 2, mean) 
names(P.postmean) <- paste0('Class',1:K.est)
dim.names <- list(paste0('Class',1:K.est), paste0('Item',1:J))
mu.reorder <- permute.mcmc(mu.array, neworder)
mu.postmean <- apply(mu.reorder$output, c(2,3), mean)
tau.reorder <- permute.mcmc(tau.array, neworder)
var.postmean <- 1/apply(tau.reorder$output, c(2,3), mean) # convert precision to variance
dimnames(mu.postmean) <- dimnames(var.postmean) <- dim.names
pars.postmean <- list(Mix.weights = P.postmean, Means = t(mu.postmean), Vars = t(var.postmean))
pars.postmean 


##=================================================================
## 4. BNP-LCA with mixed types of variables: binary & continuous 
##=================================================================

# Simulation conditions
K <- 3   # number of classes 
J <- 4   # number of binary variables
L <- 4   # number of continuous variables
N <- 500 # sample size
P <- c(0.5, 0.3, 0.2) # class weights
ip <- matrix(c(rep(0.9,J), rep(c(0.9,0.1),c(J/2,J/2)), rep(0.1,J)), # IRPs
             nrow = K, ncol = J, byrow = T)
mu1 <- rep(1, L); mu2 <- rep(0, L); mu3 <- rep(-1, L)
mus <- list(mu1, mu2, mu3) # item means
Sig <- diag(1, nrow = L)
vars <- replicate(K, Sig, simplify = FALSE) # class-specific covariance matrix

# Simulate a dataset
set.seed(1234)
simdat <- mixdat.sim(J, L, K, N, props = P.U, ips = ip, mus, vars)

# DPM-LCA starts here
# Nimble model: CRP construction
mixcode <- nimbleCode ({
  # specify priors for responses
  for (i in 1:N){
    # a Bernoulli distribution for binary response
    for (j in 1:J) {
      y[i, j] ~ dbern(ip[z[i], j])
    }
    # a normal distribution for continuous response
    for (l in 1:L) {
      x[i, l] ~ dnorm(mu[z[i], l], tau = tau[z[i], l])
    }
  } 
  # specify priors for parameters
  for (i in 1:N) {
    for (j in 1:J) {
      # set a beta prior B(a1, a2)
      ip[i, j] ~ dbeta(shape1 = a1, shape2 = a2)
    }
    for (l in 1:L) {
      # set independent priors to item mean and precision
      mu[i, l] ~ dnorm(nu1, tau = nu2)
      tau[i, l] ~ dgamma(b1, b2) # tau = 1/variance
    }
  }
  alpha ~ dgamma(shape = 2, rate = 2) 
  # CRP construction of DP
  z[1:N] ~ dCRP(conc = alpha, size = N) 
})

# Run the MCMC
nMCMC <- 10000
consts <- list(N = N, J = J, L = L, a1 = 1, a2 = 1, b1 = 1, b2 = 1, nu1 = 0, nu2 = 1/1000)
set.seed(123)
mutau.int <- rnormgamma(N*L, 0,2, 2,1) # generate initial values for mu and tau
inits <- list(z = 1:N, alpha = 1, 
              ip = matrix(rbeta(N*J, shape1=1, shape2=1), nrow = N, ncol = J),
              mu = matrix(mutau.int$mu, nrow = N, ncol = L), 
              tau = matrix(mutau.int$tau, nrow = N, ncol = L))
dat.mix <- list(y = simdat$dat[, 1:J], x = simdat$dat[, (J+1):(J+L)])
model <- nimbleModel(mixcode, data = dat.mix, inits = inits, constants = consts)
cmodel <- compileNimble(model)
conf <- configureMCMC(model, monitors = c('z','alpha','ip','mu','tau'), print = TRUE)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)
samples.crp <- runMCMC(cmcmc, niter = nMCMC, nburnin = 0, thin = 1, setSeed = TRUE)
zSamples <- samples.crp[, grep('z', colnames(samples.crp))]  
aSamples <- samples.crp[, grep('alpha', colnames(samples.crp))]  
ipSamples <- samples.crp[, grep('ip', colnames(samples.crp))]
muSamples <- samples.crp[, grep('mu', colnames(samples.crp))]
tauSamples <- samples.crp[, grep('tau', colnames(samples.crp))]
nGroups <- apply(zSamples, 1, function(x)  length(unique(x)))

# Perform burning and thinning based on the convergence and auto-correlation plots
nburn <- 1000 
nthin <- 2
z.post <- zSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 
a.post <- aSamples[nburn+seq(1, (nMCMC-nburn), nthin)] 
ip.post <- ipSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 
mu.post <- muSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 
tau.post <- tauSamples[nburn+seq(1, (nMCMC-nburn), nthin), ] 

# Align the within-cluster parameters with the cluster weights 
nMCMC.post <- (nMCMC-nburn)/nthin 
nGroups.post <- nGroups[nburn+seq(1, (nMCMC-nburn), nthin)] 
w.list <- ip.list <- mu.list <- tau.list <- list() 
for(iter in 1:nMCMC.post) { 
  w.list[[iter]] <- as.numeric(table(z.post[iter,]))/N 
  ids <- sort(unique(z.post[iter, ])) 
  ips <- matrix(nrow = J, ncol = nGroups.post[iter]) 
  ip.m <- matrix(ip.post[iter,], nrow = J, ncol = N, byrow = T) 
  mus <- taus <- matrix(nrow = L, ncol = nGroups.post[iter]) 
  mu.m <- matrix(mu.post[iter,], nrow = L, ncol = N, byrow = T) 
  tau.m <- matrix(tau.post[iter,], nrow = L, ncol = N, byrow = T) 
  ip.list[[iter]] <- sapply(ids, function(i) ip.m[, i]) 
  mu.list[[iter]] <- sapply(ids, function(i) mu.m[, i]) 
  tau.list[[iter]] <- sapply(ids, function(i) tau.m[, i]) 
} 

# Fill the gap
z.nogap <- matrix(0, nrow = nMCMC.post, ncol = N)
for (iter in 1:nMCMC.post) {
  ids <- sort(unique(as.numeric(z.post[iter, ])))
  z.nogap[iter, ] <- match(z.post[iter, ], ids)
} 

# Find the "best" partition using VI
psm <- comp.psm(z.nogap) 
x11() 
plotpsm(psm) # heat map of posterior similarity matrix
out.VI <- minVI(psm, z.nogap, method = ('all'), include.greedy = TRUE) 
K.est <- length(unique(out.VI$cl[1,])) 
table(out.VI$cl[1,])/N 

# Subset the draws with the estimated number of clusters
nClusters <- apply(z.nogap, 1, function(x) length(unique(x))) 
indx <- which(nClusters == K.est) 
ndraw <- length(indx) # number of draws filtered
z.filter <- z.nogap[indx,] 
w.filter <- w.list[indx] 
ip.filter <- ip.list[indx] 
mu.filter <- mu.list[indx] 
tau.filter <- tau.list[indx] 
w.matrix <- matrix(unlist(w.filter), nrow = ndraw, ncol = K.est, byrow = T) 
ip.array <- abind(ip.filter, along = 3) 
ip.array <- aperm(ip.array, c(3,2,1)) 
mu.array <- abind(mu.filter, along = 3) 
mu.array <- aperm(mu.array, c(3,2,1)) 
tau.array <- abind(tau.filter, along = 3) 
tau.array <- aperm(tau.array, c(3,2,1)) 

# Relabeling across MCMC samples
run <- ecr.iterative.1(z = z.filter, K = K.est) 
neworder <- run$permutations 
w.reorder <- matrix(0, nrow = ndraw, ncol = K.est) 
for(i in 1:ndraw) {
  w.reorder[i,] <- w.matrix[i,][neworder[i,]] 
} 
ip.reorder <- permute.mcmc(ip.array, neworder)
mu.reorder <- permute.mcmc(mu.array, neworder)
tau.reorder <- permute.mcmc(tau.array, neworder)

# Summarize the posteriors
P.postmean <- apply(w.reorder, 2, mean) # posterior means of proportions
names(P.postmean) <- paste0("Class",1:K.est) 
dim.names <- list(paste0("Class",1:K.est), paste0("Item",1:J)) 
ip.postmean <- apply(ip.reorder$output, c(2,3), mean) 
mu.postmean <- apply(mu.reorder$output, c(2,3), mean) 
var.postmean <- 1/apply(tau.reorder$output, c(2,3), mean) # convert precision to variance
dimnames(ip.postmean) <- dimnames(mu.postmean) <- dimnames(var.postmean) <- dim.names
pars.postmean <- list(Mix.weights = P.postmean, IRPs = ip.postmean, 
                      Means = t(mu.postmean), Vars = t(var.postmean)) 
pars.postmean 

## End
