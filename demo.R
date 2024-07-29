## load packages and R functions
library(wavelets)
library(doParallel)
library(foreach)
source('fgraph_cclasso.R')
source('Cmat_update.R')
source('Dmat_update.R')
source('Lam_update.R')


#### specify dynamic graphical model for innovation errors
p = 10
G0 <- diag(0,p)
G0[1,2] <- G0[1,3] <- G0[2,3] <- G0[7,8] <- G0[7,9] <- G0[8,9] <- 1 
G1 <- G2 <- G0
G1[4,5] <- G1[5,6] <- G1[4,6] <- 1
G2[7,10] <- G2[8,10] <- G2[9,10] <- 1
G0 <- G0 + t(G0); G1 <- G1 + t(G1); G2 <- G2 + t(G2)

## specify precision matrix of error 
T = 128
tt = 1:T
r.max <- -0.99
# monotonic change of correlation
r1 <- sapply(tt,function(t) min(max(0,t-T/4),T/2)*(r.max/(T/2))) 
# # periodic change of correlation
# r1 <- sapply(tt,function(t) if(t%%(T/2)<=(T/8)) r.max else if(t%%(T/2)<=(T/4)) r.max/(T/8)*(T/4-t%%(T/2)) else if(t%%(T/2)<=(T*3/8)) 0 else r.max/(T/8)*(t%%(T/2)-(T*3/8)))

Sigma.err <- rep(list(diag(p)),T)
for(t in 1:T) {
    Omega.tmp <- diag(c(rep(2,6),rep(3,4)))
    Omega.tmp[which(G1==1)] <- r1[t]
    Omega.tmp[which(G2==1)] <- Omega.tmp[which(G2==1)] + r2[t]
    Sigma.tmp <- solve(Omega.tmp)
    Sigma.err[[t]] <- diag(diag(Sigma.tmp)^-0.5) %*% Sigma.tmp %*% diag(diag(Sigma.tmp)^-0.5)*2
}

# standard deviation of white noise
sig.e <-  sqrt(0.4)


######## Setting 1&2: generate data from AR(1) dynamic model
n <- 50; p=10; T=128
y <- array(NA,c(n,p,T))

for(t in 1:T){
  err <- matrix(rnorm(n*p),n) %*% chol(Sigma.err[[t]]) + matrix(rnorm(n*p),n)*sig.e
  if(t == 1) y[,,t] <- err else
    y[,,t] <- y[,,t-1]*sqrt(0.5) + err*sqrt(0.5)
}


## calculate true conditional cross-covariance graph for AR(1) dynamic model
Sigmas.true <- G.true <- G.neg <- rep(list(diag(p)),T)
for(t in 1:T){
    if(t == 1) {
        Sigmas.true[[t]] <- Sigma <- Sigma.err[[1]] 
    } else {
        Sigmas.true[[t]] <- Sigma <- Sigma * 0.5 + Sigma.err[[t]] * 0.5
    }
    G.true[[t]] <- (abs(solve(Sigma)) > 0.1)*1
    G.neg[[t]] <- (abs(solve(Sigma)) < 0.001)*1
}
G.true <- abind(G.true,along=3)
G.neg <- abind(G.neg,along=3)

    
#### run Bayesian functional graph
  
MCMCspecs = list(B=1000,thin=5,burnin=1000,update=1000);
D_prior = list(a=0.1,b=0.1); lam_prior = list(a=0.1,b=1);
## wavelet basis functions
tmp = diag(T)
W = matrix(NA,T,T)
for(t in 1:T) {
  a = dwt(tmp[t,],filter='haar',n.levels=5)
  W[t,] = c(a@V[[5]],unlist(a@W[5:1]))
}
Phi <- solve(W)

fgbay.samp <- fgraph_ccc(y,Phi,MCMCspecs,D_prior,lam_prior,ncor=1)

fgbay.est <- array(diag(p),c(p,p,T))
for(t in 1:T) fgbay.est[,,t][upper.tri(diag(p))] <- apply(fgbay.samp$ccc[,t,],1,function(x) (quantile(x,0.025)>0 | quantile(x,0.975)<0)*1)


#### calculate TPs and FPs

tmp <- matrix(FALSE,p,p); tmp[upper.tri(tmp)] <- TRUE
uptri.ind <- array(tmp,c(p,p,T))
sum((fgbay.est+G.true==2)[uptri.ind])/sum(G.true[uptri.ind]) 
sum((fgbay.est+G.neg==2)[uptri.ind])/sum(G.neg[uptri.ind])




######## Setting 3&4: generate data from change-point dynamic model
n <- 50; p=10; T=128
y <- array(NA,c(n,p,T))

## generate X1 and X2 from two graphical models A1 and A2
A1 <- A2 <- G0; A1[1,4] <- A1[4,1] <- A2[4,7] <- A2[7,4] <- 1
Omega.A1 <- Omega.A2 <- diag(c(3,2,2,rep(1,3),3,2,2,1)); 
Omega.A1[which(A1==1)] <- Omega.A2[which(A2==1)] <- -0.9;
SS.A1 <- solve(Omega.A1); SS.A2 <- solve(Omega.A2)
Sigma.A1 <- diag(diag(SS.A1)^-0.5) %*% SS.A1 %*% diag(diag(SS.A1)^-0.5)*2
Sigma.A2 <- diag(diag(SS.A2)^-0.5) %*% SS.A2 %*% diag(diag(SS.A2)^-0.5)*2

x1 <- matrix(rnorm(n*p),n) %*% chol(Sigma.A1*0.3)
x2 <- matrix(rnorm(n*p),n) %*% chol(Sigma.A2*0.3)

err <- matrix(rnorm(n*p),n) %*% chol(Sigma.err[[t]]) + matrix(rnorm(n*p),n)*sig.e
for(t in 1:T){
  if(t <= T/2) y[,,t] <- x1 + err*sqrt(0.7) else y[,,t] <- x2 + err*sqrt(0.7)
}


## calculate true conditional cross-covariance graph
Sigmas.true <- G.true <- G.neg <- rep(list(diag(p)),T)
for(t in 1:T){
  if(t <= T/2) {
    Sigmas.true[[t]] <- Sigma <- Sigma.A1*0.3 + Sigma.err[[t]] * 0.7
  } else {
    Sigmas.true[[t]] <- Sigma <- Sigma.A2*0.3 + Sigma.err[[t]] * 0.7
  }
  G.true[[t]] <- (abs(solve(Sigma)) > 0.1)*1
  G.neg[[t]] <- (abs(solve(Sigma)) < 0.001)*1
}

#### run Bayesian functional graph

MCMCspecs = list(B=1000,thin=5,burnin=1000,update=1000);
D_prior = list(a=0.1,b=0.1); lam_prior = list(a=0.1,b=1);
## wavelet basis functions
tmp = diag(T)
W = matrix(NA,T,T)
for(t in 1:T) {
  a = dwt(tmp[t,],filter='haar',n.levels=5)
  W[t,] = c(a@V[[5]],unlist(a@W[5:1]))
}
Phi <- solve(W)

fgbay.samp <- fgraph_ccc(y,Phi,MCMCspecs,D_prior,lam_prior,ncor=1)

fgbay.est <- array(diag(p),c(p,p,T))
for(t in 1:T) fgbay.est[,,t][upper.tri(diag(p))] <- apply(fgbay.samp$ccc[,t,],1,function(x) (quantile(x,0.025)>0 | quantile(x,0.975)<0)*1)


#### calculate TPs and FPs

tmp <- matrix(FALSE,p,p); tmp[upper.tri(tmp)] <- TRUE
uptri.ind <- array(tmp,c(p,p,T))
sum((fgbay.est+G.true==2)[uptri.ind])/sum(G.true[uptri.ind]) 
sum((fgbay.est+G.neg==2)[uptri.ind])/sum(G.neg[uptri.ind])



######## Setting 5&6: dynamic mixture models
## specify five graphical models
p = 10
G0 <- diag(0,p)
G0[1,2] <- G0[1,3] <- G0[2,3] <- G0[7,8] <- G0[7,9] <- G0[8,9] <- 1
G1 <- G2 <- G3 <- G4 <- G0
G1[4,5] <- G1[4,6] <- 1
G2[7,10] <- G2[8,10] <- 1
G3[5,6] <- G3[4,6] <- 1
G4[7,10] <- G4[9,10] <- 1
G.all <- list(G0,G1,G2,G3,G4)

## generate precision matrix for each graph
Sigma.all <- Omega.all <- rep(list(NA),5)
for(m in 1:5){
  Omega.tmp <-  matrix(runif(p*p,-1,-0.8),p,p)*G.all[[m]]
  Omega.tmp <- Omega.tmp + t(Omega.tmp)
  diag(Omega.tmp) <- sapply(rowSums(abs(Omega.tmp))*1.1,function(x) max(x,1))
  Sigma.tmp <- solve(Omega.tmp)
  Sigma.all[[m]] <- diag(diag(Sigma.tmp)^-0.5) %*% Sigma.tmp %*% diag(diag(Sigma.tmp)^-0.5)*2
  Omega.all[[m]] <- solve(Sigma.all[[m]])
}


## generate mixing functions using Fourier basis functions (bounded below or above 0)
library(fda)
T = 128

bfor <- create.fourier.basis(rangeval=c(0,T),nbasis=5)
## Low frequency
phi.tmp <- eval.basis(1:T,bfor)[,c(1,2,3)]
# ## high frequency
# phi.tmp <- eval.basis(1:T,bfor)[,c(1,4,5)]
phi.true <- cbind(phi.tmp*(phi.tmp>0), phi.tmp[,2:3]*(phi.tmp[,2:3]<0))*sqrt(T/3)


## generate data
n <- 50; p=10; T=128
y <- array(NA,c(n,p,T))

ystar <- sapply(1:5,function(m) matrix(rnorm(n*p),n,p) %*% chol(Sigma.all[[m]]),simplify=FALSE)
for(t in 1:T){
  y[,,t] <- Reduce('+',sapply(1:5, function(m) ystar[[m]]*phi.true[t,m], simplify=FALSE)) + matrix(rnorm(n*p),n)*sig.e
}

## calculate true conditional cross-covariance graph for the dynamic mixture model
Omega.true <- G.true <- G.neg <- rep(list(diag(p)),T)
for(t in 1:T){
  Omega <- Reduce('+',sapply(1:5, function(m) phi.true[t,m]^2*Omega.all[[m]],simplify=FALSE))
  Omega.true[[t]] <- Omega
  G.true[[t]] <- (abs(Omega) > 0.1)*1
  G.neg[[t]] <- (abs(Omega) < 0.001)*1
}


#### run Bayesian functional graph

MCMCspecs = list(B=1000,thin=5,burnin=1000,update=1000);
D_prior = list(a=0.1,b=0.1); lam_prior = list(a=0.1,b=1);
## B-spline basis functions
L = 5 # L = 11 for high-frequency
bspl <- create.bspline.basis(rangeval=c(0,T),nbasis=L,norder=4)
Phi <- t(eval.basis(1:T,bspl))

fgbay.samp <- fgraph_ccc(y,Phi,MCMCspecs,D_prior,lam_prior,ncor=1)

fgbay.est <- array(diag(p),c(p,p,T))
for(t in 1:T) fgbay.est[,,t][upper.tri(diag(p))] <- apply(fgbay.samp$ccc[,t,],1,function(x) (quantile(x,0.025)>0 | quantile(x,0.975)<0)*1)


#### calculate TPs and FPs

tmp <- matrix(FALSE,p,p); tmp[upper.tri(tmp)] <- TRUE
uptri.ind <- array(tmp,c(p,p,T))
sum((fgbay.est+G.true==2)[uptri.ind])/sum(G.true[uptri.ind]) 
sum((fgbay.est+G.neg==2)[uptri.ind])/sum(G.neg[uptri.ind])

    
