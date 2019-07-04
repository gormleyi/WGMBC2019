rm(list=ls())
gc()
library(mclust)
library("MASS")
library(mvtnorm)
library(clusterGeneration)
library(Rfast)

# sim data
# set.seed(29042019)
# G.true = G = 2
# d = 12
# n = 90
# mu_true = matrix(c(rep(-2, d), rep(2, d)), G, d, byrow=T)
# mu_true = mu_true
# #irisBIC <- mclustBIC(iris[,-c(4,5)])
# #irisModel <- mclustModel(iris[,-c(4,5)], irisBIC)
# #Sigma_true = irisModel$parameters$variance$sigma
# Sigma_true = array(0, c(d,d,G))
# for(g in 1:G)
# {
#   Sigma_true[,,g] = genPositiveDefMat("unifcorrmat",dim=d, rangeVar=c(1,10))$Sigma
# }
# x = matrix(NA, n, d)
# x[1:(n/3),] = mvrnorm(n/3, mu_true[1,], Sigma_true[,,1])
# x[((n/3)+1):n,] = mvrnorm((2*n)/3, mu_true[2,], Sigma_true[,,2])
# cl.true = c(rep(1, n/3), rep(2, (2*n)/3))
# pairs(x[,1:5], col=cl.true)
# # How does mclust do?
# res = Mclust(x, G, modelName="VVV")
# table(cl.true, res$cl)
# adjustedRandIndex(cl.true, res$cl)

# Iris data
# x = as.matrix(iris[,-c(5)])
# cl.true = iris[,5]
# G.true = length(table(iris[,5]))
# irisBIC = mclustBIC(iris[,-c(5)], G=3)
# irisModel = mclustModel(iris[,-c(5)], irisBIC)
# Sigma_true = irisModel$parameters$variance$sigma
# mu_true = irisModel$parameters$mean
# tau_true = irisModel$parameters$pro
# # # How does mclust do?
# res = Mclust(x, length(table(cl.true)), modelName="VVV")
# table(cl.true, res$cl)
# adjustedRandIndex(cl.true, res$cl)

# Bank note data
# data(banknote)
# x = as.matrix(banknote[,-1])
# cl.true = banknote[,1]
# G.true = 2
# res = Mclust(x, G.true, "VVV")
# table(res$cl, cl.true)
# adjustedRandIndex(res$cl, cl.true)

# Kayee Yeung gene expression data
# setwd("/Volumes/GoogleDrive/My Drive/MastersThesis/Despoina Stamatopoulou/April2019")
# dat.all = read.table("GeneExpressionData.csv",header=TRUE,sep=",",dec=".")
# x = dat.all[,-c(1,2)]
# cl.true = dat.all[,2]
# # # How does mclust do?
# res = Mclust(x, length(table(cl.true)), modelName="VVV")
# table(cl.true, res$cl)
# adjustedRandIndex(cl.true, res$cl)

# CHIME GBM data
setwd("/Volumes/GoogleDrive/My Drive/Talks&Posters/ClassificationSocietyMeeting_2019/Code/CHIME_paper_data")
x = read.table("reducedx.txt")
cl.true = read.table("cltrue.txt")
cl.true = c(t(cl.true))
ind = c(which(cl.true==2), which(cl.true==3))
cl.true = cl.true[ind] - 1
G.true = length(table(cl.true))
x = as.matrix(x[ind,1:50])
#x = as.matrix(x[ind,])
# # # How does mclust do?
res = Mclust(x, length(table(cl.true)))
table(cl.true, res$cl)
adjustedRandIndex(cl.true, res$cl)

# Small wine data
# library(rattle.data)
# data(wine)
# x = as.matrix(wine[,-c(1:2)])
# G.true = 3
# cl.true = wine[,1]
# # # # # # How does mclust do?
# res = Mclust(x, length(table(cl.true)), modelNames = "VVV")
# table(cl.true, res$cl)
# adjustedRandIndex(cl.true, res$cl)


###############################################################################
## Set up.
n = nrow(x) 
d = ncol(x) 
Bs = 2   # block size (for now all the blocks same size)
nB = d/Bs  # number of blocks
blocks = matrix(1:d, d/Bs, Bs, byrow = TRUE)
Gmax = 4
Rmax = 100 # Number of random starts to use.
loglik_store_init = rep(NA, Rmax)
loghoodB.A = matrix(NA, n, nB)
cCLCstore = rep(NA, Gmax)
cBICstore = rep(NA, Gmax)
ICLstore = rep(NA, Gmax)
AWEstore = rep(NA, Gmax)
adjRstore = rep(NA, Gmax)
bestModel  = list()
cbestModel = list()
adjRbestModel = list()
awebestModel =  list()
nnB = n*nB
d1 = d-1
epsilon = 10^(-10)
itermax = 1000
setwd("/Volumes/GoogleDrive/My Drive/Talks&Posters/WGMBC2019/WGMBC2019")
source("CompositeConditionalApproach_Functions.R")

for(G in 1:Gmax)
{
  print(data.frame(G))
  tau = rep(NA, G)
  mu.partition = array(0, c(G, Bs, nB))
  Sigma.partition = array(0, c(2*Bs, 2*Bs, nB,G))
  E = array(NA, c(Bs, Bs, G, nB))
  muCond = array(NA, c(n, Bs, G, nB))
  Z = array(0, dim=c(n, G, nB)) 
  Z.init = array(0, dim=c(n, G, nB, Rmax))
  loglik_store_init = rep(NA, Rmax)
  cfit = array(0, c(n, G, nB))
  
  #pass = 1
  # while(pass <= nB)
  # {
  # if(pass == 1)
  # {
  #initialization of the "latent" matrices, one for each conditional block 
  if(G==1){ Z = array(1, dim=c(n,G,nB)); r = Rmax; rm(Z.init)}else{
    for(r in 1:Rmax) # If doing random starts run several: pick z giving highest ll
    {
      #initialization of the "latent" matrices, one for each bivariate margin:Check value of log likelihood and go with the associated Z 
      #res = kmeans(x,G)$cl
      res = cl.true
      for(i in 1:n){Z.init[i, res[i], ,r] = 1}# Use kmeans starts
    } #r
    dimZ.init = dim(Z.init)[4]
    k=0 # pyramid burnin scheme factor
    while(dimZ.init > 1)
    {
      for(r in 1:dimZ.init)
      {
        for(step in 1:(3^k))
        {
          tau = mstep.tau(Z.init[,,,r], nnB)
          # # Estimation of means and covs: do single marginal block first, then loop over others
          out = emstep.marginalblock(x, blocks, Z.init[,,,r], tau, mu.partition, Sigma.partition, loghoodB.A, n, Bs, G)
          mu.partition = out$mu.partition 
          Sigma.partition = out$Sigma.partition 
          loghoodB.A = out$loghoodB.A 
          densA = out$densA
          
          out = emstep.conditionalblocks(x, Z.init[,,,r], tau, mu.partition, muCond, Sigma.partition, E, loghoodB.A, n, Bs, G, iter=1)
          Z.init[,,,r] = out$Z
          tau = out$tau 
          mu.partition = out$mu.partition 
          muCond = out$muCond 
          Sigma.partition = out$Sigma.partition 
          loghoodB.A = out$loghoodB.A 
          E = out$E
          
          loglik_store_init[r] = sum(loghoodB.A)
        } 
      } #close r
      Z.init = Z.init[,,,order(loglik_store_init, decreasing=T)[1:floor(dimZ.init/2)]]
      k = k+1 # turn on or off pyramid burnin
      if(dimZ.init > 1) loglik_store_init = loglik_store_init[1:floor(dimZ.init/2)]
      dimZ.init = floor(dimZ.init/2)
    } # while dimZ.init
    Z = Z.init
    rm(Z.init, loglik_store_init)
  }
  # }else{
  #  blocks = blocks[sample(1:nB, nB, replace=F),]
  # } # which pass
  
  ######################  
  #### EM algorithm ####
  ######################
  loglik_store = NULL
  llinf = NULL
  changeinll = epsilon*2
  iter = 0
  #while(5*changeinll>=epsilon)
  while(iter<itermax)
  {
    iter  = iter+1
    if(iter%%100 == 0) print(iter)
    
    # M step # 
    # Estimation of tau #
    tau = mstep.tau(Z, nnB)
    
    if(iter == 1)
    {
     # # Estimation of means and covs: do single marginal block first, then loop over others
     out = emstep.marginalblock(x, blocks, Z, tau, mu.partition, Sigma.partition, loghoodB.A, n, Bs, G)
     Z = out$Z; mu.partition = out$mu.partition; Sigma.partition = out$Sigma.partition;
     loghoodB.A = out$loghoodB.A; densA = out$densA
    }
    
    out = emstep.conditionalblocks(x, Z, tau, mu.partition, muCond, Sigma.partition, E, loghoodB.A, n, Bs, G, iter)
    Z = out$Z; tau = out$tau; mu.partition = out$mu.partition; muCond = out$muCond;
    Sigma.partition = out$Sigma.partition; loghoodB.A = out$loghoodB.A; E = out$E
    loglik_store = c(loglik_store, sum(loghoodB.A))
    
    if(G == 1) iter = itermax
    # if(iter>3)
    # {
    #   const = (loglik_store[iter] - loglik_store[iter-1])/(loglik_store[iter-1] - loglik_store[iter-2])
    #   llinf[iter] = loglik_store[iter-1] + (loglik_store[iter] - loglik_store[iter-1])/(1-const)
    #   if(iter>4)changeinll = abs(llinf[iter] - llinf[iter-1])
    #   if(const == 1) changeinll = 0
    # }
    # if(G==1)changeinll = epsilon-1
  }# end while iter/convergence
  
  if(G==1){
    ccl = rep(1, n)
  }else{
    cfit[,,1] = densA
    for(j in 2:nB)
    {
      if(G == 1){
        cholE = array(chol(E[,,,j]), c(Bs, Bs, G))
        cfit[,,j] = matrix(sapply(1:n, function(i) cdensVVV(x[i,blocks[j,],drop = FALSE], log = TRUE,
         parameters = list(mean = muCond[i,,], variance = list(cholsigma = cholE)))),n, G)
      }else{
        cholE = array(apply(E[,,,j], 3, chol), c(Bs, Bs, G))
        cfit[,,j] = t(sapply(1:n, function(i) cdensVVV(x[i,blocks[j,],drop = FALSE], log = TRUE,
          parameters = list(mean = muCond[i,,,j], variance = list(cholsigma = cholE)))))
      }
      cfit[,,j] = sweep(matrix(cfit[,,j], n, G), 2, log(tau), "+")
    }
    cfit = apply(cfit, c(1,2), sum)
    cfit = sweep(cfit, 1, apply(cfit, 1, sum), "-")
    ccl = apply(cfit, 1, which.max)
  }
  print(table(ccl, cl.true))
  adjRstore[G] = adjustedRandIndex(ccl, cl.true)
  print(adjRstore[G])
  
  cnoparam = (G-1) + G*Bs*nB + G*((Bs*(Bs+1))/2)*nB
  cBICstore[G] = cBIC = -2*loglik_store[length(loglik_store)] + cnoparam*log(n)
  EN = Z*log(Z)
  EN[which(is.na(EN), arr.ind = T)] = 0
  cCLCstore[G] = -2*loglik_store[length(loglik_store)] + 2*sum(EN)
  ICLstore[G] = -2*loglik_store[length(loglik_store)] + 2*sum(EN) + cnoparam*log(n)
  AWEstore[G] = -2*loglik_store[length(loglik_store)] + 2*sum(EN) + (2*cnoparam)*(3/2 + log(n))
  
  #pass = pass+1
  #} # Close passes
  
  # Save best model overall
  if(G == 1){bestcCLC = cCLCstore[G]+1; cbestBIC = cBICstore[G]+1; bestadjR = adjRstore[G]-1; bestAWE = AWEstore[G]+1}
  if(cCLCstore[G] < bestcCLC)
  {
    print(paste("New best group via cCLC: G = ", G))
    bestcCLC = cCLCstore[G]
    bestModel[["Z"]] = Z
    bestModel[["tau"]] = tau
    bestModel[["mu"]] = mu.partition
    bestModel[["Sigma"]] = Sigma.partition
    bestModel[["ccl"]] = ccl
    bestModel[["G"]] = G
  }
  if(cBICstore[G] < cbestBIC)
  {
    print(paste("New best group via cBIC: G = ", G))
    cbestBIC = cBICstore[G]
    cbestModel[["Z"]] = Z
    cbestModel[["tau"]] = tau
    cbestModel[["mu"]] = mu.partition
    cbestModel[["Sigma"]] = Sigma.partition
    cbestModel[["ccl"]] = ccl
    cbestModel[["G"]] = G
  }
  if(AWEstore[G] < bestAWE)
  {
    print(paste("New best group via AWE: G = ", G))
    bestAWE = AWEstore[G]
    awebestModel[["Z"]] = Z
    awebestModel[["tau"]] = tau
    awebestModel[["mu"]] = mu.partition
    awebestModel[["Sigma"]] = Sigma.partition
    awebestModel[["ccl"]] = ccl
    awebestModel[["G"]] = G
  }
  if(adjRstore[G] > bestadjR)
  {
    bestadjR = adjRstore[G]
    print(paste("New best group based on adjR", adjRstore[G], "G = ", G))
    adjRbestModel[["Z"]] = Z
    adjRbestModel[["tau"]] = tau
    adjRbestModel[["mu"]] = mu.partition
    adjRbestModel[["Sigma"]] = Sigma.partition
    adjRbestModel[["ccl"]] = ccl
    adjRbestModel[["G"]] = G
  }
} # Close Gmax loop

# How does mclust do?
res = Mclust(x, 1:Gmax, "VVV")
res
table(cl.true, res$cl)
adjustedRandIndex(cl.true, res$cl)

# Pick best model
plot(cBICstore, type="b", main = paste("Best model: G = ", which.min(cBICstore[which(cBICstore!=0)])))
plot(cCLCstore, type="b", main = paste("Best model: G = ", which.min(cCLCstore)))
plot(ICLstore, type="b", main = paste("Best model: G = ", which.min(ICLstore)))
plot(AWEstore, type="b", main = paste("Best model: G = ", which.min(AWEstore)))


table(cl.true, bestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, bestModel[["ccl"]])
AdjRand

table(cl.true, cbestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, cbestModel[["ccl"]])
AdjRand

table(cl.true, adjRbestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, adjRbestModel[["ccl"]])
AdjRand

table(cl.true, awebestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, awebestModel[["ccl"]])
AdjRand

