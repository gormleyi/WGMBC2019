rm(list=ls())
library(mclust)
library("MASS")
library(mvtnorm)
library(clusterGeneration)
library(Rfast)

# sim data
# set.seed(29042019)
# G.true = G = 2
# d = 100
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
# res = Mclust(x, G=2)
# table(cl.true, res$cl)
# adjustedRandIndex(cl.true, res$cl)

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
x = as.matrix(x[ind,1:150])
#x = as.matrix(x[ind,])
# # # How does mclust do?
res = Mclust(x, length(table(cl.true)), "VVV")
table(cl.true, res$cl)
adjustedRandIndex(cl.true, res$cl)

# Large wine data
# setwd("/Volumes/GoogleDrive/My Drive/Talks&Posters/ClassificationSocietyMeeting_2019/Code/CHIME_paper_data")
# load("/Users/clairegormley/Dropbox/STCO_MBCbigP_Paper/Data/winedata.Rdata")
# x = wine
# classes = gsub("[0-9]", "", winelabels)
# cl.true = rep(0, nrow(x))
# cl.true[which(classes == "ERA")] = 1
# cl.true[which(classes == "GRI")] = 2
# cl.true[which(classes == "OLO")] = 3
# G.true = 3
# # # # # How does mclust do?
# res = Mclust(x, length(table(cl.true)), modelNames = "VVV")
# table(cl.true, res$cl)
# adjustedRandIndex(cl.true, res$cl)

# Small wine data
# library(rattle.data)
# data(wine)
# x = as.matrix(wine[,-1])
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
block = 2
margins = combn(d,block)
p = ncol(margins) # number of combinations
denom = matrix(NA, n, p)
Gmax = 4
Rmax = 100 # Number of random starts to use.
loglik_store_init = rep(NA, Rmax)
cCLCstore = rep(NA, Gmax)
cBICstore = rep(NA, Gmax)
ICLstore = rep(NA, Gmax)
AWEstore = rep(NA, Gmax)
adjRstore = rep(NA, Gmax)
bestModel = list()
cbestModel =  list()
adjRbestModel =  list()
awebestModel =  list()
np = n*p
d1 = d-1
epsilon = 10^(-2)
itermax = 3000
beta.anneal.min = 1
beta.anneal.const = 1.05
beta.anneal.vals = rep(beta.anneal.min, 100)
for(j in 2:length(beta.anneal.vals))
{
  beta.anneal.vals[j] = beta.anneal.vals[j-1]*beta.anneal.const
}
beta.anneal.vals[which(beta.anneal.vals > 1)] = 1
beta.anneal.vals = beta.anneal.vals[1:(length(beta.anneal.vals) - sum(beta.anneal.vals==1) +1)]

for(G in 1:Gmax)
{
  print(data.frame(G))
  tau =rep(NA, G)
  mu.partition = array(0, c(G, block, p))
  Sigma.partition = array(0, c(block,block,p,G))
  Z = array(0, dim=c(n,G,p)) 
  Z.init = array(0, dim=c(n, G, p, Rmax))
  loglik_store_init = rep(NA, Rmax)
  cfit = array(0, c(n, G, p))
  
  if(G==1){
    Z = array(1, dim=c(n,G,p)); r = Rmax; rm(Z.init);
  }else{
    for(r in 1:Rmax) # If doing random starts run several
    {
      #initialization of the "latent" matrices, one for each bivariate margin:Check value of log likelihood and go with the associated Z 
      res = kmeans(x,G)$cl
      #res = cl.true
      #res = sample(1:G, n, replace=TRUE)
      #res = Mclust(x, G, "VVV")$cl
      for(i in 1:n){Z.init[i, res[i], ,r] = 1} # Use lots of random starts, same across all blocks
      # for(j in 1:p)
      # {
      #   for(i in 1:n)
      #   {
      #     res = Mclust(x[,margins[,j]], G, "VVV", verbose=FALSE)$cl # Bivariate mclust starts, different across all blocks
      #     Z.init[i, res[i], j,r] = 1
      #   }
      # }
    } #close Rmax. # Pull out the best starting Z and go from there
    dimZ.init = dim(Z.init)[4]
    k=0 # pyramid burnin scheme factor
    while(dimZ.init > 1)
    {
      for(r in 1:dimZ.init)
      {
        for(step in 1:(2^k))
        {
          # M step # 
          tau = apply(Z.init[,,,r], 2, sum)/np
          for(j in 1:p)
          {
            njg = colsums(as.matrix(Z.init[,,j,r]))
            for(g in 1:G)
            {
              mu.partition[g,,j] = apply(Z.init[,g,j,r]*x[,margins[,j]], 2, sum)/njg[g]
              mu.mat = matrix(mu.partition[g, , j], n, block, byrow=T)
              mult = x[,margins[,j]] - mu.mat
              Sigma.partition[,,j,g] = (t(Z.init[,g,j,r]*mult)%*%mult)/njg[g]
            }
            dens = cdensVVV(x[,margins[,j]], logarithm = TRUE, parameters = 
            list(mean = t(mu.partition[,,j]), variance = list(
             cholsigma = array(apply(Sigma.partition[,,j,], 3, chol), c(block,block,G)))))
            dens = sweep(dens, 2, log(tau), "+")
            zMax = apply(dens, 1, max)
            denom[,j] = zMax + log(rowSums(exp(dens - zMax)))
          } # closes j
        }
        loglik_store_init[r] = sum(denom)
      } #close for r
      Z.init = Z.init[,,,order(loglik_store_init, decreasing=T)[1:floor(dimZ.init/2)]]
      k = k+1 # turn on or off pyramid burnin
      if(dimZ.init > 1){ loglik_store_init = loglik_store_init[1:floor(dimZ.init/2)]}
      dimZ.init = floor(dimZ.init/2)
    }
    Z = Z.init
    rm(Z.init, loglik_store_init)
  }
  
  ######################
  #### EM algorithm ####
  ######################
  loglik_store = NULL
  llinf = NULL
  
  for(b in 1:length(beta.anneal.vals))
  {
    beta.anneal = beta.anneal.vals[b]
    if(G==1){beta.anneal = beta.anneal.vals[length(beta.anneal.vals)]; b = length(beta.anneal.vals)}
    changeinll = epsilon*2
    iter = 0
    
    #while(changeinll>=epsilon)
    while(iter<itermax)
    {
      iter  = iter+1
      if(iter%%100 == 0) print(iter)
      
      # M step: tau.# 
      tau = apply(Z, 2, sum)/np
      
      # Estimation of mu_g, Sigma_g and E-step#
      for(j in 1:p)
      {
        njg = colsums(as.matrix(Z[,,j]))
        for(g in 1:G)
        {
          mu.partition[g,,j] = apply(Z[,g,j]*x[,margins[,j]], 2, sum)/njg[g]
          mu.mat = matrix(mu.partition[g,,j], n, block, byrow=T)
          mult = x[,margins[,j]] - mu.mat
          Sigma.partition[,,j,g] = (t(Z[,g,j]*mult)%*%mult)/njg[g]
        } #closes g
        
        if(G>1){
          dens = cdensVVV(x[,margins[,j]], logarithm = TRUE, parameters = 
          list(mean = t(mu.partition[,,j]), variance = list(
           cholsigma = array(apply(Sigma.partition[,,j,], 3, chol), c(block,block,G)))))
        }else{
          dens = cdensVVV(x[,margins[,j]], logarithm = TRUE, parameters = 
           list(mean = mu.partition[,,j], variance = list(
           cholsigma = array(chol(Sigma.partition[,,j,1]), c(block,block,G)))))
        }
        dens = sweep(dens, 2, log(tau), "+")
        dens = sweep(dens, 2, beta.anneal, "*")
        zMax = apply(dens, 1, max)
        denom[,j] = zMax + log(rowSums(exp(dens - zMax)))
        Z[,,j] = exp(dens - denom[,j])
      } # closes j
      
      loglik_store = c(loglik_store, sum(denom))
      if(G == 1) iter = itermax
      # changeinll = (loglik_store[iter] - loglik_store[iter-1])/loglik_store[iter]
      # changeinll = loglik_store[iter] - loglik_store[iter-1]
      # if(iter>3)
      # {
      #  const = (loglik_store[iter] - loglik_store[iter-1])/(loglik_store[iter-1] - loglik_store[iter-2])
      #  llinf[iter] = loglik_store[iter-1] + (loglik_store[iter] - loglik_store[iter-1])/(1-const)
      #  if(iter>4)changeinll = abs(llinf[iter] - llinf[iter-1])
      #  if(const == 1) changeinll = 0
      # }
      #if(G==1)changeinll = epsilon-1
    }# end while convergence check
    print(paste("Finished beta.anneal = ", beta.anneal))
    
    #plot(1:iter, loglik_store, type="l", main=paste("G = ", G))
    
    if(G==1){
      ccl = rep(1, n)
    }else{
      # Classify observations based on CMAP fit a la Ranalli and Rocci (2017)
      for(j in 1:p)
      {
        for(g in 1:G)
        {
          cfit[,g,j] = log(tau[g]) + dmvnorm(x[,margins[,j]], mu.partition[g,,j], Sigma.partition[,,j,g], logged=TRUE)
        }
      }
      cfit = apply(cfit, c(1,2), sum)
      #cfit = sweep(exp(cfit), 1, apply(exp(cfit), 1, sum), "/")
      ccl = apply(cfit, 1, which.max)
    }
    print(table(ccl, cl.true))
    adjRstore[G] = adjustedRandIndex(ccl, cl.true)
    print(adjRstore[G])
  } # end for beta.anneal.vals
  
  cnoparam = (G-1) + G*block*p + G*((block*(block+1))/2)*p
  cBICstore[G] = -2*loglik_store[length(loglik_store)] + cnoparam*log(n)
  EN = Z*log(Z)
  EN[which(is.na(EN), arr.ind = T)] = 0
  cCLCstore[G] = -2*loglik_store[length(loglik_store)] + 2*sum(EN)
  ICLstore[G] = -2*loglik_store[length(loglik_store)] + 2*sum(EN) + cnoparam*log(n)
  AWEstore[G] = -2*loglik_store[length(loglik_store)] + 2*sum(EN) + (2*cnoparam)*(3/2 + log(n))
  rm(EN)
  
  # Save best model overall
  if(G == 1){bestcCLC = cCLCstore[G]+1; cbestBIC = cBICstore[G]+1; bestadjR = adjRstore[G]-1; bestAWE = AWEstore[G]+1}
  if(!is.na(cCLCstore[G])){
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
    }}
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
res = Mclust(x, G.true, "VVV")
res
table(cl.true, res$cl)
adjustedRandIndex(cl.true, res$cl)

table(cl.true, bestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, bestModel[["ccl"]])
AdjRand

table(cl.true, cbestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, cbestModel[["ccl"]])
AdjRand

table(cl.true, awebestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, awebestModel[["ccl"]])
AdjRand

table(cl.true, adjRbestModel[["ccl"]])
AdjRand = adjustedRandIndex(cl.true, adjRbestModel[["ccl"]])
AdjRand

# Pick best model
plot(cBICstore, type="b", main = paste("Best model: G = ", which.min(cBICstore[which(cBICstore!=0)])))
plot(cCLCstore, type="b", main = paste("Best model: G = ", which.min(cCLCstore)))
plot(ICLstore, type="b", main = paste("Best model: G = ", which.min(ICLstore)))
plot(AWEstore, type="b", main = paste("Best model: G = ", which.min(AWEstore)))




library(DescTools)
ents = apply(Z, 3, Entropy)
plot(sort(ents), type="b")
# Classify observations based on CMAP fit a la Ranalli and Rocci (2017)
cfit = array(0, c(n, G, 50))
ind = order(ents, decreasing=F)
for(j in 1:10)
{
  k = ind[j]
  for(g in 1:G)
  {
    cfit[,g,j] = log(tau[g]) + dmvnorm(x[,margins[,k]], mu.partition[g, margins[,k], k], Sigma.partition[margins[,k],margins[,k],k,g], logged=TRUE)
  }
  # if(j > 1)
  # {
  #  temp = matchClasses(table(apply(cfit[,,1], 1, which.max),apply(cfit[,,j], 1, which.max)), verbose=FALSE)
  #  cfit[,,j] = cfit[,as.numeric(temp),j]
  # }
}
cfit = apply(cfit, c(1,2), sum)
#cfit = sweep(exp(cfit), 1, apply(exp(cfit), 1, sum), "/")
reduced.ccl = apply(cfit, 1, which.max)
table(reduced.ccl, cl.true)

library(pracma)
cZ = apply(awebestModel$Z, 3, map)
cZ = apply(cZ, 1, Mode)
