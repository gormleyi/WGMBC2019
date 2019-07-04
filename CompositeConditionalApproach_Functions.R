#################################################################################
########## Functions required for conditional composite likelihood approach #####
#################################################################################

### Overall M step for tau
mstep.tau = function(Z, nnB)
{
  tau = apply(Z, 2, sum)/nnB
  tau
}

### EM segement for the initial marginal block
emstep.marginalblock = function(x, blocks, Z, tau, mu.partition, Sigma.partition, loghoodB.A, n, Bs, G)
{
  xA = x[,blocks[1,]]
  setA = 1:Bs
  setB = (Bs + 1):(2*Bs)
  ng = colSums(as.matrix(Z[,,1]))
  for(g in 1:G)
  {
    mu.partition[g,,1] = apply(Z[,g,1]*xA, 2, sum)/ng[g]
    mu.mat = matrix(mu.partition[g,,1], n, Bs, byrow=T)
    mult = xA - mu.mat
    Sigma.partition[setB,setB,1,g] = (t(Z[,g,1]*mult)%*%mult)/ng[g]
    Sigma.partition[setA,setA,2,g] = Sigma.partition[setB,setB,1,g]
  }
  if(G == 1){
    cholE = array(chol(Sigma.partition[setB,setB,1,]), c(Bs,Bs,1))
    densA = matrix(sapply(1:n, function(i) cdensVVV(xA[i,,drop = FALSE], log = TRUE,
    parameters = list(mean = mu.partition[,,1], variance = list(cholsigma = cholE)))), n, G)
  }else{
    cholE = array(apply(Sigma.partition[setB,setB,1,], 3, chol), c(Bs, Bs, G))
    densA = t(sapply(1:n, function(i) cdensVVV(xA[i,,drop = FALSE], log = TRUE,
    parameters = list(mean = t(mu.partition[,,1]), variance = list(cholsigma = cholE)))))
  }
  densA = sweep(densA, 2, log(tau), "+")
  tmp = apply(densA, 1, max)
  loghoodB.A[,1] = tmp + log(rowSums(exp(densA - tmp)))
  Z[,,1] = exp(densA - loghoodB.A[,1])
  list(Z = Z, mu.partition = mu.partition, Sigma.partition = Sigma.partition, loghoodB.A = loghoodB.A, densA = densA)
} # end emsteps.marginal block


### EM segment for the conditional blocks
emstep.conditionalblocks = function(x, Z, tau, mu.partition, muCond, Sigma.partition, E, loghoodB.A, n, Bs, G, iter)
{
  if(iter  == 1){
    start = 2
    xA <- x[,blocks[1,]]
  }else{
    start = 1
    xA <- x[,blocks[nB,]]
  }
  for(j in start:nB)
  {
    xB = x[,blocks[j,]]
    xAB = as.matrix(cbind(xA, xB))
    JAB = ncol(xAB)
    setA = 1:(JAB - Bs)
    setB = (JAB - Bs + 1):JAB
    
    O = covw(xAB, Z[,,j], normalize = FALSE)
    W = O$W[setA, setA,]
    V = O$W[setA, setB,]
    U = O$W[setB, setB,]
    
    if(j == 1){prev = nB}else{prev = j-1}
    
    for(g in 1:G)
    {
      #Efficient approach
      reg = lm.wfit(model.matrix(xAB[,setB] ~ xAB[,setA], data = as.data.frame(xAB)), xAB[,setB], Z[,g,j])
      beta = as.matrix(reg$coefficients)[-1,, drop = FALSE]       # t( t(C) %*% inv )
      alpha = as.matrix(reg$coefficients)[1,, drop = FALSE]
      C = crossprod(Sigma.partition[setA,setA,j,g], beta)
      mu.partition[g,,j] = colSums(sweep(beta, 1, mu.partition[g,,prev], "*")) + alpha
      E[,,g,j] = covw(residuals(reg), Z[,g,j], normalize = FALSE )$S[,,1]
      Sigma.partition[setB, setB, j, g] = E[,,g,j] + crossprod(C, solve(Sigma.partition[setA,setA,j,g])) %*% C
      Sigma.partition[setA, setB, j, g] = C
      Sigma.partition[setB, setA, j, g] = t(C)
      
      # mean of conditional distribution B|A
      tmp = sweep(xAB[,setA], 2, mu.partition[g,,prev], "-")
      muCond[,,g,j] = sweep(tmp %*% beta, 2, mu.partition[g,,j], "+")
    } # g loop
    
    
    # E step #
    # We create one latent-value matrix Z for every marginal composite likelihood component. #
    if(G == 1){
      cholE = array(chol(E[,,,j]), c(Bs, Bs, G))
      densB.A = matrix(sapply(1:n, function(i) cdensVVV(xB[i,,drop = FALSE], log = TRUE,
       parameters = list(mean = muCond[i,,,j], variance = list(cholsigma = cholE)))), n, G)
    }else{
      cholE = array(apply(E[,,,j], 3, chol), c(Bs, Bs, G))
      densB.A = t(sapply(1:n, function(i) cdensVVV(xB[i,,drop = FALSE], log = TRUE,
         parameters = list(mean = muCond[i,,,j], variance = list(cholsigma = cholE)))))
    }
    densproB.A = sweep(densB.A, 2, log(tau), "+")
    zMax = apply(densproB.A, 1, max)
    loghoodB.A[,j] = zMax + log(rowSums(exp(densproB.A - zMax)))
    Z[,,j] = exp(densproB.A - loghoodB.A[,j])
    
    xA = xB
    if(j<nB){Sigma.partition[setA,setA,(j+1),] = Sigma.partition[setB, setB, j, ]}else{
      Sigma.partition[setA,setA,1,] = Sigma.partition[setB, setB, j, ]}
  } # close j 
  # Once have passed through then run conditional of block nB on block 1
  
  
  list(Z = Z, tau = tau, mu.partition = mu.partition, muCond = muCond, Sigma.partition = Sigma.partition, loghoodB.A = loghoodB.A, E = E)
} # Close em step conditional blocks


### Classification function
classify = function(x, blocks, cfit, densA, tau, muCond, E, n, nB, G)
{
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
      parameters = list(mean = muCond[i,,], variance = list(cholsigma = cholE)))))
    }
    cfit[,,j] = sweep(matrix(cfit[,,j], n, G), 2, log(tau), "+")
  }
  cfit = apply(cfit, c(1,2), sum)
  #cfit = sweep(cfit, 1, apply(cfit, 1, sum), "-")
  ccl = apply(cfit, 1, which.max)
  list(ccl)
}


