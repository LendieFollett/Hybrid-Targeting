#--calculate pair.comp_{ij} = 1{Y_i < Y_j} --------------------------------
#' Compute Pairwise Comparison Matrix for Full Ranking List of the Entities
#'
#' Compute the pairwise comparison matrix from the ranking lists of the ranked entities.
#' @param rank.vec A full ranking list containing the ranks of all the \eqn{N} entities. Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @return An \eqn{N} by \eqn{N} pairwise comparison for all \eqn{N} entities, where the (\eqn{i},\eqn{j}) element equals 1 if \eqn{i} is ranked higher than \eqn{j}, and 0 if \eqn{i} is ranker lower than \eqn{j}. Note that the diagonal elements (\eqn{i},\eqn{i})'s are set to NA.
#' @export
FullRankToPairComp <- function( rank.vec, n = length(rank.vec) ){
  pair.comp <- matrix(NA, n, n)
  for(i in 1:n){
    j = which(rank.vec == i)
    pair.comp[j,  rank.vec > i] = 1
    pair.comp[j,  rank.vec < i] = 0
  }
  return(pair.comp)
}


### Gibbs update Z given beta--------------------------
### pair.comp.ten[,,j]: pairwise comparison matrix for jth ranker ###
### Z.mat[,j]: latent variable vector for jth ranker ###
### mu: shared mean vector for testing sample N1 x 1 ###
### weight.vec[j]: weight for jth ranker ###
GibbsUpLatentGivenRankGroup <- function(pair.comp.ten, Z, mu, omega_rank = rep(1, ncol(Z)), R = ncol(Z) ){
  for(r in 1:R){ #loop over rankers (e.g., CBT, geography, etc...)
    Z[,r] = GibbsUpLatentGivenRankInd(pair.comp.ten[,,r], Z[,r], mu, weight = omega_rank[r])
  }
  return(Z)
}


GibbsUpLatentGivenRankInd <- function(pair.comp, Z, mu, weight){
  up.order = sort( rowSums( pair.comp, na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix
  for(i in up.order){
    
    set1 = which( pair.comp[i, ] == 1)
    set0 = which( pair.comp[i, ] == 0)
    
    if(length(set1) > 0){
      upper = min(Z[set1])
    }else{
      upper = Inf
    }
    
    if(length(set0) > 0){
      lower = max(Z[set0])
    }else{
      lower = -Inf
    }
    
    Z[i] = rtruncnorm( 1, lower, upper, mean = mu[i], sd = 1/sqrt(weight) )
  }
  return(Z)
}



# Gibbs update for beta (includes intercept) given Z, y_micro, y_comm---------------------

### Gibbs update for the shared mean mu ###
### Z.mat is a N1 x R matrix ###
### Z.mat[,j]: latent variable vector for jth ranker j = 1, ..., R ###
### X_comm is a KxP matrix (training dataset, aggregated by community)###
### X_micro is a N0xP matrix (training dataset, on the household level) ###
### X_micro is a N1xP matrix (testing dataset, on the household level) ###
### weight.vec: (A + M + R)x1 vector of weights omega_rank, omega_comm, omega_micro###
### omega_rank[r] = weight of r^th ranker###
### sigma2.alpha, sigma2.beta: prior parameters for mu = (alpha, beta) ###
### para.expan: whether use parameter expansion  LRF: FALSE???###
GibbsUpMuGivenLatentGroup <- function(Z, Y_comm=NA, Y_micro=NA, #<-- 3 "response" matrices
                                      X_comm=NA, X_micro0=NA, X_micro1=NA, 
                                      omega_comm = rep(1,ncol(Y_comm)), 
                                      omega_micro = rep(1, ncol(Y_micro)),
                                      omega_rank = rep(1, ncol(Z)),
                                      sigma2.beta = 1){
  
  #LRF TO ADDRESS: Y_comm might be missing, Y_micro might be missing, ...assuming ranking will be there...
  #                same logic for corresponding x matrices
  #X_comm isn't just for training... need to fix that
  #dimensions
  A <- ncol(Y_comm)
  M <- ncol(Y_micro)
  R <- ncol(Z)
  K <- nrow(Y_comm)
  N0 <- nrow(Y_micro)
  N1 <- nrow(Z)
  if(all(X_comm[,1]==1)){
    P <- ncol(X_micro)-1
  }else{
    P <- ncol(X_micro)
  }
  
  #prior on beta vector - mean = 0, variance = sigma2.beta
  sigma2.beta <- 2.5 
  
  #Complete 'data' vector
  u <- c(as.vector(Y_comm), #KxA --> (AK)x1
             as.vector(Y_micro), #N0xM --> (M*N0)x1
             as.vector(Z)) #N1xR --> (R*N1)x1
#length(u) == A*K + M*N0 +R*N1 check
  
  ### X.mat it full, standardized X matrix with training first, then testing - constructed by kroneker products ###
  ### X.mat is a (A*K + M*N0 +R*N1)x(P+1) matrix ###
  X <- rbind(kronecker(rep(1, A), X_comm), #(A*K)xP
                 kronecker(rep(1, M), X_micro0),#(M*N0)xP
                 kronecker(rep(1, R), X_micro1))#(R*N1)xP
  
  Sigma_inv_diag <- 1/c(rep(omega_comm, each = K), 
                      rep(omega_micro, each = N0),
                      rep(omega_rank, each = N1))
  
  #A<-1x(A*K + M*N0 +R*N1)%*%square(A*K + M*N0 +R*N1)%*%(A*K + M*N0 +R*N1)xP --> 1xP
  pt1 <- u^T%*%diag(Sigma_inv_diag)%*%X
  
  pt2 <- t(X)%*%diag(Sigma_inv_diag)%*%X + diag(P+1)/sigma2.beta
  
  pt2_inv <- solve(pt2)
  
  alpha_beta <- mvrnorm(1, mu = pt1%*%pt2_inv, Sigma = pt2_inv)

  return(alpha_beta)
}



### Gibbs update for information quality weights omega_comm, omega_micro, omega_rank---------
#y is either Y_comm (KxA), Y_micro (N0xM), or Z (N1xR)
GibbsUpQualityWeights <- function(y, mu, weight.prior.value = c(0.5, 1, 2), weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)), Col = ncol(y), Row = nrow(y) ){
  n.prior.value <- length(weight.prior.value)
  weight_samp <- rep(NA, Col)

for( col in 1:Col){ #over information source within y
  log.post.prob = rep(0, n.prior.value) #re-initialize for next information source
  for(k in 1:n.prior.value){ #over potential values
    log.post.prob[k] = log.post.prob[k]+
      log( weight.prior.prob[k] ) + Row/2 * log( weight.prior.value[k] ) - weight.prior.value[k]/2 * sum( (y[,col] - mu)^2 )
    #log(prior value) - log(sigma) -(1/(2*sigma*sigma))*sum[(y-mu)^2]=
    #log(prior value) - .5log(1/w) -(w/2)*sum[(y-mu)^2]
  }

  log.post.prob = log.post.prob - max(log.post.prob)
  post.prob = exp(log.post.prob)
  
  weight_samp[col] <- weight.prior.value[ as.vector( rmultinom(1, 1, prob = post.prob) ) == 1 ]
  
}

  return(weight_samp)
}


#--to run BCRank model MCMC--------------------------------

#' Bayesian Consensus Targeting
#'
#' Implement the Bayesian model for mixed data inputs with covariate information and subjective + objective information with varying qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @import MASS
#' @param pair.com.ten An \eqn{N1} by \eqn{N1} by \eqn{R} pairwise comparison array for all \eqn{N1} entities and \eqn{R} rankers, where the (\eqn{i},\eqn{j},\eqn{r}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{r}, 0 if \eqn{i} is ranker lower than \eqn{j}, and NA if the relation between \eqn{i} and \eqn{j} is missing. Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{r})'s for all rankers should be set to NA as well.
#' @param X_micro0 An \eqn{N0} by \eqn{P+1} covariate matrix for the \eqn{N0} entities with \eqn{P} covariates. Assumes 1st col is 1's for intercept.
#' @param X_micro1 An \eqn{N1} by \eqn{P+1} covariate matrix for the \eqn{N1} entities with \eqn{P} covariates. Assumes 1st col is 1's for intercept.
#' @param X_comm An \eqn{K} by \eqn{P+1} covariate matrix. An aggregate of X_micro0 and X_micro 1. 
#' @param sigma_beta Currently given/fixed. Prior variance on beta.
#' @param weight.prior.value A vector for the support of the discrete prior on weight parameter.
#' @param weight.prior.prob A vector for the prior probability mass of the discrete prior on weight parameter.
#' @param iter.max Number of iterations for Gibbs sampler.
#' @return A list containing posterior samples of mu, the shared 'wellness' mean, conditional on the test X_micro1.
#' @export
BCTarget <- function(pair.comp.ten, X_micro0, X_micro1, X_comm,
                               Y_comm, Y_micro,
                               sigma_beta = 5^2,
                               weight.prior.value = c(0.5, 1, 2), 
                               weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                               N1 = dim(pair.comp.ten)[1], 
                               R = dim(pair.comp.ten)[3], 
                               iter.max = 5000, para.expan = TRUE, print.opt = 100,
                               initial.list = NULL){
  
  if(all(X_comm[,1]==1)){
    P <- ncol(X_comm)-1
  }else{
    P <- ncol(X_comm)
  }

  
  A <- ncol(Y_comm)
  M <- ncol(Y_micro)
  R <- dim(pair.comp.ten)[3]
  K <- nrow(Y_comm)
  N0 <- nrow(Y_micro)
  N1 <- dim(pair.comp.ten)[1]
  

  ## store MCMC draws
  draw = list(
    Z = array(NA, dim = c( iter.max,N1, R)),
    beta = array(NA, dim = c(iter.max,P+1)),
    mu = array(NA, dim = c(iter.max,N1)),
    omega_comm = array(NA, dim = c(iter.max, A) ),
    omega_micro = array(NA, dim = c(iter.max, M) ),
    omega_rank = array(NA, dim = c(iter.max, R) ),
    sigma2.beta = rep(NA, iter.max)
  )
  
  if(is.null(initial.list)){
    ## initial values for Z
    Z = matrix(NA, nrow = N1, ncol = R)
    for(j in 1:R){
      Z[sort( rowSums( pair.comp.ten[,,j], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix, j] = (c(N1 : 1) - (1+N1)/2)/sd(c(N1 : 1))
    }
    
    ## initial values for alpha, beta and thus mu
    beta = rep(0, P+1)
    mu = as.vector(X_micro1 %*% beta )
    
    ## initial values for weights
    omega_comm = rep(1, A) 
    omega_micro = rep(1, M) 
    omega_rank = rep(1, R)
    
    ## initial values for sigma2
    sigma2.beta = 2.5
  }else{
    
    Z = initial.list$Z
    beta = initial.list$beta
    mu = as.vector( X_micro1 %*% beta )
    omega_comm = initial.list$omega_comm
    omega_micro = initial.list$omega_micro
    omega_rank = initial.list$omega_rank
    sigma2.beta = initial.list$sigma2.beta
    
  }
  
  ## store initial value
  draw$Z[1,,] = Z
  draw$beta[1,] = beta
  draw$mu[1,] = mu
  draw$omega_comm[1,] = omega_comm
  draw$omega_micro[1,] = omega_micro
  draw$omega_rank[1,] = omega_rank
  
  
  ## Gibbs iteration
  for(iter in 2:iter.max){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, Z = Z, mu = mu, omega_rank = omega_rank, R = R )
    
    # update beta or equivalently mu given Z.mat
    beta = GibbsUpMuGivenLatentGroup(Z = Z, Y_comm = Y_comm, Y_micro = Y_micro,
                                                 X_comm = X_comm, X_micro0 = X_micro0, X_micro1 = X_micro1,
                                                 omega_comm=omega_comm, omega_micro = omega_micro, omega_rank = omega_rank,
                                                 sigma2.beta = 2.5)
    
    # update quality weights
    omega_comm <- GibbsUpQualityWeights(y=Y_comm, mu=X_comm %*% beta, weight.prior.value = c(0.5, 1, 2 ))
    omega_micro <-GibbsUpQualityWeights(y=Y_micro, mu=X_micro0 %*% beta, weight.prior.value = c(0.5, 1, 2 ))
    omega_rank <- GibbsUpQualityWeights(y=Z, mu=X_micro1 %*% beta, weight.prior.value = c(0.5, 1, 2 ))

    
    #LRF TO ADDRESS: this is to be computed with the 'connections' dummy 0'd out
    mu = as.vector( X_micro1 %*% beta )
    
    ### update weight
    #weight.vec = GibbsUpWeightGroup(Z.mat = Z.mat, mu = mu, weight.prior.value = weight.prior.value, weight.prior.prob = weight.prior.prob, n.item = n.item, n.ranker = n.ranker)
    
    ### update sigma2
    #sigma2.alpha = GibbsUpsigma2(alpha, nu.alpha, tau2.alpha)
    #if(p.cov > 0){
    #  sigma2.beta = GibbsUpsigma2(beta, nu.beta, tau2.beta)
    #}
    
    # store value at this iteration
    draw$Z[iter,,] = Z
    draw$beta[iter,] = beta
    draw$mu[iter,] = mu
    draw$omega_micro[iter,] = omega_micro
    draw$omega_comm[iter,] = omega_comm
    draw$omega_rank[iter,] = omega_rank
    # print iteration number
    if(iter %% print.opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(weight.vec))
      # print(c(sigma2.alpha, sigma2.beta))
    }
  }
  return(draw)
}


