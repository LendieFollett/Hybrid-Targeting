#--calculate pair.comp_{ij} = 1{Y_i < Y_j} --------------------------------
library(MASS) #for mvrnorm
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


### Gibbs update Z.mat given (alpha, beta)--------------------------
### pair.comp.ten[,,j]: pairwise comparison matrix for jth ranker ###
### Z.mat[,j]: latent variable vector for jth ranker ###
### mu: shared mean vector for this group of rankers ###
### weight.vec[j]: weight for jth ranker ###
GibbsUpLatentGivenRankGroup <- function(pair.comp.ten, Z.mat, mu, weight.vec = rep(1, ncol(Z.mat)), n.ranker = ncol(Z.mat) ){
  for(j in 1:n.ranker){ #loop over rankers (e.g., CBT, geography, etc...)
    Z.mat[,j] = GibbsUpLatentGivenRankInd(pair.comp.ten[,,j], Z.mat[,j], mu, weight = weight.vec[j])
  }
  return(Z.mat)
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



# Gibbs update for (alpha, beta) given Z.mat---------------------

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
GibbsUpMuGivenLatentGroup <- function(Z.mat, Y_comm=NA, Y_micro=NA, #<-- 3 "response" matrices
                                      X_comm=NA, X_micro=NA, X_1, 
                                      omega_comm = rep(1,ncol(Y_comm)), 
                                      omega_micro = rep(1, ncol(Y_micro)),
                                      omega_rank = ncol(Z.mat),
                                      sigma2.alpha = 2, sigma2.beta = 1, R = ncol(Z.mat), 
                                      n.item =nrow(Z.mat), p.cov = ncol(X.mat), para.expan = FALSE){
  
  #LRF TO ADDRESS: Y_comm might be missing, Y_micro might be missing, ...assuming ranking will be there...
  #                same logic for corresponding x matrices
  #X_comm isn't just for training... need to fix that
  #dimensions
  A <- ncol(Y_comm)
  M <- ncol(Y_micro)
  R <- ncol(Z.mat)
  K <- nrow(Y_comm)
  N0 <- nrow(Y_micro)
  N1 <- nrow(Z.mat)
  
  #Complete 'data' vector
  u <- rbind(as.vector(Y_comm), #KxA --> (AK)x1
             as.vector(Y_micro), #N0xM --> (M*N0)x1
             as.vector(Z.mat)) #N1xR --> (R*N1)x1

  
  ### X.mat it full, standardized X matrix with training first, then testing - constructed by kroneker products ###
  ### X.mat is a (A*K + M*N0 +R*N1)x(P+1) matrix ###
  X.mat <- rbind(kronecker(rep(1, A), X_comm), #(A*K)xP
                 kronecker(rep(1, M), X_micro),#(M*N0)xP
                 kronecker(rep(1, R), X_1))#(R*N1)xP
  
  Sigma_inv_diag <- 1/c(rep(omega_comm, each = K), 
                      rep(omega_micro, each = N0),
                      rep(omega_rank, each = N1))
  
  #A<-1x(A*K + M*N0 +R*N1)%*%square(A*K + M*N0 +R*N1)%*%(A*K + M*N0 +R*N1)xP --> 1xP
  pt1 <- u^T%*%diag(Sigma_inv_diag)%*%X.mat
  
  pt2 <- t(X.mat)%*%diag(Sigma_inv_diag)%*%X.mat + diag(A*K + M*N0 +R*N1)/sigma2.beta
  
  pt2_inv <- solve(pt2)
  
  alpha_beta <- mvrnorm(1, mu = pt1%*%pt2_inv, Sigma = pt2_inv)
  
  alpha <- alpha_beta[1]
  beta <- alpha_beta[-1]
  

  return(list(alpha = alpha, beta = beta))
}




#--to run BARCW model MCMC--------------------------------

#' Bayesian Analysis of Rank data with entities' Covariates and rankers' Weights.
#'
#' Implement the Bayesian model for rand data with ranked entities' covariates information and rankers' with varying qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @param pair.com.ten An \eqn{N} by \eqn{N} by \eqn{M} pairwise comparison tensor for all \eqn{N} entities and \eqn{M} rankers, where the (\eqn{i},\eqn{j},\eqn{m}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{m}, 0 if \eqn{i} is ranker lower than \eqn{j}, and NA if the relation between \eqn{i} and \eqn{j} is missing. Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{m})'s for all rankers should be set to NA as well.
#' @param X.mat An \eqn{N} by \eqn{L} covariate matrix for the \eqn{N} entities with \eqn{L} covariates.
#' @param tau2.alpha The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param nu.alpha The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_alpha}.
#' @param tau2.beta The scale parameter for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param nu.beta The d.f. for the scaled inverse chi-squared prior on \eqn{\sigma^2_beta}.
#' @param weight.prior.value A vector for the support of the discrete prior on weight parameter.
#' @param weight.prior.prob A vector for the probability mass of the discrete prior on weight parameter.
#' @param iter.max Number of iterations for Gibbs sampler.
#' @param para.expan Logical variable for whether using parameter expansion in the Gibbs sampler.
#' @return A list containing posterior samples of all the missing evaluation scores for all rankers and all the model parameters.
#' @export
BayesRankCovWeight <- function(pair.comp.ten, X.mat = matrix(NA, nrow =dim(pair.comp.ten)[1], ncol = 0),
                               tau2.alpha = 5^2, nu.alpha = 3,
                               tau2.beta = 5^2, nu.beta = 3,
                               weight.prior.value = c(0.5, 1, 2), weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                               n.item = dim(pair.comp.ten)[1], n.ranker = dim(pair.comp.ten)[3], p.cov = ncol(X.mat),
                               iter.max = 5000, para.expan = TRUE, print.opt = 100,
                               initial.list = NULL){
  ## store MCMC draws
  draw = list(
    Z.mat = array(NA, dim = c(n.item, n.ranker, iter.max)),
    alpha = array(NA, dim = c(n.item, iter.max)),
    beta = array(NA, dim = c(p.cov, iter.max)),
    mu = array(NA, dim = c(n.item, iter.max)),
    weight.vec = array(NA, dim = c(n.ranker, iter.max) ),
    sigma2.alpha = rep(NA, iter.max),
    sigma2.beta = rep(NA, iter.max)
  )
  
  if(is.null(initial.list)){
    ## initial values for Z
    Z.mat = matrix(NA, nrow = n.item, ncol = n.ranker)
    for(j in 1:n.ranker){
      Z.mat[sort( rowSums( pair.comp.ten[,,j], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix, j] = (c(n.item : 1) - (1+n.item)/2)/sd(c(n.item : 1))
    }
    
    ## initial values for alpha, beta and thus mu
    alpha = rep(0, n.item)
    beta = rep(0, p.cov)
    mu = as.vector( alpha + X.mat %*% beta )
    
    ## initial values for weights
    weight.vec = rep(1, n.ranker)
    
    ## initial values for sigma2
    sigma2.alpha = tau2.alpha
    sigma2.beta = tau2.beta
  }else{
    
    Z.mat = initial.list$Z.mat
    alpha = initial.list$alpha
    beta = initial.list$beta
    mu = as.vector( alpha + X.mat %*% beta )
    weight.vec = initial.list$weight.vec
    sigma2.alpha = initial.list$sigma2.alpha
    sigma2.beta = initial.list$sigma2.beta
    
  }
  
  ## store initial value
  draw$Z.mat[,,1] = Z.mat
  draw$alpha[,1] = alpha
  draw$beta[,1] = beta
  draw$mu[,1] = mu
  draw$weight.vec[,1] = weight.vec
  
  ## Gibbs iteration
  for(iter in 2:iter.max){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z.mat = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, Z.mat = Z.mat, mu = mu, weight.vec = weight.vec, n.ranker = n.ranker )
    
    # update (alpha, beta) or equivalently mu given Z.mat
    mean.para.update = GibbsUpMuGivenLatentGroup(Z.mat = Z.mat, X.mat = X.mat, weight.vec = weight.vec, sigma2.alpha = sigma2.alpha, sigma2.beta = sigma2.beta, n.ranker = n.ranker, n.item = n.item, p.cov = p.cov, para.expan = para.expan)
    
    ### for check only
    #Z.mat = Z.mat/mean.para.update$theta
    
    alpha = mean.para.update$alpha
    beta = mean.para.update$beta
    mu = as.vector( alpha + X.mat %*% beta )
    
    ### update weight
    weight.vec = GibbsUpWeightGroup(Z.mat = Z.mat, mu = mu, weight.prior.value = weight.prior.value, weight.prior.prob = weight.prior.prob, n.item = n.item, n.ranker = n.ranker)
    
    ### update sigma2
    sigma2.alpha = GibbsUpsigma2(alpha, nu.alpha, tau2.alpha)
    if(p.cov > 0){
      sigma2.beta = GibbsUpsigma2(beta, nu.beta, tau2.beta)
    }
    
    # store value at this iteration
    draw$Z.mat[,,iter] = Z.mat
    draw$alpha[,iter] = alpha
    draw$beta[,iter] = beta
    draw$mu[,iter] = mu
    draw$weight.vec[, iter] = weight.vec
    draw$sigma2.alpha[iter] = sigma2.alpha
    draw$sigma2.beta[iter] = sigma2.beta
    
    # print iteration number
    if(iter %% print.opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(weight.vec))
      # print(c(sigma2.alpha, sigma2.beta))
    }
  }
  return(draw)
}


