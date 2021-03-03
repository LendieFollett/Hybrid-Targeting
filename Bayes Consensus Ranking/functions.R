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
### Z.mat[,j]: latent variable vector for jth ranker ###
###X.mat it full, standardized X matrix with training first, then testing
### weight.vec[j]: weight for jth ranker ###
### sigma2.alpha, sigma2.beta: prior parameters for mu = (alpha, beta) ###
### para.expan: whether use parameter expansion ###
GibbsUpMuGivenLatentGroup <- function(Z.mat0, Z.mat1, X.mat = matrix(NA, nrow = nrow(Z.mat0) + nrow(Z.mat1), ncol = 0), weight.vec = rep(1, ncol(Z.mat0)), sigma2.alpha = 2, sigma2.beta = 1, n.ranker = ncol(Z.mat0), n.item =nrow(Z.mat0) + nrow(Z.mat1), p.cov = ncol(X.mat), para.expan = FALSE){
  Z.mat <-Z.mat0
  diagLambda = c( rep(sigma2.alpha, n.item), rep(sigma2.beta, p.cov) )
  V <- cbind( diag(n.item), X.mat )
  
  # Sigma.old = solve( diag(1/diagLambda, nrow = n.item + p.cov) + sum(weight.vec) * t(V) %*% V )
  
  Sigma.inv.eigen = eigen( diag(1/diagLambda, nrow = n.item + p.cov) + sum(weight.vec) * t(V) %*% V )
  Sigma = Sigma.inv.eigen$vectors %*% diag(1/Sigma.inv.eigen$values, nrow = n.item + p.cov, ncol = n.item + p.cov) %*% t(Sigma.inv.eigen$vectors)
  
  lambda = t(V) %*% rowSums( t( t(Z.mat) * weight.vec ) )
  eta = Sigma %*%  lambda
  
  if(para.expan){
    S = sum( colSums(Z.mat^2) * weight.vec ) - as.vector( t(lambda) %*% Sigma %*% lambda )
    theta = as.vector( sqrt( S/rchisq(1, df = n.item * n.ranker) ) )
  }else{
    theta = 1
  }
  
  # alpha.beta = as.vector( rmvnorm(1, mean = eta/theta, sigma = Sigma) )
  alpha.beta = as.vector( eta/theta + Sigma.inv.eigen$vectors %*% diag(1/sqrt(Sigma.inv.eigen$values), nrow = n.item + p.cov, ncol = n.item + p.cov) %*% rnorm(n.item + p.cov) )
  
  
  alpha = alpha.beta[c(1:n.item)]
  beta = alpha.beta[-c(1:n.item)]
  
  ### parameter move
  # alpha = alpha - mean(alpha) + mean( rnorm(n.item, mean = 0, sd = sqrt(sigma2.alpha)) )
  
  return(list(alpha = alpha, beta = beta, theta = theta))
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


