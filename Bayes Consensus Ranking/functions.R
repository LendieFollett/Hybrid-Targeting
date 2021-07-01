#--calculate pair.comp_{ij} = 1{Y_i < Y_j} --------------------------------
#' Compute Pairwise Comparison Matrix for Full Ranking List of the Entities
#'
#' Compute the pairwise comparison matrix from the ranking lists of the ranked entities.
#' @param rank.vec A full ranking list containing the ranks of all the \eqn{N} entities. 
#' Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score 
#' of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, 
#' the rank of an entity equals \eqn{N+1} minus its ranked position.
#' ranked higher = lower Z latent score
#' @return An \eqn{N} by \eqn{N} pairwise comparison for all \eqn{N} entities, where the (\eqn{i},\eqn{j}) element equals 1 if \eqn{i} is ranked higher than \eqn{j}, and 0 if \eqn{i} is ranker lower than \eqn{j}. Note that the diagonal elements (\eqn{i},\eqn{i})'s are set to NA.
#' @export
FullRankToPairComp <- function( rank.vec, n = length(rank.vec) ){
  pair.comp <- matrix(NA, n, n)
  for(i in 1:n){
    j = which(rank.vec == i)
    pair.comp[j,  rank.vec > i] = 1 # 
    pair.comp[j,  rank.vec < i] = 0
  }
  return(pair.comp)
}


### Gibbs update Z given beta--------------------------
### pair.comp.ten[[r]] pairwise rank comparison matrix from rth community###
### Z.mat[,j]: latent variable vector for jth ranker ###
### mu: shared mean vector for testing sample N1 x 1 ###
### weight.vec[j]: weight for jth ranker ###
GibbsUpLatentGivenRankGroup <- function(pair.comp.ten, Z, mu, omega_rank , R = ncol(Z) ){
  for(r in 1:R){ #loop over rankers (in simple case, ranker = community)
    rcases <- which(!is.na(Z[,r]) )
    #order from lowest ranked to highest (best well being to worst well being)
    up.order = sort( rowSums( pair.comp.ten[[r]], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix
    #Z needs to be only individuals ranked by ranker r
    Z[rcases,r] = GibbsUpLatentGivenRankInd(pair.comp.ten[[r]], Z[rcases,r],up.order, mu[rcases], weight = omega_rank)
  }
  return(Z)
}


GibbsUpLatentGivenRankInd <- function(pair.comp, Z_sub,up.order, mu_sub, weight){

  for(i in up.order){
    set1 = which( pair.comp[i, ] == 1)
    set0 = which( pair.comp[i, ] != 1)
    
    if(length(set1) > 0){
      upper = min(Z_sub[set1])
    }else{
      upper = Inf
    }
    
    if(length(set0) > 0){
      lower = max(Z_sub[set0])
    }else{
      lower = -Inf
    }
    
    Z_sub[i] = rtruncnorm( 1, lower, upper, mean = mu_sub[i], sd = 1/sqrt(weight) )
  }
  return(Z_sub)
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
GibbsUpMuGivenLatentGroup <- function(X ,
                                      Y ,
                                      omega ,
                                      mu_beta,
                                      rank = FALSE){
  
  #LRF TO ADDRESS: Y_comm might be missing, Y_micro might be missing, ...assuming ranking will be there...
  #                same logic for corresponding x matrices

P <- ncol(X) - 1
  
  #Complete 'data' vector
# u <- c(as.vector(Y_comm), #KxA --> (AK)x1
#             as.vector(Y_micro), #N0xM --> (M*N0)x1
#             as.vector(Z)) #N1xR --> (R*N1)x1
  if(!rank){
    u <- as.vector(Y) 
    
    n <- length(u)
    c <- ncol(Y)
  } else{

    if (any(apply(Y, 1, function(x){sum(!is.na(x))}) > 1)){
      #LRF needs to address: condition for when a person is ranked by multiple sources
    }else{
      u <- apply(Y, 1, function(x){sum(x, na.rm=TRUE)}) #basically take the only non-NA element
      n <- length(u)
      c <- 1
    }
    
  }

  X <-kronecker(rep(1, c), X) #(A*K)xP
  
  Sigma_inv_y<-omega*diag(n)
  
  Sigma_inv_beta <- omega*diag(P+1) #prior covariance matrix
  
  #A<-1x(A*K + M*N0 +R*N1)%*%square(A*K + M*N0 +R*N1)%*%(A*K + M*N0 +R*N1)xP --> 1xP
  pt1 <- u^T%*%Sigma_inv_y%*%X + t(mu_beta)%*%Sigma_inv_beta
  
  pt2 <- t(X)%*%Sigma_inv_y%*%X + Sigma_inv_beta
  
  pt2_inv <- solve(pt2)
  
  alpha_beta <- mvrnorm(1, mu = t(pt1%*%pt2_inv), Sigma = pt2_inv)

  return(alpha_beta)
}



GibbsUpGlobalMuGivenMu<- function(beta_rank = NULL,
                                  beta_comm = NULL,
                                  beta_micro = NULL,
                                  omega_rank = NULL,
                                  omega_comm = NULL,
                                  omega_micro = NULL){
  
  P <- max(length(beta_rank), length(beta_comm), length(beta_micro)) - 1

  if (!is.null(omega_rank)){
    Omega_rank <- diag(P + 1)*1/omega_rank
  } else{
    Omega_rank <- diag(rep(0, P+1))
    beta_rank <- rep(0, P+1)
  }
  
  if (!is.null(omega_comm)){
    Omega_comm <- diag(P + 1)*1/omega_comm
  } else{
    Omega_comm <- diag(rep(0, P+1))
    beta_comm <- rep(0, P+1)
  }
  
  if (!is.null(omega_micro)){
    Omega_micro <- diag(P + 1)*1/omega_micro
  } else{
    Omega_micro <- diag(rep(0, P+1))
    beta_micro <- rep(0, P+1)
  }
  
  
  post_Sigma <- solve(solve(Omega_rank) + solve(Omega_comm) + solve(Omega_micro) + diag(P+1))
  
  post_mu <- (t(beta_rank)%*%solve(Omega_rank) + t(beta_comm)%*%solve(Omega_comm) + t(beta_micro)%*%solve(Omega_micro))%*%post_Sigma
  
  mu_beta <- mvrnorm(1, mu = t(post_mu), Sigma = post_Sigma)
  
  return(mu_beta)
}



### Gibbs update for information quality weights omega_comm, omega_micro, omega_rank---------
#y is either Y_comm (KxA), Y_micro (N0xM), or Z (N1xR)
GibbsUpQualityWeights <- function(y, mu, beta, mu_beta, weight_prior_value = c(0.5, 1, 2), prior_prob = rep(1/length(weight_prior_value), length(weight_prior_value))){
  Col <- ncol(y)
  n.prior.value <- length(weight_prior_value)
  weight_samp <- rep(NA, Col)

  log.post.prob = rep(0, n.prior.value) #re-initialize for next information source   
  for(k in 1:n.prior.value){ #over potential values
    for( col in 1:Col){ #over information source within y
      idx <- which(!is.na(y[,col])) #in the case of omega_rank
      Row <- length(idx)
    log.post.prob[k] <-  log.post.prob[k] +sum(dnorm(y[idx,col], mu[idx], sqrt(1/weight_prior_value[k]), log = TRUE))
    }
    log.post.prob[k] <- log.post.prob[k] + log(prior_prob[k])+ sum(dnorm(beta, mean = mu_beta, sd =sqrt(1/weight_prior_value[k]), log = TRUE ))
  }
#note: w = 1/sigma^2; sigma^2 = 1/w; sigma = 1/sqrt(w)
  log.post.prob = log.post.prob - max(log.post.prob)
  post.prob = exp(log.post.prob)
  
  weight_samp <- weight_prior_value[sample(c(1,2,3),size = 1, prob= post.prob)]
  
  return(weight_samp)
}

#' Bayesian Consensus Targeting WITHOUT random effects for testing subjects (no sigma^2_rank)
#'
#' Implement the Bayesian model for mixed data inputs with covariate information and subjective + objective information with varying qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @import MASS
#' @param X_micro0 An \eqn{N0} by \eqn{P+1} covariate matrix for the \eqn{N0} 'training sample' entities with \eqn{P} covariates. If 1st column isn't 1s, an intercept is added.
#' @param X_micro1 An \eqn{N1} by \eqn{P+1} covariate matrix for the \eqn{N1} 'testing sample' entities with \eqn{P} covariates. If 1st column isn't 1s, an intercept is added.
#' @param X_comm An \eqn{K} by \eqn{P+1} covariate matrix. Usually an aggregate of X_micro0. Must have same columns, in same order, as X_micro0/1.
#' @param X_elite is an optional character string specifying column name of binary elite connection indicator OR the numeric column position. Used for debiasing.
#' @param Y_comm A \eqn{K} by \eqn{A} numeric matrix of community-level response variables. Each column represents a distinct response.
#' @param Y_micro A \eqn{N0} by \eqn{A} numeric matrix of micro-level response variables for 'training sample'. Each column represents a distinct response.
#' @param Tau A \eqn{N1} by \eqn{R} integer matrix, possibly containing NA values, describing ranks of individuals. Each column is distinct 'ranker', each row is distinct individual in 'testing sample'.  
#' @param weight_prior_value A vector for the support of the discrete prior on weight parameter.
#' @param prior_prob_rank A vector for the prior probability mass of the discrete prior on weight parameter for rank. Same for micro and comm.
#' @param iter.keep Number of iterations kept for Gibbs sampler after burn-in.
#' @param iter.burn Number of iterations for burn in (discarded)
#' @return A list containing posterior samples of mu, the shared 'wellness' mean, conditional on the test X_micro1.
#' @export
BCTarget<- function(Tau, X_micro0=NULL, X_micro1=NULL, X_comm=NULL,
                                X_elite = NULL,
                                Y_comm=NULL, Y_micro=NULL,
                                weight_prior_value = c(0.5, 1, 2), 
                                prior_prob_rank = rep(1/length(weight_prior_value), length(weight_prior_value)),
                                prior_prob_micro = rep(1/length(weight_prior_value), length(weight_prior_value)),
                                prior_prob_comm = rep(1/length(weight_prior_value), length(weight_prior_value)),
                                N1 = dim(X_micro1)[1], #how many people in test set
                                R = ncol(Tau), #how many rankers. often will be equal to K
                                iter.keep = 5000,
                                iter.burn = 5000,
                                print.opt = 100,
                                initial.list = NULL){
  #pair.com.ten An \eqn{N1} by \eqn{N1} by \eqn{R} pairwise comparison array for all \eqn{N1} entities and \eqn{R} rankers, 
  #where the (\eqn{i},\eqn{j},\eqn{r}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{r}, 
  #0 if \eqn{i} is ranker lower than \eqn{j}, 
  #and NA if the relation between \eqn{i} and \eqn{j} is missing. 
  #Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{r})'s for all rankers should be set to NA as well.
  #create pair.comp.ten matrix
  pair.comp.ten = list()#array(NA, dim = c(N1, N1, R)) ## get pairwise comparison matrices from the ranking lists
  for(r in 1:R){
    #pair.comp.ten[!is.na(Tau[,r]),!is.na(Tau[,r]),r] = as(FullRankToPairComp( Tau[,r][!is.na(Tau[,r])] ), "dgTMatrix")
    pair.comp.ten[[r]] = FullRankToPairComp( Tau[!is.na(Tau[,r]),r] )
  }
  
  
  
  if (!is.null(X_elite)){
    X_micro1_noelite <-  X_micro1
    X_micro1_noelite[,X_elite] <- 0
  }
  
  # intercept? P includes variables only
  if(all(X_micro1[,1]==1)){
    P <- ncol(X_micro1)-1
  }else{
    X_micro1 <- cbind(1, X_micro1)
    X_micro0 <- cbind(1, X_micro0)
    X_comm <- cbind(1, X_comm)
    P <- ncol(X_micro1)-1
  }
  
  
  
  A <- ncol(Y_comm)
  M <- ncol(Y_micro)
  R <- length(pair.comp.ten)
  K <- nrow(Y_comm)
  N0 <- nrow(Y_micro)
  N1 <- dim(X_micro1)[1]
  
  #will need to add conditional logic here to account for multiple rankers of individuals
  Z.len <- N1
  
  ## store MCMC draws
  draw = list(
   # Z = array(NA, dim = c( iter.keep,N1)),
    mu_beta = array(NA, dim = c(iter.keep,P+1)),
    beta_rank = array(NA, dim = c(iter.keep,P+1)),
    beta_comm = array(NA, dim = c(iter.keep,P+1)),
    beta_micro = array(NA, dim = c(iter.keep,P+1)),
    mu = array(NA, dim = c(iter.keep,N1)),
    mu_noelite = array(NA, dim = c(iter.keep,N1)), #for debiasing
    omega_comm = array(NA, dim = c(iter.keep, 1) ),
    omega_micro = array(NA, dim = c(iter.keep, 1) ),
    omega_rank = array(NA, dim = c(iter.keep, 1) )
  )
  
  ## set initial values for parameters, where given
  if(is.null(initial.list$Z)){
    Z = matrix(NA, nrow = N1, ncol = R)
    for(j in 1:R){
      rcases <- which(!is.na(Tau[,j]))
      nranked <- length(rcases)
      Z[rcases,j][sort( rowSums( pair.comp.ten[[j]], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix]= (c(nranked : 1) - (1+nranked)/2)/sd(c(nranked : 1))
    }
    }else{
      Z <- initial.list$Z
    }
  
  
  if(is.null(initial.list$beta_comm)|is.null(Y_comm)){beta_comm <- rep(0, P+1)}else{  beta_comm <-  initial.list$beta_comm } 
  if(is.null(initial.list$beta_rank)){beta_rank <- rep(0, P+1)}else{  beta_rank <-  initial.list$beta_rank } 
  if(is.null(initial.list$beta_micro)|is.null(Y_micro)){beta_micro <- rep(0, P+1)}else{  beta_micro <-  initial.list$beta_micro } 
  mu_beta <- cbind(beta_comm, beta_rank, beta_micro) %>%apply(1, mean)
  
  ## initial values for weights
  omega_comm = 1
  omega_micro = 1
  omega_rank = 1
  
  ## Gibbs iteration
  for(iter in 1:(iter.burn + iter.keep)){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, Z = Z, mu = X_micro1 %*% beta_rank, omega_rank = omega_rank, R = R )
    
    # ----> update beta_rank
    beta_rank = GibbsUpMuGivenLatentGroup(Y = Z,
                                          X = X_micro1,
                                     omega = omega_rank,
                                     mu_beta = mu_beta,
                                     rank=TRUE)
    
    # ----> update quality weights
    omega_rank <- GibbsUpQualityWeights(y=Z , 
                                        mu=X_micro1 %*% beta_rank, 
                                        beta_rank,  mu_beta,
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob = prior_prob_rank)
    # ----> update beta_comm
    if(!is.null(Y_comm)){
    beta_comm = GibbsUpMuGivenLatentGroup(Y = Y_comm ,
                                          X = X_comm ,
                                          omega = omega_comm,
                                          mu_beta = mu_beta)
    
    # ----> update quality weights
    omega_comm <- GibbsUpQualityWeights(y=Y_comm, 
                                        mu=X_comm %*% beta_comm,
                                        beta_comm, mu_beta,
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob = prior_prob_comm)
    }
    if(is.null(Y_micro)){
    # ----> update beta_micro 
    beta_micro = GibbsUpMuGivenLatentGroup(Y = Y_micro,
                                           X = X_micro0,
                                           omega = omega_micro,
                                           mu_beta = mu_beta)
    
    # ----> update quality weights    
    omega_micro <-GibbsUpQualityWeights(y=Y_micro, 
                                        mu=X_micro0 %*% beta_micro,
                                        beta_micro, mu_beta, 
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob = prior_prob_micro)
    }


    mu_beta <- GibbsUpGlobalMuGivenMu(beta_rank,  beta_comm,  beta_micro,
                           omega_rank, omega_comm, omega_micro )
    

    #LRF TO ADDRESS: this is to be computed with the 'connections' dummy 0'd out
    mu = as.vector( X_micro1 %*% mu_beta  )
    
    if(!is.null(X_elite)){
    mu_noelite = as.vector( X_micro1_noelite %*% mu_beta )
    }else{
      mu_noelite = mu
    }
    
    if(iter > iter.burn){
      j = iter - iter.burn
      # store value at this iteration
      #draw$Z[j,,] = Z
      draw$mu_beta[j,] = mu_beta
      draw$beta_rank[j,] = beta_rank
      draw$beta_comm[j,] = beta_comm
      draw$beta_micro[j,] = beta_micro
      draw$mu[j,] = mu
      draw$mu_noelite[j,] = mu_noelite
      draw$omega_micro[j] = omega_micro
      draw$omega_comm[j] = omega_comm
      draw$omega_rank[j] = omega_rank
    }
    # print iteration number
    if(iter %% print.opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(weight.vec))
      # print(c(sigma2.alpha, sigma2.beta))
    }
  }
  
  
  return(draw)
}

