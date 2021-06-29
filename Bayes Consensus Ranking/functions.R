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
  return(as(pair.comp, "dgTMatrix"))
}


### Gibbs update Z given beta--------------------------
### pair.comp.ten[[r]] pairwise rank comparison matrix from rth community###
### Z.mat[,j]: latent variable vector for jth ranker ###
### mu: shared mean vector for testing sample N1 x 1 ###
### weight.vec[j]: weight for jth ranker ###
GibbsUpLatentGivenRankGroup <- function(pair.comp.ten, Z, mu, omega_rank = 1, R = ncol(Z) ){
  for(r in 1:R){ #loop over rankers (in simple case, ranker = community)
    print(r)
    rcases <- which(!is.na(Z[,r]) )
    #order from lowest ranked to highest (best well being to worst well being)
    up.order = sort( rowSums( pair.comp.ten[[r]], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix
    #Z needs to be only individuals ranked by ranker r
    Z[rcases,r] = GibbsUpLatentGivenRankInd2(pair.comp.ten[[r]], Z[rcases,r],up.order, mu, weight = omega_rank)
  }
  return(Z)
}


GibbsUpLatentGivenRankInd <- function(pair.comp, Z,up.order, mu, weight){

  for(i in up.order){
    set1 = which( pair.comp[i, ] == 1)
    set0 = which( pair.comp[i, ] != 1)
    
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
                                      sigma2_beta = 5^2){
  
  #LRF TO ADDRESS: Y_comm might be missing, Y_micro might be missing, ...assuming ranking will be there...
  #                same logic for corresponding x matrices
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
  
  Sigma_inv_diag <- c(rep(omega_comm, each = K), 
                      rep(omega_micro, each = N0),
                      rep(omega_rank, each = N1))
  
  #A<-1x(A*K + M*N0 +R*N1)%*%square(A*K + M*N0 +R*N1)%*%(A*K + M*N0 +R*N1)xP --> 1xP
  pt1 <- u^T%*%(Sigma_inv_diag*Diagonal(length(Sigma_inv_diag)))%*%X
  
  pt2 <- t(X)%*%(Sigma_inv_diag*Diagonal(length(Sigma_inv_diag)))%*%X + diag(P+1)/sigma2_beta
  
  pt2_inv <- solve(pt2)
  
  alpha_beta <- mvrnorm(1, mu = t(pt1%*%pt2_inv), Sigma = pt2_inv)

  return(alpha_beta)
}


# Update random effects for comm, micro, Z
#y is the appropriate matrix
#mu is the corresponding mean
#omega is the error variance of y
#Sigma_gamma is the standard deviation on that source's random effect
#LRF ADDRESS: add sampling steps for sigma_gamma parameters
GibbsUpGammaGivenLatentGroup <- function(y, xbeta, Xr, omega, sigma_gamma = 2.5){
  
  N <- nrow(y) #number of random effects to estimate = number of rows
  #ISSUE: ONLY APPLIES WHEN MULTIPLE SOURCES OF SAME KIND AVAILABLE..
  #Update conditionally if ncol > 1
  Col <- ncol(y)
  
  #Complete 'data' vector
  u <- as.vector(y)
#LRF - 
  Sigma_inv_diag <- c(rep(omega, each =nrow(y)))
  
  Xf <- rep(xbeta, Col)

  pt1 <- (u-Xf)^T%*%(Sigma_inv_diag*Diagonal(length(Sigma_inv_diag)))%*%Xr
  
  pt2 <- t(Xr)%*%(Sigma_inv_diag*Diagonal(length(Sigma_inv_diag)))%*%Xr + diag(N)/(sigma_gamma^2)
  
  pt2_inv <- solve(pt2)
  
  gamma <- mvrnorm(1, mu = t(pt1%*%pt2_inv), Sigma = pt2_inv)
  
  return(gamma)
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
      log( weight.prior.prob[k] ) + (Row/2) * log( weight.prior.value[k] ) - (weight.prior.value[k]/2) * sum( (y[,col] - mu)^2 )
    #log(prior value) - log(sigma) -(1/(2*sigma*sigma))*sum[(y-mu)^2]=
    #log(prior value) + .5log(w) -(w/2)*sum[(y-mu)^2]
  }
#note: w = 1/sigma^2; sigma^2 = 1/w; sigma = 1/sqrt(w)
  log.post.prob = log.post.prob - max(log.post.prob)
  post.prob = exp(log.post.prob)
  
  weight_samp[col] <- weight.prior.value[ as.vector( rmultinom(1, 1, prob = post.prob) ) == 1 ]
  
}

  return(weight_samp)
}

### Gibbs update for variances on random effects
### Gibbs update for sigma2, given prior sigma2 ~ Scale-Inv-chi2(nu, tau2) and data iid ~ N(0, sigma2)
GibbsUpsigma2 <- function(x, nu, tau2){
  if(nu < Inf){
    n.x = length(x)
    
    nu.post = nu + n.x
    tau2.post = ( nu * tau2 + sum(x^2) )/(nu + n.x)
    
    sigma2 = tau2.post * nu.post/rchisq(1, df = nu.post)
  }else{
    sigma2 = tau2
  }
  
  return(sigma2)
}


#' Bayesian Consensus Targeting WITHOUT random effects for testing subjects (no sigma^2_rank)
#'
#' Implement the Bayesian model for mixed data inputs with covariate information and subjective + objective information with varying qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @import MASS
#' @param X_micro0 An \eqn{N0} by \eqn{P+1} covariate matrix for the \eqn{N0} entities with \eqn{P} covariates. Assumes 1st col is 1's for intercept.
#' @param X_micro1 An \eqn{N1} by \eqn{P+1} covariate matrix for the \eqn{N1} entities with \eqn{P} covariates. Assumes 1st col is 1's for intercept.
#' @param X_comm An \eqn{K} by \eqn{P+1} covariate matrix. An aggregate of X_micro0 and X_micro 1. 
#' @param elite is an optional character string specifying column name of binary elite connection indicator OR the column position. used for debiasing
#' @param sigma_beta Currently given/fixed. Prior variance on beta.
#' @param weight.prior.value A vector for the support of the discrete prior on weight parameter.
#' @param weight.prior.prob A vector for the prior probability mass of the discrete prior on weight parameter.
#' @param iter.keep Number of iterations kept for Gibbs sampler after burn-in.
#' @param iter.burn Number of iterations for burn in (discarded)
#' @return A list containing posterior samples of mu, the shared 'wellness' mean, conditional on the test X_micro1.
#' @export
BCTarget<- function(Tau, X_micro0, X_micro1, X_comm,
                                X_elite = NULL,
                                Y_comm, Y_micro,
                                sigma2_beta = 5^2,
                                weight.prior.value = c(0.5, 1, 2), 
                                weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                                N1 = dim(X_micro1)[1], #how many people in test set
                                R = length(pair.comp.ten), #how many rankers. often will be equal to K
                                iter.keep = 5000,
                                iter.burn = 5000,
                                para.expan = TRUE, print.opt = 100,
                                initial.list = NULL){
  #pair.com.ten An \eqn{N1} by \eqn{N1} by \eqn{R} pairwise comparison array for all \eqn{N1} entities and \eqn{R} rankers, 
  #where the (\eqn{i},\eqn{j},\eqn{r}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{r}, 
  #0 if \eqn{i} is ranker lower than \eqn{j}, 
  #and NA if the relation between \eqn{i} and \eqn{j} is missing. 
  #Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{r})'s for all rankers should be set to NA as well.
  #create pair.comp.ten matrix
  pair.comp.ten = list()#array(NA, dim = c(N1, N1, R)) ## get pairwise comparison matrices from the ranking lists
  for(r in 1:R){
    print(r)
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
    P <- ncol(X_micro1)
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
    Z = array(NA, dim = c( iter.keep,N1)),
    beta = array(NA, dim = c(iter.keep,P+1)),
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
  
  if(is.null(initial.list$beta)){beta <- rep(0, P+1)}else{  beta <-  initial.list$beta} 
  if(is.null(initial.list$gamma_comm)){gamma_comm <- rep(0, K)}else{  gamma_comm <-  initial.list$gamma_comm} 
  if(is.null(initial.list$gamma_micro)){gamma_micro <- rep(0, N0)}else{  gamma_micro <-  initial.list$gamma_micro} 
  
  
  ## initial values for alpha, beta and thus mu
  mu <- as.vector(X_micro1 %*% beta )
  
  ## initial values for weights
  omega_comm = rep(1, A) 
  omega_micro = rep(1, M) 
  omega_rank = rep(1, R)
  
  ## initial values for sigma2
  

  ## Gibbs iteration
  for(iter in 1:(iter.burn + iter.keep)){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, Z = Z, mu = mu, omega_rank = omega_rank, R = R )
    
    # update beta (includes intercept) or equivalently mu given Z
    beta_rank = GibbsUpMuGivenLatentGroup(Y = ,
                                          X = ,
                                     omega = ,
                                     sigma2_beta = 5^2)
    
    beta_comm = GibbsUpMuGivenLatentGroup(Y = ,
                                          X = ,
                                          omega = ,
                                          sigma2_beta = 5^2)
    
    beta_micro = GibbsUpMuGivenLatentGroup(Y = ,
                                           X = ,
                                           omega = ,
                                          sigma2_beta = 5^2)
    
  
    # update quality weights
    # LRF - INCLUDE RANDOM EFFECTS
    omega_comm <- GibbsUpQualityWeights(y=Y_comm-kronecker(t(rep(1, A)), gamma_comm), mu=X_comm %*% beta, weight.prior.value = c(0.5, 1, 2 ))
    omega_micro <-GibbsUpQualityWeights(y=Y_micro-kronecker(t(rep(1, M)), gamma_micro), mu=X_micro0 %*% beta, weight.prior.value = c(0.5, 1, 2 ))
    omega_rank <- GibbsUpQualityWeights(y=Z , mu=X_micro1 %*% beta, weight.prior.value = c(0.5, 1, 2 ))


    #LRF TO ADDRESS: this is to be computed with the 'connections' dummy 0'd out
    mu = as.vector( X_micro1 %*% beta  )
    
    if(!is.null(X_elite)){
    mu_noelite = as.vector( X_micro1_noelite %*% beta )
    }else{
      mu_noelite = mu
    }
    
    if(iter > iter.burn){
      j = iter - iter.burn
      # store value at this iteration
      draw$Z[j,,] = Z
      draw$beta[j,] = beta
      draw$mu[j,] = mu
      draw$mu_noelite[j,] = mu_noelite
      draw$omega_micro[j,] = omega_micro
      draw$omega_comm[j,] = omega_comm
      draw$omega_rank[j,] = omega_rank
      draw$gamma_comm[j,] = gamma_comm
      draw$gamma_micro[j,] = gamma_micro
      draw$sigma2_comm[j] = sigma2_comm
      draw$sigma2_micro[j] = sigma2_micro
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

#--Bayesian Consensus Targeting with random effects--------------------------------

#' Bayesian Consensus Targeting with random effects
#'
#' Implement the Bayesian model for mixed data inputs with covariate information and subjective + objective information with varying qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @import MASS
#' @param pair.com.ten An \eqn{N1} by \eqn{N1} by \eqn{R} pairwise comparison array for all \eqn{N1} entities and \eqn{R} rankers, 
#'        where the (\eqn{i},\eqn{j},\eqn{r}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{r}, 
#'        0 if \eqn{i} is ranker lower than \eqn{j}, 
#'        and NA if the relation between \eqn{i} and \eqn{j} is missing. 
#'        Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{r})'s for all rankers should be set to NA as well.
#' @param X_micro0 An \eqn{N0} by \eqn{P+1} covariate matrix for the \eqn{N0} entities with \eqn{P} covariates. Assumes 1st col is 1's for intercept.
#' @param X_micro1 An \eqn{N1} by \eqn{P+1} covariate matrix for the \eqn{N1} entities with \eqn{P} covariates. Assumes 1st col is 1's for intercept.
#' @param X_comm An \eqn{K} by \eqn{P+1} covariate matrix. An aggregate of X_micro0 and X_micro 1. 
#' @param elite is an optional character string specifying column name of binary elite connection indicator OR the column position. used for debiasing
#' @param sigma_beta Currently given/fixed. Prior variance on beta.
#' @param weight.prior.value A vector for the support of the discrete prior on weight parameter.
#' @param weight.prior.prob A vector for the prior probability mass of the discrete prior on weight parameter.
#' @param iter.keep Number of iterations kept for Gibbs sampler after burn-in.
#' @param iter.burn Number of iterations for burn in (discarded)
#' @return A list containing posterior samples of mu, the shared 'wellness' mean, conditional on the test X_micro1.
#' @export
BCRTarget_inc_gamma <- function(pair.comp.ten, X_micro0, X_micro1, X_comm,
                     X_elite = NULL,
                               Y_comm, Y_micro,
                               sigma2_beta = 5^2,
                               weight.prior.value = c(0.5, 1, 2), 
                               weight.prior.prob = rep(1/length(weight.prior.value), length(weight.prior.value)),
                               N1 = dim(pair.comp.ten)[1], 
                               R = dim(pair.comp.ten)[3], 
                               iter.keep = 5000, 
                               iter.burn = 5000,
                               para.expan = TRUE, print.opt = 100,
                               initial.list = NULL){

  if (!is.null(X_elite)){
    X_micro1_noelite <-  X_micro1[,X_elite] <- 0
  }
  
  if(all(X_micro1[,1]==1)){
    P <- ncol(X_micro1)-1
  }else{
    P <- ncol(X_micro1)
  }

  
  A <- ncol(Y_comm)
  M <- ncol(Y_micro)
  R <- dim(pair.comp.ten)[3]
  K <- nrow(Y_comm)
  N0 <- nrow(Y_micro)
  N1 <- dim(pair.comp.ten)[1]
  
  #construct random effect matrices
  Xr_micro <-  kronecker(rep(1, M),Diagonal(N0)) #for training micro set
  Xr_comm<-  kronecker(rep(1, A),Diagonal(K)) #for community agg set  
  Xr_rank<-  kronecker(rep(1, R),Diagonal(N1)) #for rank set  
  

  
  ## store MCMC draws
  draw = list(
    Z = array(NA, dim = c( iter.keep,N1, R)),
    beta = array(NA, dim = c(iter.keep,P+1)),
    mu = array(NA, dim = c(iter.keep,N1)),
    mu_noelite = array(NA, dim = c(iter.keep,N1)), #for debiasing
    omega_comm = array(NA, dim = c(iter.keep, A) ),
    omega_micro = array(NA, dim = c(iter.keep, M) ),
    omega_rank = array(NA, dim = c(iter.keep, R) ),
    gamma_rank = array(NA, dim = c(iter.keep, N1) ),
    gamma_comm = array(NA, dim = c(iter.keep, K) ),
    gamma_micro = array(NA, dim = c(iter.keep, N0) ),
    sigma2_rank = rep(NA, iter.keep),
    sigma2_comm = rep(NA, iter.keep),
    sigma2_micro = rep(NA, iter.keep)
  )
  
  ## set initial values for parameters, where given
  if(is.null(initial.list$Z)){
    Z = matrix(NA, nrow = N1, ncol = R)
    for(j in 1:R){
      Z[sort( rowSums( pair.comp.ten[,,j], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix, j] = (c(N1 : 1) - (1+N1)/2)/sd(c(N1 : 1))
    }}else{
      Z <- initial.list$Z
    }
    
  if(is.null(initial.list$beta)){beta <- rep(0, P+1)}else{  beta <-  initial.list$beta} 
  if(is.null(initial.list$gamma_rank)){gamma_rank <- rep(0, N1)}else{  gamma_rank <-  initial.list$gamma_rank} 
  if(is.null(initial.list$gamma_comm)){gamma_comm <- rep(0, K)}else{  gamma_comm <-  initial.list$gamma_comm} 
  if(is.null(initial.list$gamma_micro)){gamma_micro <- rep(0, N0)}else{  gamma_micro <-  initial.list$gamma_micro} 

    
    ## initial values for alpha, beta and thus mu
    mu <- as.vector(X_micro1 %*% beta )
    
    ## initial values for weights
    omega_comm = rep(1, A) 
    omega_micro = rep(1, M) 
    omega_rank = rep(1, R)
    
    sigma2_comm= 1
    sigma2_micro= 1
    sigma2_rank= 1
    ## initial values for sigma2

  
  ## Gibbs iteration
  for(iter in 1:(iter.burn + iter.keep)){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, Z = Z, mu = mu, omega_rank = omega_rank, R = R )

    # update beta or equivalently mu given Z, random effects
    #random effects subtracted off response matrices prior to posterior computations
    beta = GibbsUpMuGivenLatentGroup(Z = Z - kronecker(t(rep(1, R)), gamma_rank), 
                                     Y_comm = Y_comm-kronecker(t(rep(1, A)), gamma_comm), 
                                     Y_micro = Y_micro-kronecker(t(rep(1, M)), gamma_micro),
                                                 X_comm = X_comm, X_micro0 = X_micro0, X_micro1 = X_micro1,
                                                 omega_comm=omega_comm, omega_micro = omega_micro, omega_rank = omega_rank,
                                                 sigma2_beta = 5^2)
    
    
    #update gamma, the random effects
    gamma_rank <- GibbsUpGammaGivenLatentGroup(Z,      X_micro1 %*% beta, Xr_rank, omega_rank, sigma_gamma = sigma2_rank)
    gamma_comm <- GibbsUpGammaGivenLatentGroup(Y_comm, X_comm %*% beta,   Xr_comm, omega_comm, sigma_gamma = sigma2_comm)
    gamma_micro <-GibbsUpGammaGivenLatentGroup(Y_micro,X_micro0 %*% beta, Xr_micro,omega_micro, sigma_gamma = sigma2_micro)
      
    # update quality weights
    # LRF - INCLUDE RANDOM EFFECTS
    omega_comm <- GibbsUpQualityWeights(y=Y_comm-kronecker(t(rep(1, A)), gamma_comm), mu=X_comm %*% beta, weight.prior.value = c(0.5, 1, 2 ))
    omega_micro <-GibbsUpQualityWeights(y=Y_micro-kronecker(t(rep(1, M)), gamma_micro), mu=X_micro0 %*% beta, weight.prior.value = c(0.5, 1, 2 ))
    omega_rank <- GibbsUpQualityWeights(y=Z - kronecker(t(rep(1, R)),gamma_rank), mu=X_micro1 %*% beta, weight.prior.value = c(0.5, 1, 2 ))
    
    # update random effect variances
    sigma2_comm <- GibbsUpsigma2(gamma_comm, nu=3, tau2=25)
    sigma2_micro <- GibbsUpsigma2(gamma_micro, nu=3, tau2=25)
    sigma2_rank <- GibbsUpsigma2(gamma_rank, nu=3, tau2=25)
    
    #LRF TO ADDRESS: this is to be computed with the 'connections' dummy 0'd out
    mu = as.vector( X_micro1 %*% beta + Xr_rank[1:N1,]%*%gamma_rank )
    
    mu_noelite = as.vector( X_micro1_noelite %*% beta + Xr_rank[1:N1,]%*%gamma_rank )

    if(iter > iter.burn){
      j = iter - iter.burn
    # store value at this iteration
    draw$Z[j,,] = Z
    draw$beta[j,] = beta
    draw$mu[j,] = mu
    draw$mu_noelite[j,] = mu_noelite
    draw$omega_micro[j,] = omega_micro
    draw$omega_comm[j,] = omega_comm
    draw$omega_rank[j,] = omega_rank
    draw$gamma_rank[j,] = gamma_rank
    draw$gamma_comm[j,] = gamma_comm
    draw$gamma_micro[j,] = gamma_micro
    draw$sigma2_comm[j] = sigma2_comm
    draw$sigma2_micro[j] = sigma2_micro
    draw$sigma2_rank[j] = sigma2_rank
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
