library(DirichletReg)
#--calculate pair.comp_{ij} = 1{Y_i < Y_j} --------------------------------
#' Compute Pairwise Comparison Matrix for Full Ranking List of the Entities
#'
#' Compute the pairwise comparison matrix from the ranking lists of the ranked entities.
#' @param rank.vec A full ranking list containing the ranks of all the entities. 
#' Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score 
#' of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, 
#' the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @return An \eqn{N} by \eqn{N} pairwise comparison for all \eqn{N} entities, 
#' where the (\eqn{i},\eqn{j}) element equals 1 if \eqn{i} is ranked higher than \eqn{j}, 
#' and 0 if \eqn{i} is ranker lower than \eqn{j}. Note that the diagonal elements (\eqn{i},\eqn{i})'s are set to NA.
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
    Z[rcases,r] = GibbsUpLatentGivenRankInd(pair.comp.ten[[r]], Z[rcases,r],up.order, mu[rcases], weight = omega_rank[r])
  }
  return(Z)
}


GibbsUpLatentGivenRankInd <- function(pair.comp, Z_sub,up.order, mu_sub, weight){

  for(i in up.order){
    set1 = which( pair.comp[i, ] == 1) #who has a LARGER rank than person i
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



# Gibbs update for beta (includes intercept) given Z, y_micro---------------------

### Gibbs update for the shared mean mu ###
### Z.mat is a N1 x R matrix ###
### Z.mat[,j]: latent variable vector for jth ranker j = 1, ..., R ###
### X_micro is a N0xP matrix (training dataset, on the household level) ###
### X_micro is a N1xP matrix (testing dataset, on the household level) ###
### weight.vec: (A + M + R)x1 vector of weights omega_rank, omega_micro###
### omega_rank[r] = weight of r^th ranker###
### sigma2.alpha, sigma2.beta: prior parameters for mu = (alpha, beta) ###
GibbsUpMuGivenLatentGroup <- function(X ,
                                      Y ,
                                      omega ,
                                      mu_beta,
                                      rank = FALSE,
                                      con){
  
  #LRF TO ADDRESS: Y_comm might be missing, Y_micro might be missing, ...assuming ranking will be there...
  #                same logic for corresponding x matrices

P <- ncol(X) - 1
  

  if(!rank){ #if it's micro
    u <- as.vector(Y)
    
    n <- length(u)
    c <- ncol(Y)
    Sigma_inv_y<-omega*diag(n) 
  } else{ #if it's rank

    if (any(apply(Y, 1, function(x){sum(!is.na(x))}) > 1)){
      #LRF needs to address: condition for when a person is ranked by multiple sources
    }else{
      u <- apply(Y, 1, function(x){sum(x, na.rm=TRUE)}) #basically take the only non-NA element
      n <- length(u)
      c <- 1
    }
    Sigma_inv_y<-rep(omega, times = apply(Y, 2, function(x){sum(!is.na(x))}))*diag(n) 
  }

  if(!is.null(c)){
  X <-kronecker(rep(1, c), X) #(A*K)xP
  } 



if(!rank){  
  Sigma_inv_beta <- diag(c(1/2.5^2, (1/con^2)*rep(1, P))) #prior variance on beta is con^2 beta = N(mu_beta, con^2) (removed omega)
  
  #A<-1x(A*K + M*N0 +R*N1)%*%square(A*K + M*N0 +R*N1)%*%(A*K + M*N0 +R*N1)xP --> 1xP

  pt1 <- u^T%*%Sigma_inv_y%*%X + t(c(0,mu_beta))%*%Sigma_inv_beta
  
  pt2 <- t(X)%*%Sigma_inv_y%*%X + Sigma_inv_beta
  
  pt2_inv <- solve(pt2)
  
  alpha_beta <- mvrnorm(1, mu = t(pt1%*%pt2_inv), Sigma = pt2_inv)
  
}else{
  Sigma_inv_beta <- diag((1/con^2)*rep(1, P)) #prior covariance matrix (removed omega)
  
  #A<-1x(A*K + M*N0 +R*N1)%*%square(A*K + M*N0 +R*N1)%*%(A*K + M*N0 +R*N1)xP --> 1xP
  
  pt1 <- u^T%*%Sigma_inv_y%*%X[,-1] + t(mu_beta)%*%Sigma_inv_beta
  
  pt2 <- t(X[,-1])%*%Sigma_inv_y%*%X[,-1] + Sigma_inv_beta
  
  pt2_inv <- solve(pt2)
  
  alpha_beta <-c(0,mvrnorm(1, mu = t(pt1%*%pt2_inv), Sigma = pt2_inv)) #intercept is 0
  
}
  
  return(alpha_beta)
}



GibbsUpGlobalMuGivenMu<- function(beta_rank = NULL,
                                  beta_micro = NULL,
                                  omega_rank = NULL,
                                  omega_micro = NULL,
                                  con){
  
  P <- max(length(beta_rank), length(beta_micro)) - 1

  Omega_rank <-  diag(rep(1, P))*con^2#diag(P + 1)*1/omega_rank  (removed omega)
  Omega_micro <- diag(rep(1, P))*con^2#diag(P + 1)*1/omega_micro  (removed omega)
  
  post_Sigma <- solve(solve(Omega_rank) + solve(Omega_micro) + diag(P)/1^2) #prior sd on mu_beta = 1
  
  post_mu <- (t(beta_rank[-1])%*%solve(Omega_rank) + t(beta_micro[-1])%*%solve(Omega_micro))%*%post_Sigma
  
  mu_beta <- mvrnorm(1, mu = t(post_mu), Sigma = post_Sigma)
  
  return(mu_beta)
}



### Gibbs update for information quality weights omega_micro, omega_rank---------
#y is either Y_comm (KxA), Y_micro (N0xM), or Z (N1xR)
#prior_prob will be current iteration's dirichlet weights in the case of updating omega_rank
GibbsUpQualityWeights <- function(y, mu, beta, mu_beta,con, weight_prior_value = c(0.5, 1, 2), prior_prob = rep(1/length(weight_prior_value), length(weight_prior_value)),rank=FALSE){
  
  if(!rank){ #if it's not omega_rank
  Col <- ncol(y)
  if(is.null(Col)){Col <- 1}
  n.prior.value <- length(weight_prior_value)
  weight_samp <- rep(NA, Col)

  log.post.prob = rep(0, n.prior.value) #re-initialize for next information source   
  for(k in 1:n.prior.value){ #over potential values
    for( col in 1:Col){ #over information source within y
      idx <- which(!is.na(y[,col])) #in the case of omega_rank
      log.post.prob[k] <-  log.post.prob[k] +sum(dnorm(y[idx,col], mu[idx], sqrt(1/weight_prior_value[k]), log = TRUE))
    }
    log.post.prob[k] <- log.post.prob[k] + log(prior_prob[k])#+ sum(dnorm(beta[-1], mean = mu_beta, sd =sqrt(con/weight_prior_value[k]), log = TRUE ))
  }
#note: w = 1/sigma^2; sigma^2 = 1/w; sigma = 1/sqrt(w)
  log.post.prob = log.post.prob - max(log.post.prob)
  post.prob = exp(log.post.prob)
  
  weight_samp <- weight_prior_value[sample(c(1,2,3),size = 1, prob= post.prob)]
  
  return(weight_samp)
  }else{ #if it is omega_rank
    Col <- ncol(y) #need to update an omega for every oclumn

    n.prior.value <- length(weight_prior_value)
    weight_samp <- rep(NA, Col)
    
  for( col in 1:Col){ #over information source within y
    log.post.prob = rep(0, n.prior.value) #re-initialize for next information source   
    idx <- which(!is.na(y[,col])) #in the case of omega_rank
    for(k in 1:n.prior.value){ #over potential values
        log.post.prob[k] <-  log.post.prob[k] +sum(dnorm(y[idx,col], mu[idx], sqrt(1/weight_prior_value[k]), log = TRUE))
        log.post.prob[k] <- log.post.prob[k] + log(prior_prob[k])#+ sum(dnorm(beta[-1], mean = mu_beta, sd =sqrt(con/weight_prior_value[k]), log = TRUE ))
    }
    
    log.post.prob = log.post.prob - max(log.post.prob)
    post.prob = exp(log.post.prob)
    
    weight_samp[col] <- weight_prior_value[sample(c(1,2,3),size = 1, prob= post.prob)]
    
    }
    #note: w = 1/sigma^2; sigma^2 = 1/w; sigma = 1/sqrt(w)
    return(weight_samp)
  }
}


GibbsUpConstant <- function(beta_rank, beta_micro, mu_beta, omega_rank, omega_micro,con_old){
  
  if(!is.null(beta_micro)){
    y <- c(beta_rank[-1], beta_micro[-1])
    mu <- rep(mu_beta, 2)
  }else{
    y <- c(beta_rank[-1])
    mu <- mu_beta
  }
  
  con_prop <- con_old + rnorm(1, 0, .05) #random walk on standard deviation scale
  while(con_prop < 0){
    con_prop <- con_old + rnorm(1, 0, .05)
  }
  
 lik_old <-  dnorm(y, mu, con_old, log=TRUE) %>%sum
 lik_prop <- dnorm(y, mu, con_prop, log=TRUE) %>%sum 
  
  prior_old <- dnorm(con_old, 0, 2.5, log=TRUE)#prior constant ~ N^+(0,2.5)
  prior_prop<- dnorm(con_prop, 0, 2.5, log = TRUE)
  
  alpha <- lik_prop + prior_prop - lik_old - prior_old
  
  if(alpha > log(runif(1,0,1))){
    con <- con_prop
  }else {
    con <- con_old
  }
  
  return(con)
}

GibbsUpLambda <- function(omega, prior_prob, prior_weights){
  
  sum1 <- sum(omega == prior_weights[1])
  sum2 <- sum(omega == prior_weights[2])
  sum3 <- sum(omega == prior_weights[3])
  
 return(rdirichlet(n=1,alpha = c(sum1 + prior_prob[1],sum2 + prior_prob[2],sum3 + prior_prob[3])))
  
}


#' Bayesian Consensus Targeting WITHOUT random effects for testing subjects (no sigma^2_rank)
#'
#' Implement the Bayesian model for mixed data inputs with covariate information and subjective + objective information with varying qualities or weights.
#' @import truncnorm
#' @import mvtnorm
#' @import MASS
#' @param X_PMT An \eqn{N0} by \eqn{P+1} covariate matrix for the \eqn{N0} 'training sample' entities with \eqn{P} covariates. If 1st column isn't 1s, an intercept is added.
#' @param X_CBT An \eqn{N1} by \eqn{P+1} covariate matrix for the \eqn{N1} 'testing sample' entities with \eqn{P} covariates. If 1st column isn't 1s, an intercept is added.
#' @param X_program An \eqn{N} by \eqn{P+1} covariate matrix for the program area. Will be used for predictions. 
#' @param X_elite is an optional character string specifying column name of binary elite connection indicator OR the numeric column position. Used for debiasing.
#' @param Y_micro A \eqn{N0} by \eqn{A} numeric matrix of micro-level response variables for 'training sample'. Each column represents a distinct response.
#' @param Tau A \eqn{N1} by \eqn{R} integer matrix, possibly containing NA values, describing ranks of individuals. Each column is distinct 'ranker', each row is distinct individual in 'testing sample'.  
#' @param weight_prior_value A vector for the support of the discrete prior on weight parameter.
#' @param prior_prob_rank A vector for the prior probability mass of the discrete prior on weight parameter for rank. Same for micro and comm.
#' @param iter_keep Number of iterations kept for Gibbs sampler after burn-in.
#' @param iter_burn Number of iterations for burn in (discarded)
#' @param print_opt Frequency of printing MCMC sampling progress.
#' @return A list containing posterior samples of mu, the shared 'wellness' mean, conditional on the test X_CBT.
#' @export
HybridTarget<- function(Tau, X_PMT=NULL, X_CBT=NULL, X_program=NULL,
                                X_elite = NULL,
                                Y_micro=NULL,
                                weight_prior_value = c(0.5, 1, 2), 
                                prior_prob_rank = rep(1/length(weight_prior_value), length(weight_prior_value)),
                                prior_prob_micro = rep(1/length(weight_prior_value), length(weight_prior_value)),
                                N1 = dim(X_CBT)[1], #how many people in test set
                                R = ncol(Tau), #how many rankers. often will be equal to K
                                iter_keep = 5000,
                                iter_burn = 5000,
                                print_opt = 100,
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
    X_program_noelite <-  X_program
    X_program_noelite[,X_elite] <- 0
  }
  
  # intercept? P includes variables only
  if(all(X_CBT[,1]==1)){
    P <- ncol(X_CBT)-1
  }else{
    X_CBT <- cbind(1, X_CBT)
    X_PMT <- cbind(1, X_PMT)
    P <- ncol(X_CBT)-1
  }
  
  
  M <- ncol(Y_micro)
  R <- length(pair.comp.ten)
  N0 <- nrow(Y_micro)
  N1 <- dim(X_CBT)[1]
  
  ## store MCMC draws
  draw = list(
    lambda = array(NA, dim = c(iter_keep, length(prior_prob_rank))),
    Z = array(NA, dim = c( iter_keep,N1)),
    mu_beta = array(NA, dim = c(iter_keep,P)),
    beta_rank = array(NA, dim = c(iter_keep,P+1)),
    beta_micro = array(NA, dim = c(iter_keep,P+1)),
    mu = array(NA, dim = c(iter_keep,nrow(X_program))),
    mu_noelite = array(NA, dim = c(iter_keep,nrow(X_program))), #for debiasing
    omega_micro = array(NA, dim = c(iter_keep, 1) ),
    omega_rank = array(NA, dim = c(iter_keep, R) ),
    con = array(NA, dim = c(iter_keep, 1) )
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
  
  
  if(is.null(initial.list$beta_rank)){beta_rank <- rep(0, P+1)}else{  beta_rank <-  initial.list$beta_rank } 
  if(is.null(initial.list$beta_micro)|is.null(Y_micro)){beta_micro <- rep(0, P+1)}else{  beta_micro <-  initial.list$beta_micro } 
  if(is.null(initial.list$con)){con <- .5}else{ con <- initial.list$con  } 
  mu_beta <- cbind(beta_rank[-1], beta_micro[-1]) %>%apply(1, mean)
  
  ## initial values for weights
  omega_micro = 1
  omega_rank = rep(1, R)
  ## Gibbs iteration
  for(iter in 1:(iter_burn + iter_keep)){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, 
                                    Z = Z, 
                                    mu = X_CBT %*% beta_rank, 
                                    omega_rank = omega_rank, 
                                    R = R )
    #Z <- (Z -mean(Z, na.rm = TRUE))/sd(Z, na.rm = TRUE)
    
    # ----> update beta_rank
    beta_rank = GibbsUpMuGivenLatentGroup(Y = Z,
                                          X = X_CBT,
                                     omega = omega_rank,
                                     mu_beta = mu_beta,
                                     con = con,
                                     rank=TRUE)
    
    # ----> update quality weights
    lambda <- GibbsUpLambda(omega_rank, prior_prob =prior_prob_rank, prior_weights=c(0.5, 1, 2 ))
    omega_rank <- GibbsUpQualityWeights(y=Z , 
                                        mu=X_CBT %*% beta_rank, 
                                        beta_rank,  mu_beta,
                                        con = con,
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob = lambda,
                                        rank=TRUE)


    # ----> update beta_micro 
    beta_micro = GibbsUpMuGivenLatentGroup(Y = Y_micro,
                                           X = X_PMT,
                                           omega = omega_micro,
                                           mu_beta = mu_beta, 
                                           con = con)
    
    # ----> update quality weights    
    omega_micro <-GibbsUpQualityWeights(y=Y_micro, 
                                        mu=X_PMT %*% beta_micro,
                                        beta_micro, mu_beta, 
                                        con = con,
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob = prior_prob_micro)
    


    mu_beta <- GibbsUpGlobalMuGivenMu(beta_rank,  beta_micro,
                           omega_rank, omega_micro ,con)
    
    
    con <- GibbsUpConstant(beta_rank, beta_micro, mu_beta, omega_rank, omega_micro,con)

    #LRF TO ADDRESS: this is to be computed with the 'connections' dummy 0'd out
    mu = as.vector( X_program %*% c(omega_micro[1],mu_beta ))
    if(!is.null(X_elite)){
    mu_noelite = as.vector( X_program_noelite %*% c(omega_micro[1],mu_beta ) )
    }else{
      mu_noelite = mu
    }
    
    if(iter > iter_burn){
      j = iter - iter_burn
      # store value at this iteration
      draw$Z[j,] = apply(Z, 1, function(x){sum(x,na.rm=TRUE)})
      draw$mu_beta[j,] = mu_beta
      draw$lambda[j,] = lambda
      draw$beta_rank[j,] = beta_rank
      draw$beta_micro[j,] = beta_micro
      draw$mu[j,] = mu
      draw$mu_noelite[j,] = mu_noelite
      draw$omega_micro[j] = omega_micro
      draw$omega_rank[j,] = omega_rank
      draw$con[j] = con
    }
    # print iteration number
    if(iter %% print_opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(weight.vec))
      # print(c(sigma2.alpha, sigma2.beta))
    }
  }
  
  
  return(draw)
}



#CBT RANKING ONLY
CBTarget<- function(Tau, X_CBT=NULL, X_program=NULL,
                    X_elite = NULL,
                    weight_prior_value = c(0.5, 1, 2), 
                    prior_prob_rank = rep(1/length(weight_prior_value), length(weight_prior_value)),
                    N1 = dim(X_CBT)[1], #how many people in test set
                    R = ncol(Tau), #how many rankers. often will be equal to K
                    iter_keep = 5000,
                    iter_burn = 5000,
                    print_opt = 100,
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
    X_program_noelite <-  X_program
    X_program_noelite[,X_elite] <- 0
  }
  
  # intercept? P includes variables only
  if(all(X_CBT[,1]==1)){
    P <- ncol(X_CBT)-1
  }else{
    X_CBT <- cbind(1, X_CBT)
    P <- ncol(X_CBT)-1
  }
  
  R <- length(pair.comp.ten)
  N1 <- dim(X_CBT)[1]
  
  ## store MCMC draws
  draw = list(
    Z = array(NA, dim = c( iter_keep,N1)),
    beta_rank = array(NA, dim = c(iter_keep,P+1)),
    mu = array(NA, dim = c(iter_keep,nrow(X_program))),
    mu_noelite = array(NA, dim = c(iter_keep,nrow(X_program))), #for debiasing
    omega_rank = array(NA, dim = c(iter_keep, R) ),
    con = array(NA, dim = c(iter_keep, 1) )
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
  
  
  if(is.null(initial.list$beta_rank)){beta_rank <- rep(0, P+1)}else{  beta_rank <-  initial.list$beta_rank } 
  if(is.null(initial.list$con)){con <- .5}else{ con <- initial.list$con  } 
  
  ## initial values for weights
  omega_rank = rep(1, R)
  lambda <- prior_prob_rank
  #prior mean on beta_rank
  mu_beta = rep(0, length(beta_rank)-1)
  ## Gibbs iteration
  for(iter in 1:(iter_burn + iter_keep)){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, 
                                    Z = Z, 
                                    mu = X_CBT %*% beta_rank, 
                                    omega_rank = omega_rank, 
                                    R = R )
    #Z <- (Z -mean(Z, na.rm = TRUE))/sd(Z, na.rm = TRUE)
    
    # ----> update beta_rank
    beta_rank = GibbsUpMuGivenLatentGroup(Y = Z,
                                          X = X_CBT,
                                          omega = omega_rank,
                                          mu_beta = mu_beta,
                                          con = con,
                                          rank=TRUE)
    
    # ----> update quality weights
    omega_rank <- GibbsUpQualityWeights(y=Z , 
                                        mu=X_CBT %*% beta_rank, 
                                        beta_rank,  mu_beta,
                                        con = con,
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob = lambda,rank=TRUE)
    lambda <- GibbsUpLambda(omega_rank, prior_prob = rep(1/3, 3), prior_weights=c(0.5, 1, 2 ))
    
    # ----> update con    
    con <- GibbsUpConstant(beta_rank, NULL, mu_beta, omega_rank, NULL,con)
    
    #LRF TO ADDRESS: this is to be computed with the 'connections' dummy 0'd out
    mu = as.vector( X_program %*% beta_rank )
    if(!is.null(X_elite)){
      mu_noelite = as.vector( X_program_noelite %*% beta_rank  )
    }else{
      mu_noelite = mu
    }
    
    if(iter > iter_burn){
      j = iter - iter_burn
      # store value at this iteration
      draw$Z[j,] = apply(Z, 1, function(x){sum(x,na.rm=TRUE)})
      draw$beta_rank[j,] = beta_rank
      draw$mu[j,] = mu
      draw$mu_noelite[j,] = mu_noelite
      draw$omega_rank[j,] = omega_rank
      draw$con[j] = con
    }
    # print iteration number
    if(iter %% print_opt == 0){
      print(paste("Gibbs Iteration", iter))
      # print(table(weight.vec))
      # print(c(sigma2.alpha, sigma2.beta))
    }
  }
  
  
  return(draw)
}

