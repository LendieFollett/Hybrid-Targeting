
#' Bayesian Consensus Targeting With Multiple Rankers
#' Includes random effect in score model, the possibility for
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
                        prior_prob_rank = list(rep(1/length(weight_prior_value), length(weight_prior_value))), #override if heterogeneous
                        groups = rep(1, ncol(Tau)), #Defaults to homogeneous weights
                        N1 = dim(X_CBT)[1], #how many people in test set
                        R = ncol(Tau), #how many rankers. often will be equal to K
                        iter_keep = 5000,
                        iter_burn = 5000,
                        print_opt = 100,
                        initial.list = NULL){
  #pair.com.ten A list of length R with elements: \eqn{N1[k]} by \eqn{N1[k]}  pairwise comparison array
  #where the (\eqn{i},\eqn{j}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{r}, 
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
  
  #If using random effects
  #ISSUE: ASSUMING EACH PERSON IS RANKED BY THE SAME NUMBER OF RANKERS
  if (any(apply(Z, 1, function(x){sum(!is.na(x))}) > 1)){ #evaluates to TRUE only when multiple ranks per household
    nrank = apply(Z, 1, function(x){sum(!is.na(x))}) #how many times was each of the N1 households ranked
    if(sd(nrank)>0){
      stop("Differing number of rankers per household")
    }
    X_RAND<- kronecker(diag(N1),rep(1, nrank[1])) #
    Z_bin <- apply(Z, 1:2, function(x){ifelse(x != 0 & !is.na(x), 1, 0)})
  }
  
  if(is.null(initial.list$beta_rank)){beta_rank <- rep(0, P+1)}else{  beta_rank <-  initial.list$beta_rank } 
  if(is.null(initial.list$beta_micro)|is.null(Y_micro)){beta_micro <- rep(0, P+1)}else{  beta_micro <-  initial.list$beta_micro } 
  if(is.null(initial.list$con)){con <- .5}else{ con <- initial.list$con  } 
  mu_beta <- cbind(beta_rank[-1], beta_micro[-1]) %>%apply(1, mean)
  
  ## initial values for weights
  omega_micro = 1
  omega_rank = rep(1, R)
  ## Gibbs iteration
  
  ## initial values for random effects
  alpha <- rep(0, N1)
  alpha_mat <- array(0, dim = dim(Z))
  sigma2_alpha <- 2.5^2
  
  ## store MCMC draws
  draw = list(
    Z = array(NA, dim = c( iter_keep,N1)),
    mu_beta = array(NA, dim = c(iter_keep,P)),
    beta_rank = array(NA, dim = c(iter_keep,P+1)),
    beta_micro = array(NA, dim = c(iter_keep,P+1)),
    mu = array(NA, dim = c(iter_keep,nrow(X_program))),
    mu_noelite = array(NA, dim = c(iter_keep,nrow(X_program))), 
    omega_micro = array(NA, dim = c(iter_keep, 1) ),
    omega_rank = array(NA, dim = c(iter_keep, R) ),
    con = array(NA, dim = c(iter_keep, 1) ),
    alpha = array(NA, dim = c(iter_keep, N1)),
    sigma2_alpha = array(NA, dim = c(iter_keep, 1))
  )
  
  for(iter in 1:(iter_burn + iter_keep)){
    
    # ----> update Z 
    Z <- GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, 
                                    Z = Z, 
                                    mu = X_CBT %*% beta_rank + alpha, 
                                    omega_rank = omega_rank, 
                                    R = R )

    # ----> update beta_rank
    beta_rank <- GibbsUpMuGivenLatentGroup(Y = Z -alpha_mat,
                                          X = X_CBT,
                                          omega = omega_rank,
                                          mu_beta = mu_beta,
                                          con = con,
                                          rank=TRUE)
    
    # ----> update quality weights, potentially heterogeneous
    omega_rank <- GibbsUpQualityWeightsHeter(y=Z , 
                                             groups = groups,
                                        mu=X_CBT %*% beta_rank, 
                                        beta_rank, 
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob =prior_prob_rank,
                                        rank=TRUE)
    
    
    # ----> update beta_micro 
    beta_micro <- GibbsUpMuGivenLatentGroup(Y = Y_micro,
                                           X = X_PMT,
                                           omega = omega_micro,
                                           mu_beta = mu_beta, 
                                           con = con)
    
    # ----> update quality weights    
    omega_micro <- 1/GibbsUpsigma_alpha(Y_micro-X_PMT %*% beta_micro, nu=3, tau2=25)  
    
    
    # ----> update mu_beta
    mu_beta <- GibbsUpGlobalMuGivenMu(beta_rank,  beta_micro,
                                      omega_rank, omega_micro ,con)
    
    # ----> update Omega (shared variance of delta, gamma around mu_beta)
    con <- GibbsUpsigma_alpha(c(beta_rank[-1], beta_micro[-1]) - c(mu_beta, mu_beta), nu=3, tau2=25)#GibbsUpConstant(beta_rank, beta_micro, mu_beta, omega_rank, omega_micro,con)
    
    
    # ----> update random effect parameters IF multiple rankers per household
    # (this is kind of slow....)
    if (any(apply(Z, 1, function(x){sum(!is.na(x))}) > 1)){ #evaluates to TRUE only when multiple ranks per household
      
    alpha <- GibbsUpGammaGivenLatentGroup(Z,      X_CBT %*% beta_rank, X_RAND, omega_rank, sigma2_alpha = sigma2_alpha) 
    sigma2_alpha <- GibbsUpsigma_alpha(alpha, nu=3, tau2=25)  
    
    alpha_mat <- Z_bin*alpha #reformatted alpha
    }
    
    
    #LRF TO ADDRESS: this is to be computed with the 'connections' dummy 0'd out
    mu = as.vector( X_program[,-1] %*% beta_rank[-1] )
    if(!is.null(X_elite)){
      mu_noelite = as.vector( X_program_noelite[,-1] %*% beta_rank[-1] ) 
    }else{
      mu_noelite = mu
    }
    
    if(iter > iter_burn){
      j = iter - iter_burn
      # store value at this iteration
      draw$Z[j,] = apply(Z, 1, function(x){sum(x,na.rm=TRUE)})
      draw$mu_beta[j,] = mu_beta
     # draw$lambda[j,] = lambda
      draw$beta_rank[j,] = beta_rank
      draw$beta_micro[j,] = beta_micro
      draw$mu[j,] = mu
      draw$mu_noelite[j,] = mu_noelite
      draw$omega_micro[j] = omega_micro
      draw$omega_rank[j,] = omega_rank
      draw$con[j] = con
      draw$alpha[j,] = alpha
      draw$sigma2_alpha[j,] = sigma2_alpha
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
