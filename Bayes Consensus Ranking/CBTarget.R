

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
  
  #If using random effects
  #ISSUE: ASSUMING EACH PERSON IS RANKED BY THE SAME NUMBER OF RANKERS
  if (any(apply(Z, 1, function(x){sum(!is.na(x))}) > 1)){ #evaluates to TRUE only when multiple ranks per household
    nrank = apply(Z, 1, function(x){sum(!is.na(x))}) #how many times was each of the N1 households ranked
    if(sd(nrank)>0){
      stop("Differing number of rankers per household")
    }
    X_RAND<- kronecker(diag(N1),rep(1, nrank[1])) #
  }
  
  if(is.null(initial.list$beta_rank)){beta_rank <- rep(0, P+1)}else{  beta_rank <-  initial.list$beta_rank } 
  if(is.null(initial.list$con)){con <- .5}else{ con <- initial.list$con  } 
  
  ## initial values for weights
  omega_rank = rep(1, R)
  #lambda <- prior_prob_rank
  #prior mean on beta_rank
  mu_beta = rep(0, length(beta_rank)-1)
  ## Gibbs iteration
  
  
  ## initial values for random effects
  alpha <- rep(0, N1)
  alpha_mat <- array(0, dim = dim(Z))
  sigma_alpha <- 2.5
  
  for(iter in 1:(iter_burn + iter_keep)){
    
    # update Z.mat given (alpha, beta) or equivalently mu
    Z = GibbsUpLatentGivenRankGroup(pair.comp.ten = pair.comp.ten, 
                                    Z = Z, 
                                    mu = X_CBT %*% beta_rank+ alpha, 
                                    omega_rank = omega_rank, 
                                    R = R )

    
    # ----> update beta_rank
    beta_rank = GibbsUpMuGivenLatentGroup(Y = Z - alpha_mat,
                                          X = X_CBT,
                                          omega = omega_rank,
                                          mu_beta = mu_beta,
                                          con = con,
                                          rank=TRUE)
    
    # ----> update quality weights
    omega_rank <- GibbsUpQualityWeights(y=Z , 
                                        mu=X_CBT %*% beta_rank, 
                                        beta_rank,  
                                        weight_prior_value = c(0.5, 1, 2 ), prior_prob = rep(1/3, 3),
                                        rank=TRUE)


    # ----> update con    
    con <- GibbsUpConstant(beta_rank, NULL, mu_beta, omega_rank, NULL,con)
    
    
    # ----> update random effect parameters IF multiple rankers per household
    # (this is kind of slow....)
    if (any(apply(Z, 1, function(x){sum(!is.na(x))}) > 1)){ #evaluates to TRUE only when multiple ranks per household
      
      alpha <- GibbsUpGammaGivenLatentGroup(Z,      X_CBT %*% beta_rank, X_RAND, omega_rank, sigma2_alpha = sigma2_alpha) 
      sigma2_alpha <- GibbsUpsigma_alpha(alpha, nu=3, tau2=25)  
      
      alpha_mat <- apply(Z, 1:2, function(x){ifelse(x != 0 & !is.na(x), 1, 0)})*alpha #reformatted alpha
    }
    
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

