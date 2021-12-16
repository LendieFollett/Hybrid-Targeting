
#If multiple rankers and heterogeneous qualities:
# - groups will index ranker type (length(groups) = ncol(Tau))
# - prior_prob_rank will be a list with # of elements equal to the # of distinct ranker types, in order

#CBT RANKING ONLY
CBTarget<- function(Tau, X_CBT=NULL, X_program=NULL,
                    X_elite = NULL,
                    weight_prior_value = c(0.5, 1, 2), 
                    groups = rep(1, ncol(Tau)), #Defaults to homogeneous weights
                    N1 = dim(X_CBT)[1], #how many people in test set
                    R = ncol(Tau), #how many rankers. often will be equal to K
                    iter_keep = 5000,
                    iter_burn = 5000,
                    print_opt = 100,
                    initial.list = NULL,
                    delta_prior_mean = rep(0, length(beta_rank)-1), #FOR DYNAMIC UPDATING
                    delta_prior_var = 2.5^2, #FOR DYNAMIC UPDATING
                    prior_prob_rank = list(rep(1/length(weight_prior_value), length(weight_prior_value)))#FOR DYNAMIC UPDATING
                    ){
  # DIMENSIONS ######################
  
  #pair.com.ten A list of length R with elements: \eqn{N_CBT[k]} by \eqn{N_CBT[k]}  pairwise comparison array where N_CBT[k] is the number of people in community k (ranked by ranker r)
  #where the (\eqn{i},\eqn{j}) element equals 1 if \eqn{i} is ranked higher than \eqn{j} by ranker \eqn{r}, 
  #0 if \eqn{i} is ranker lower than \eqn{j}, 
  #and NA if the relation between \eqn{i} and \eqn{j} is missing. 
  #Note that the diagonal elements (\eqn{i},\eqn{i},\eqn{r})'s for all rankers should be set to NA as well.
  N_CBT <- dim(X_CBT)[1] #how many people in test set 
  R <- ncol(Tau) #how many rankers total (not per community, but overall)

  
  # STRUCTURES ######################
  
  #Matrix of latent scores
  Z = matrix(NA, nrow = N_CBT, ncol = R)
  
  #create pair.comp.ten matrix
  pair.comp.ten = list()
  for(r in 1:R){
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
  
  
  
  # INITIAL VALUES ###################### 
  
  ## set initial values for parameters, where given
  if(is.null(initial.list$Z)){
    for(j in 1:R){
      rcases <- which(!is.na(Tau[,j]))
      nranked <- length(rcases)
      Z[rcases,j][sort( rowSums( pair.comp.ten[[j]], na.rm = TRUE ), decreasing = FALSE, index.return = TRUE )$ix]= (c(nranked : 1) - (1+nranked)/2)/sd(c(nranked : 1))
    }
  }else{
    Z <- initial.list$Z
  }
  
  #If using random effects
  #for conditional random effect logic
  multiple_rankers <- any(apply(Z, 1, function(x){sum(!is.na(x))}) > 1)
  
  #ISSUE: ASSUMING EACH PERSON IS RANKED BY THE SAME NUMBER OF RANKERS
  if (multiple_rankers){ #evaluates to TRUE only when multiple ranks per household
    nrank = apply(Z, 1, function(x){sum(!is.na(x))}) #how many times was each of the N_CBT households ranked
    if(sd(nrank)>0){
      stop("Differing number of rankers per household")
    }
    #Random effect design matrix: ordered by household, then ranker within household
    X_RAND<- kronecker(diag(N_CBT),
                       rep(1, nrank[1])) #
    #Binary matrix indicating positions of ranker*household observations
    Z_bin <- apply(Z, 1:2, function(x){ifelse( !is.na(x), 1, 0)}) #x != 0 & removed
    n_non_na <- sum(Z_bin)
  }

  
  if(is.null(initial.list$beta_rank)){beta_rank <- rep(0, P+1)}else{  beta_rank <-  initial.list$beta_rank } 

  
  ## initial values for weights
  omega_rank = rep(1, R)
  con <- 1
  ## Gibbs iteration
  
  ## initial values for random effects
  alpha <- rep(0, N_CBT)
  alpha_mat <- array(0, dim = dim(Z))
  sigma2_alpha <-1
  
  
  ## store MCMC draws
  draw = list(
    Z = array(0, dim = c( N_CBT,R)),
    beta_rank = array(NA, dim = c(iter_keep,P+1)),
    mu = array(NA, dim = c(iter_keep,nrow(X_program))),
    mu_noelite = array(NA, dim = c(iter_keep,nrow(X_program))), 
    omega_rank = array(NA, dim = c(iter_keep, R) ),
    alpha = array(NA, dim = c(iter_keep, N_CBT)),
    sigma2_alpha = array(NA, dim = c(iter_keep, 1))
  )
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
                                          mu_beta = delta_prior_mean, #possibly using dynamic updating
                                          con = delta_prior_var, #possibly using dynamic updating
                                          rank=TRUE,
                                          multiple_rankers = multiple_rankers)
    
    # ----> update quality weights
    omega_rank <- GibbsUpQualityWeightsHeter(y=Z , 
                                             groups = groups,
                                        mu=X_CBT %*% beta_rank + alpha, 
                                        beta = beta_rank,  
                                        weight_prior_value = weight_prior_value, 
                                        prior_prob = prior_prob_rank,#possibly using dynamic updating
                                        rank=TRUE)


    # ----> update con    
    #LRF: the prior variance on delta (rank coefs) can be fixed, potentially uisng dynamic updating
   # con <- GibbsUpsigma_alpha(beta_rank[-1], nu=1, tau2=1)#GibbsUpConstant(beta_rank, beta_micro, mu_beta, omega_rank, omega_micro,con)
    
  
    
    # ----> update random effect parameters IF multiple rankers per household
    # (this is kind of slow....)
    if (multiple_rankers){ #evaluates to TRUE only when multiple ranks per household
      
      alpha  <- GibbsUpGammaGivenLatentGroupRCPP(y=Z, 
                                                 xbeta = X_CBT %*% beta_rank, 
                                                 Xr = X_RAND, 
                                                 omega_rank, 
                                                 sigma2_alpha = sigma2_alpha,
                                                 n_non_na = n_non_na) %>%c()
      
      
      sigma2_alpha <- 1#GibbsUpsigma_alpha(alpha, nu=1, tau2=1)  
      
      alpha_mat <- Z_bin*alpha #reformatted alpha
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
      draw$Z = draw$Z+Z/iter_keep#apply(Z, 1, function(x){sum(x,na.rm=TRUE)})
      draw$beta_rank[j,] = beta_rank
      draw$mu[j,] = mu
      draw$mu_noelite[j,] = mu_noelite
      draw$omega_rank[j,] = omega_rank
      #draw$con[j] = con
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

