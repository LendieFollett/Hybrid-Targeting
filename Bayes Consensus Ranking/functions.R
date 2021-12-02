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
                                      con,
                                      multiple_rankers){
  
  #LRF TO ADDRESS: Y_comm might be missing, Y_micro might be missing, ...assuming ranking will be there...
  #                same logic for corresponding x matrices

P <- ncol(X) - 1
  

  if(!rank){ #if it's micro
    u <- as.vector(Y)
    
    n <- length(u)
    c <- ncol(Y)
    Sigma_inv_y<-omega*diag(n) 
    if(!is.null(c)){
      X <-kronecker(rep(1, c), X) #(A*K)xP
    } 
  } else{ #if it's rank

    if (multiple_rankers){ 
      # condition for when a person is ranked by multiple sources
      u <- as.vector(Y) #stacked columns (ordered by ranker)
      rows <- which(!is.na(u))
      rows2 <- apply(Y, 2, function(x){which(!is.na(x))})
      u <- u[rows]
      n <- length(u)
      
        if(is.list(rows2)){ #if differing number of people per ranker*community combo
          X <- X[unlist(rows2),]
          Sigma_inv_y<-rep(omega, times =lapply(rows2, length))*diag(n)
          }else{
          X <- X[as.vector(rows2),]
          Sigma_inv_y<-rep(omega, times =rep(nrow(rows2), length(omega)))*diag(n)} 
    }else{
      u <- apply(Y, 1, function(x){sum(x, na.rm=TRUE)}) #basically take the only non-NA element
      n <- length(u)
      c <- 1
      Sigma_inv_y<-rep(omega, times = apply(Y, 2, function(x){sum(!is.na(x))}))*diag(n) 
      X <- X
    }
    #LRF:THIS ASSUMES WE'RE ORDERED BY COMMUNITY....
  }

if(!rank){  
  Sigma_inv_beta <- diag(c(1/2.5^2, (1/con)*rep(1, P))) #prior variance on non-intercept beta is sigma2_beta ; beta ~ N(mu_beta, sigma2_beta) 
  
  #A<-1x(A*K + M*N0 +R*N1)%*%square(A*K + M*N0 +R*N1)%*%(A*K + M*N0 +R*N1)xP --> 1xP

  pt1 <- u^T%*%Sigma_inv_y%*%X + t(c(0,mu_beta))%*%Sigma_inv_beta 
  
  pt2 <- t(X)%*%Sigma_inv_y%*%X + Sigma_inv_beta
  
  pt2_inv <- solve(pt2)
  
  alpha_beta <- mvrnorm(1, mu = t(pt1%*%pt2_inv), Sigma = pt2_inv)
  
}else{
  Sigma_inv_beta <- diag((1/con)*rep(1, P)) #prior covariance matrix (removed omega)
  
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

  Omega_rank <-  diag(rep(1, P))*con#diag(P + 1)*1/omega_rank  (removed omega)
  Omega_micro <- diag(rep(1, P))*con#diag(P + 1)*1/omega_micro  (removed omega)
  
  post_Sigma <- solve(solve(Omega_rank) + solve(Omega_micro) + diag(P)/1^2) #prior sd on mu_beta = 1
  
  post_mu <- (t(beta_rank[-1])%*%solve(Omega_rank) + t(beta_micro[-1])%*%solve(Omega_micro))%*%post_Sigma
  
  mu_beta <- mvrnorm(1, mu = t(post_mu), Sigma = post_Sigma)
  
  return(mu_beta)
}



### Gibbs update for information quality weights omega_micro, omega_rank---------
#y is Z (N1xR)
#prior_prob is weights in the case of updating omega_rank
GibbsUpQualityWeights <- function(y, mu, beta, weight_prior_value = c(0.5, 1, 2), prior_prob = rep(1/length(weight_prior_value), length(weight_prior_value)),rank=FALSE){
  
 #if it is omega_rank
    Col <- ncol(y) #need to update an omega for every oclumn

    n.prior.value <- length(weight_prior_value)
    weight_samp <- rep(NA, Col)
    
    log.post.prob = rep(0, n.prior.value) #re-initialize for next information source     
  for( col in 1:Col){ #over information source within y
   # log.post.prob = rep(0, n.prior.value) #re-initialize for next information source   
    idx <- which(!is.na(y[,col])) #in the case of omega_rank
    for(k in 1:n.prior.value){ #over potential values
        log.post.prob[k] <- log.post.prob[k] + sum(dnorm(y[idx,col], mu[idx], sqrt(1/weight_prior_value[k]), log = TRUE))#+ sum(dnorm(beta[-1], mean = mu_beta, sd =sqrt(con/weight_prior_value[k]), log = TRUE ))
    }
    
  }
    log.post.prob <- log.post.prob + log(prior_prob)
    log.post.prob = log.post.prob - max(log.post.prob)
    post.prob = exp(log.post.prob)
    weight_samp[1:Col] <- weight_prior_value[sample(c(1,2,3),size = 1, prob= post.prob)]
    #note: w = 1/sigma^2; sigma^2 = 1/w; sigma = 1/sqrt(w)
    return(weight_samp)

}


### Gibbs update for information quality weights omega_rank when ranks are of varying qualities---------
#y is either Y_comm (KxA), Y_micro (N0xM), or Z (N1xR)
#prior_prob will be current iteration's dirichlet weights in the case of updating omega_rank
#groups is an id vector of length ncol(y) = ncol(z), indexing unique rankers type e.g., groups = c(1,2,3,1,2,3,1,2,3,...)
#each unique grouping of rankers will share the same omega (score error variance)
GibbsUpQualityWeightsHeter <- function(y, groups, mu, beta, mu_beta,con, weight_prior_value = c(0.5, 1, 2), prior_prob = rep(1/length(weight_prior_value), length(weight_prior_value)),rank=FALSE){
  
 #if it is omega_rank
    R <- ncol(y) #need to update an omega for every column
    
    n.prior.value <- length(weight_prior_value)
    weight_samp <- rep(NA, R)
    
    for( g in unique(groups)){ #over unique ranker groups (e.g., women, leaders, etc...)
      log.post.prob = rep(0, n.prior.value) #re-initialize for next information source   
      cols <- which(groups == g) #within any group there will be multiple columns
        for(k in 1:n.prior.value){ #over potential values
          
          for (col in cols){ 
            idx <- which(!is.na(y[,col])) #in the case of omega_rank
            log.post.prob[k] <- log.post.prob[k] + sum(dnorm(y[idx,col], mu[idx], sqrt(1/weight_prior_value[k]), log = TRUE))
          }
          
        }
      log.post.prob <- log.post.prob + log(prior_prob[[g]]) #the g'th prior (specific to ranker type group)
      log.post.prob = log.post.prob - max(log.post.prob)
      post.prob = exp(log.post.prob)
      weight_samp[cols] <- sample(x=weight_prior_value,size = 1, prob= post.prob)
      
    }
    #note: w = 1/sigma^2; sigma^2 = 1/w; sigma = 1/sqrt(w)
    return(weight_samp)
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
  
  prior_old <- dnorm(con_old, 0, 1, log=TRUE)#prior constant ~ N^+(0,1)
  prior_prop<- dnorm(con_prop, 0, 1, log = TRUE)
  
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

# Update random effects for scores
#y is the Z score matrix
#mu is the corresponding mean x*delta
#omega is a vector error variance of the scores (potentially heterogeneous)
#sigma_alpha is the standard deviation on that source's random effect

GibbsUpGammaGivenLatentGroup <- function(y, xbeta, Xr, omega, sigma2_alpha = 2.5^2, Z_bin){
  
  N <- nrow(y) #number of random effects to estimate = number of rows
  #Update conditionally if ncol > 1
  Col <- ncol(y)
  
  #Complete 'data' vector
  resid <- y- xbeta%*%array(1, dim = c(1, ncol(y)))
  u <- c( t( resid) )  # stacked rows (ordered by person,then ranker within person) consistent with XRAND
  u <- u[complete.cases(u)]

  #calculate number of people per ranker
  #nperson = apply(y, 2, function(x){sum(!is.na(x))})
  Sigma_inv_diag <-c(t(Z_bin)*omega) 
  Sigma_inv_diag <- ifelse(Sigma_inv_diag == 0, NA, Sigma_inv_diag) 
  Sigma_inv_diag <- Sigma_inv_diag[complete.cases(Sigma_inv_diag)]

  pt1 <- (u*Sigma_inv_diag)^T%*%Xr
  #https://www.sciencedirect.com/topics/computer-science/diagonal-matrix
  pt2 <- (t(Xr)*Sigma_inv_diag)%*%Xr + diag(N)/(sigma2_alpha) #Inverse of posterior covariance 
  
  pt2_inv <- solve(pt2)
  
  alpha <- mvrnorm(1, mu = t(pt1%*%pt2_inv), Sigma = pt2_inv)
  
  return(alpha)
}


### Gibbs update for variances on random effects
### Gibbs update for sigma2, given prior sigma2 ~ Scale-Inv-chi2(nu, tau2) and data iid ~ N(0, sigma2)
GibbsUpsigma_alpha <- function(x, nu, tau2){
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

