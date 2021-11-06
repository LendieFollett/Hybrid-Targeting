#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppTN.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppTN)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



NumericVector GibbsUpLatentGivenRankInd2(NumericMatrix pair_comp, NumericVector Z, IntegerVector up_order, NumericVector mu, double weight ) {
  int N1 =Z.length();
  IntegerVector set1(N1);
  IntegerVector set0(N1);
  NumericVector Z1(N1);
  NumericVector Z0(N1);
  double upper;
  double lower;
  int i;

  //for(i in up.order){
  for (int idx = 0; idx < N1; idx++){
    i = up_order[idx]; //this assumes up_order is of length N1
      for (int idx2 = 0; idx2 < N1; idx2++){
        set1[idx2] = (pair_comp(i, idx2) == 1); // i is ranked higher than idx2 (i has better wellbeing)
        set0[idx2] = (pair_comp(i, idx2) != 1); // otherwise
        }
   Z1 = Z[set1]; //latent score samples of those are ranked higher (have lower true Z values)
   Z0 = Z[set0]; //latent score samples of those are ranked lower (have higher true Z values)

    if(sum(set1) > 0){ // if there is anyone ranked higher (with a lower true Z)
      upper = *std::min_element(Z1.begin(), Z1.end());
    }else{
      upper = std::numeric_limits<double>::infinity();
    }

    if(sum(set0) > 0){
      lower = *std::max_element(Z0.begin(), Z0.end());
    }else{
      lower = -std::numeric_limits<double>::infinity();
    }
    
do{
  Z[i] = R::rnorm(mu[i], 1/sqrt(weight) );
} while ( Z[i] < lower || Z[i] > upper);

  }
  return(Z);

}

NumericVector reprcpp(NumericVector x, NumericVector y) {
  int n = y.size();
  std::vector<double> myvector(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    std::fill(myvector.begin()+ind, myvector.begin()+ind+y[i], x[i]);
    ind += y[i];
  }
  return Rcpp::wrap(myvector);
}



arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  ////std::cout << sigma << std::endl;
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::vec GibbsUpGammaGivenLatentGroupRCPP(arma::mat y, arma::vec xbeta, arma::mat Xr, arma::vec omega, double sigma2_alpha ){
  int N = y.n_rows;
  int Col = y.n_cols;
  
  arma::mat temp(1, Col);
  temp.ones();
  
  
//Complete 'data' vector
  arma::mat resid = y- (xbeta*temp);
// stack columns

int n_non_na = 0;
for (int idx = 0; idx < Col; idx++){
  for (int idx2 = 0; idx2 < N; idx2++){
  n_non_na +=   (1-Rcpp::NumericVector::is_na(y(idx2, idx)));
  }
}

arma::vec u(n_non_na);

arma::mat Sigma_inv_diag(n_non_na, n_non_na);
arma::vec nperson(Col);
nperson.zeros();



int j = -1;
for (int idx = 0; idx < Col; idx++){
  for (int idx2 = 0; idx2 < N; idx2++){
    if( Rcpp::NumericVector::is_na(y(idx2, idx)) == 1 ){
      //j = j;
      }else{
    j +=1;
    u(j) = resid(idx2, idx);
    nperson(idx) = nperson(idx) + 1;
    Sigma_inv_diag(j,j) = omega(idx);
      }
    }
}

    arma::mat UnityMatrix = arma::eye(Xr.n_cols,Xr.n_cols);
std::cout << u.submat(0,0,20,0) << std::endl;
    arma::mat pt1 =arma::trans(u)*Sigma_inv_diag*Xr;
//    std::cout << "second  "  << std::endl;
//https://www.sciencedirect.com/topics/computer-science/diagonal-matrix
  arma::mat pt2_prelim(Xr.n_cols, Xr.n_rows);
  arma::mat Xr_trans = arma::trans(Xr);
  
  for (int c = 0; c < Xr_trans.n_cols; c++){
    pt2_prelim.col(c) = Xr_trans.col(c)*Sigma_inv_diag(c,c);
  }
    arma::mat pt2 = pt2_prelim*Xr + UnityMatrix/(sigma2_alpha);// #Inverse of posterior covariance 
   // std::cout << "third  "  << std::endl;
    arma::mat pt2_inv = inv(pt2);
  //  std::cout << "inv  "  << std::endl;
    arma::vec alpha = arma::trans(mvrnormArma(1, arma::trans(pt1*pt2_inv), pt2_inv));
    
  return alpha;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
