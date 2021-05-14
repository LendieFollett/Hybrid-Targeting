#include <Rcpp.h>
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

// [[Rcpp::export]]
NumericVector GibbsUpLatentGivenRankInd(NumericMatrix pair_comp, NumericVector Z, NumericVector up_order, NumericVector mu, double weight ) {
  int N1 =Z.length();
  IntegerVector set1(N1);
  IntegerVector set0(N1);
  double upper;
  double lower;
  
  //for(i in up.order){
  for (int idx = 0; idx < N1; idx++){
    int i = up_order[idx];
      for (int idx2 = 0; idx2 < N1; idx2++){
        set1[idx2] = (pair_comp[i, idx2] == 1);
        set0[idx2] = (pair_comp[i, idx2] != 1);
        }
  NumericVector Z1 = Z[set1];
  NumericVector Z0 = Z[set1];
    
    if(sum(set1) > 0){
      upper = *std::min_element(Z1.begin(), Z1.end());
    }else{
      upper = std::numeric_limits<double>::infinity();
    }
    
    if(sum(set0) > 0){
      lower = *std::max_element(Z0.begin(), Z0.end());
    }else{
      lower = -std::numeric_limits<double>::infinity();;
    }
    
    Z[i] = R::rtruncnorm( 1, lower, upper, mean = mu[i], sd = 1/sqrt(weight) )
  }
  return(Z)

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
