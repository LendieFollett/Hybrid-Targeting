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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
