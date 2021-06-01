#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//helper function to initialiZe the affinity vector for cpp_RWR

arma::mat init_p(arma::mat p, CharacterVector seed_proteins,
                 CharacterVector rownames) {

  LogicalVector ind = in(rownames, seed_proteins);
  double num_seeds = seed_proteins.size(); //define numer of seeds for p vec initialization
  int n = ind.size(); //

  for(int i = 0; i < n; ++i) {
    if(ind[i]) {
      p[i] = 1 / num_seeds; //normalize start probability by the number of seeds we have
    }
  }
  return p;
}


//helper function to check if the diff between p and pold is < eps

double check_p(arma::mat p, arma::mat pold) {

  arma::mat p_diff = p - pold;

  //prep to sum the absolute value of p_diff
  double p_sum = 0;
  int n = p_diff.n_rows;

  //iterate through the matrix
  for(int i = 0; i<n; ++i) {
    p_sum += std::abs(p_diff[i]);
  }

  return p_sum;

}

// main function to perform random walk with restarts on a sparse matrix
// [[Rcpp::export]]

Rcpp::List cpp_RWR(arma::sp_mat w, CharacterVector rownames,
                   CharacterVector seed_proteins, int tmax = 1000,
                   double gamma = 0.6, double eps = 1e-10) {

  int n = sqrt(w.size()); // return the length of our square adjacency matrix
  arma::mat p(n, 1, arma::fill::zeros); //initialize p vector

  // Normalize starting probability by the number of seeds
  p = init_p(p = p, seed_proteins = seed_proteins, rownames = rownames);

  //initialize p0
  arma::mat p0 = p;

  arma::mat pold; // initialize pold
  w = w.t(); // transpose w for calculations

  int t; //initialize iterator
  //loop through tmax
  for(t = 0; t < tmax; ++t) {
    pold = p; //reset pold to p of the previous iteration

    //do calculation for a given iteration
    p = ((1 - gamma) * (w * pold)) + gamma * p0;

    //check if exit condition has been met.
    if(check_p(p = p, pold = pold) < eps) {break;}
  }

  return Rcpp::List::create(Rcpp::Named("p") = p,
                            Rcpp::Named("seed_proteins") = seed_proteins,
                            Rcpp::Named("n.iter") = t);
}


