#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

double fcalc_np(double &c_i, NumericVector &c_j) { //use pointers for speed
  //initialize c_j_sum
  double c_j_sum = 0;
  int n = c_j.size();

  // calculate the sum of c_j
  for(int i = 0; i < n; ++i) {
    c_j_sum += c_j[i];
  }
  double g_i = c_i * log(c_i / (c_j_sum));
  return g_i;
}

//' Function to calculate the network potential for vertices v
//'
//' @param neighbors list of neighbors for every node in the graph, type Rcpp::list
//' @param vertices node list for graph, type Rcpp::StringVector
//' @param v list of nodes for which we plan to calculate network potential
//' @param exp named vector of expression for each node in vertices
//'
// [[Rcpp::export]]

NumericVector fcalc_np_all(List &neighbors, StringVector &vertices, StringVector &v,
                           NumericVector &exp) {
  int n = vertices.size(); // define size of node list
  std::vector<double> np_vec(n);
  // initialize the variables we will use in the loop- does this save time?
  NumericVector neighbors_i;
  //StringVector neighbors_named_i;
  NumericVector c_j;
  double c_i;
  int neighbors_n;

  //logicalvector for whether we want to actually do calculations
  LogicalVector ind = in(vertices, v);

  for(int i = 0; i < n; ++i) {
    //don't do anything if its not one of our targets
    if(ind[i]){
      // do all the indexing to create a named list of neighbors for a given node v
      // then created named NumericVector with the expression for those neighbors
      // guarantee this could be done in fewer lines but oh well

      neighbors_i = neighbors[i];
      neighbors_n = neighbors_i.size();

      // need to correct for zero indexing in cpp vs 1 indexing in R
      for(int j = 0; j< neighbors_n; ++j) {
        neighbors_i[j] -= 1;
      }
      //lets just used numerical indexing
      //neighbors_named_i = v[neighbors_i];

      //define c_j
      c_j = exp[neighbors_i];

      // define c_i
      c_i = exp[i];

      //calculate network potential for each node and assign to new vector
      //hopefully it will be faster using numerical indexing with std::vector instead of
      //indexing based on the names of NumericVector
      np_vec.at(i) = fcalc_np(c_i = c_i, c_j = c_j);

    } else {np_vec.at(i) = 0;}
  }

  //convert back to a NumericVector
  //no way this will be faster but lets try our best
  NumericVector np_vec2 = wrap(np_vec);
  np_vec2.names() = vertices;

  return np_vec2;

}

