// [[Rcpp::depends(RcppEigen)]]
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>


//' Rcpp function for matrix
//' @param A a matrix
//' @param B a matrix
//' @param C a matrix
//'
//' @return A matrix
//'
//' @export
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, Eigen::Map<Eigen::MatrixXd> C){
  Eigen::MatrixXd D = A * B * C;

  return Rcpp::wrap(D);
}
