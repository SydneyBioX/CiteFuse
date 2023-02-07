#include <Rcpp.h>
#include <math.h>
#include "backend.h"
using namespace Rcpp;




// Function for vlr
// [[Rcpp::export]]
NumericMatrix lr2vlr(NumericMatrix lr){
  
  // Calculate variation matrix
  NumericMatrix x = clone(lr);
  NumericMatrix X = covRcpp(x, 0);
  int nfeats = lr.ncol();
  
  // Find diagonal
  NumericVector diag(nfeats);
  for(int j = 0; j < nfeats; j++){
    diag[j] = X(j, j);
  }
  
  // Calculate vlr
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      X(i, j) = -2 * X(i, j) + diag[i] + diag[j];
    }
  }
  
  return X;
}



// Function for rho
// [[Rcpp::export]]
NumericMatrix lr2rho(NumericMatrix lr){
  
  // Make vlr from log-ratio data
  NumericMatrix x = clone(lr);
  NumericMatrix mat = lr2vlr(x);
  int nsubjs = lr.nrow();
  int nfeats = lr.ncol();
  
  // Calculate variance of the i-th lr composition
  NumericVector vars(nfeats);
  for(int i = 0; i < nfeats; i++){
    
    vars[i] = sum(pow(lr(_, i) - mean(lr(_, i)), 2.0)) / (nsubjs - 1);
  }
  
  // Calculate rho = 1 - vlr[i, j] / (var[, i] + var[, j])
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      
      if(i == j){
        mat(i, j) = 1; // Force diagonal = 1
      }else{
        mat(i, j) = 1 - mat(i, j) / (vars[i] + vars[j]);
      }
    }
  }
  
  return mat;
}
