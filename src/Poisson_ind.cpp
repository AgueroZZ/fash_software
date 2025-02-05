#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  DATA_SPARSE_MATRIX(B); // Design matrix (for random effects)
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  
  int d = P.cols(); // Number of Spline coefficients
  int d_beta = X.cols(); // Number of fixed effect coefficients

  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)
  DATA_SCALAR(sigmaIWP); // theta = -2log(sigma)
  // Type thetaIWP = -2*log(sigmaIWP);
  DATA_VECTOR(offset);

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta), eta = offset + B * U + X * beta
  int Wdim = W.size();
  int betadim = Wdim - d;
  vector<Type> U(d);
  vector<Type> beta(betadim);
  for (int i=0;i<d;i++) U(i) = W(i);
  for (int i=0;i<betadim;i++) beta(i) = W(i+d);

  // Transformations
  vector<Type> eta = offset + X * beta + B * U;

  // Log likelihood
  Type ll = 0;
  ll = sum(dpois(y, exp(eta), TRUE));
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;

  // Cross product
  if(sigmaIWP != 0){
    vector<Type> PU = P*U;
    Type UPU = (U * PU).sum();
    lpW += -0.5 * exp(-2*log(sigmaIWP)) * UPU; // U part
  }
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part

  // Log determinant
  if(sigmaIWP != 0){
    Type logdet1 = d * (-2*log(sigmaIWP)) + logPdet;
    lpW += 0.5 * logdet1; // P part
  }
  Type logdet2 = d_beta * log(betaprec);
  lpW += 0.5 * logdet2; // for fixed effect


  // wrt the dimension
  lpW += -0.5 * d * log(2*M_PI);
  lpW += -0.5 * d_beta * log(2*M_PI);


  REPORT(lpW);
  
  // Final result!
  Type logpost = -1 * (ll + lpW);
  
  return logpost;
}