/* -*- compile-command: "R CMD INSTALL .. && R --vanilla < ../tests/testthat/test-eigen.R" -*- */

#include "eigen.h"
#include <Eigen/Dense>
#include <iostream>

int testEigen(double *X, int nrow, int ncol){
  Eigen::MatrixXd mat(nrow, ncol);
  int N=nrow*ncol;
  Eigen::Map< Eigen::MatrixXd > imported(X, nrow, ncol);
  Eigen::Map< Eigen::VectorXd > vec(X, N);
  vec(2) = 5;
  imported(0,1) += 1;
  imported(2,1) = (imported.row(0) - imported.row(1)).norm();
  imported(0,0) = 3.5;
  Eigen::VectorXi t = Eigen::VectorXi::LinSpaced(N,0,N-1);
  // for(int i=0; i<N; i++){
  //   t(i)=i;
  // }
  std::cout << t << std::endl;
  //std::sort(X, X+nrow*ncol, [&X](int lhs, int rhs){return X[lhs] < X[rhs];});
  //std::sort(t.data(), t.data()+t.size(), [&X](int lhs, int rhs){return X[lhs] < X[rhs];});
  std::sort(t.data(), t.data()+t.size(), [&vec](int lhs, int rhs){return vec(lhs) < vec(rhs);});
  std::cout << t << std::endl;
  return 0;
}
