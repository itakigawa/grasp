#ifndef LOGREG_H_
#define LOGREG_H_

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp> 

#include <boost/math/distributions/chi_squared.hpp>

namespace ublas = boost::numeric::ublas;
namespace stat  = boost::math;

// aliases
typedef ublas::matrix<double> Matrix;
typedef ublas::vector<double> Vector;
typedef ublas::matrix_row<Matrix> MatrixRow;

// Solve A x = b via LU decomposition & backward substitution
template<typename MAT,typename V1,typename V2>
bool linear_solve(MAT& A, V1& b, V2& x){
  unsigned n = A.size1();
  if(n != A.size2()) return false;
  ublas::permutation_matrix<double> P(n);
  if(ublas::lu_factorize(A,P) != 0) return false;
  x.assign(b);                 // initialize
  ublas::lu_substitute(A,P,x); // overwrite x
  return true;
}

template<typename V1, typename V2>
double prob(V1& x, V2& theta){
  double val = inner_prod(x,theta);
  return exp(val)/(1+exp(val));
}

template<typename DMAT, typename V1, typename V2, typename V3>
double deviance(DMAT& X, V1& theta, V2& pos, V3& neg){
  double sum = 0.0;
  for(unsigned i=0; i<X.size1(); ++i){
    ublas::matrix_row<DMAT> x (X,i);
    double pr  = prob(x,theta);
    double num = pos(i)+neg(i);
    sum += pos(i)*(log(pos(i)/(num*pr)));
    sum += neg(i)*(log(neg(i)/(num*(1.0-pr))));
  }
  return 2.0*sum;
}

template<typename DMAT, typename V1, typename V2, typename V3>
double L(DMAT& X, V1& theta, V2& pos, V3& neg){
  double sum = 0.0;
  for(unsigned i=0; i<X.size1(); ++i){
    ublas::matrix_row<DMAT> x (X,i);
    double pr = prob(x,theta);
    sum += pos(i)*log(pr);
    sum += neg(i)*log(1.0-pr);
  }
  return sum;
}

template<typename DMAT, typename V1, typename V2>
double analysis_of_deviance(DMAT& X, V1& pos, V2& neg){

  unsigned m = X.size1();
  unsigned n = X.size2();

  Vector theta(n);
  theta.assign(ublas::zero_vector<double>(n));

  double val1 = 0.0, val2 = 0.0;
  unsigned iter = 0;
  Vector z(m), delta(n);
  Matrix D(m,m), M(n,n);

  do {
    D.assign(ublas::zero_matrix<double>(m));

    for(unsigned i=0; i<m; ++i){
      MatrixRow x (X,i);
      double pr  = prob(x,theta);
      double num = pos(i)+neg(i);
      D(i,i) = num * pr * (1.0-pr);
      z(i) = pos(i)-num*pr;
    }

    M.assign(prod(trans(X),Matrix(prod(D,X))));
    delta.assign(prod(trans(X),z));
    
    if(linear_solve(M,delta,delta)){
      val1 = L(X,theta,pos,neg);
      theta.plus_assign(delta);
      val2 = L(X,theta,pos,neg);
    }
    iter += 1;

  } while(fabs((val1-val2)/val1) > 1.0e-8 && iter < 25);
  
  return deviance(X,theta,pos,neg);
}

template<typename V1, typename V2>
double log10_p(V1& pos, V2& neg){
  for(unsigned i=0; i<pos.size(); ++i){
    if(pos(i)==0) return log10(1.0);
    if(neg(i)==0) return log10(1.0);
  }

  Matrix X1(4,3);
  X1(0,0) = 1;  X1(0,1) = 0;  X1(0,2) = 0;
  X1(1,0) = 1;  X1(1,1) = 0;  X1(1,2) = 1;
  X1(2,0) = 1;  X1(2,1) = 1;  X1(2,2) = 0;
  X1(3,0) = 1;  X1(3,1) = 1;  X1(3,2) = 1;

  double d1 = analysis_of_deviance(X1,pos,neg);

  Matrix X2(4,4);
  X2(0,0) = 1;  X2(0,1) = 0;  X2(0,2) = 0;  X2(0,3) = 0;
  X2(1,0) = 1;  X2(1,1) = 0;  X2(1,2) = 1;  X2(1,3) = 0;
  X2(2,0) = 1;  X2(2,1) = 1;  X2(2,2) = 0;  X2(2,3) = 0;
  X2(3,0) = 1;  X2(3,1) = 1;  X2(3,2) = 1;  X2(3,3) = 1;

  double d2 = analysis_of_deviance(X2,pos,neg);

  stat::chi_squared dist(1);
  double pval = stat::cdf(stat::complement(dist,fabs(d1-d2)));

  return log(pval)/log(10.0);
}

#endif /* LOGREG_H_ */
