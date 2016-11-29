#include <RcppArmadillo.h>
using namespace Rcpp;


//[[Rcpp::export]]
Rcpp::List nuclearFPCA(List yList, List BList, arma::mat Omega, double lambda1,
                      double lambda3){
  //import data
  //int n=Rcpp::as<int>(n1), q=Rcpp::as<int>(q1),maxIter=Rcpp::as<int>(maxIter1);
  //double lambda1=as<double>(lambda11),lambda2=as<double>(lambda21),lambda3=as<double>(lambda31);
  //Rcpp::NumericMatrix OmegaRcpp(Omega1);
  //arma::mat Omega(OmegaRcpp.begin(),OmegaRcpp.nrow(),OmegaRcpp.ncol());
  //import B and y
  //Rcpp::List yList(y1), BList(B1);
  
  int n = yList.size();
  int q = Omega.n_cols;
  int maxIter = 1000;
  std::vector<arma::mat> y(n);
  std::vector<arma::mat> B(n);
  int i,j;
  for(i=0;i<n;i++){
    SEXP yl=yList[i];
    SEXP Bl=BList[i];
    Rcpp::NumericVector yRcpp(yl);
    Rcpp::NumericMatrix BRcpp(Bl);
    y[i]=arma::mat(yRcpp.begin(), yRcpp.size(), 1);
    B[i]=arma::mat(BRcpp.begin(), BRcpp.nrow(), q);
  }
  //initialization
  arma::mat S(q,n,arma::fill::zeros);
  arma::mat P(S),U(S),R(S),D(S),PreS(S);
  arma::mat Ru,Rv;// store the SVD of R
  arma::vec Rd;
  arma::mat theta(q,1,arma::fill::zeros);
  arma::mat left(q,q,arma::fill::zeros);
  arma::mat right(q,1,arma::fill::zeros);
  int iter=0;bool flag=true;
  double rho=1e-4,error,error2;
  while(flag && iter<maxIter){
    iter++;
    //step 1
    //  left=lambda2*Omega;
    //  right.zeros();
    //  for(i=0;i<n;i++){
    //    left=left+B[i].t()*B[i];
    //    right=right+B[i].t()*y[i]-B[i].t()*B[i]*S.col(i);
    //  }
    //  theta=arma::solve(left,right);  
    //step 2
    //  for(i=0;i<n;i++){
    //       left=lambda3*Omega+B[i].t()*B[i]+rho*arma::eye<arma::mat>(q,q);
    //       right=B[i].t()*y[i]-B[i].t()*B[i]*theta+rho*P.col(i)+0.5*U.col(i);
    //       S.col(i)=arma::solve(left,right);
    //  }
    // step 2-2
    for(i=0;i<n;i++){
      left=lambda3*Omega+B[i].t()*B[i]+rho*arma::eye<arma::mat>(q,q);
      right=B[i].t()*y[i]+rho*P.col(i)+0.5*U.col(i);
      S.col(i)=arma::solve(left,right);
    }
    //step 3
    R=S-U/(2*rho);
    arma::svd(Ru,Rd,Rv,R);
    D.zeros();
    for(j=0;j< ((int) Rd.size());j++){
      Rd[j]=Rd[j]-lambda1/rho;
      if(Rd[j]<0.00001) Rd[j]=0;
      D(j,j)=Rd[j];
    }
    P=Ru*D*Rv.t();
    //step 4
    U=U+2*rho*(P-S);
    rho=rho*1.01;
    if(rho>1e7) rho=1e7;
    error=arma::norm(arma::vectorise(S-PreS),2);
    error2=arma::norm(arma::vectorise(S),2);
    if(error/error2<(1e-4)) {
      flag=false;
    }
    PreS=S;
  }
  return Rcpp::List::create(Rcpp::Named("theta")=theta,Rcpp::Named("S")=S,
                            Rcpp::Named("D")=iter, Rcpp::Named("E")=error);
}