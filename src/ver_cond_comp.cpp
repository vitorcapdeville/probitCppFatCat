#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

// Rcpp::NumericVector zrnorm(int n) {
//   Rcpp::NumericVector x(n);
//   for (int i=0; i<n; i++) {
//     x[i] = zigg.norm();
//   }
//   return x;
// }

int indica(bool cond){
  if(cond){
    return 1;
  } else {
    return 0;
  }
}

double innerProduct(NumericVector x,
                    NumericVector y) {
  double val = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
  return val;
}


// [[Rcpp::export]]
NumericVector probProbit(NumericVector alfa, NumericVector beta, NumericVector f, double sigma2){
  int k = alfa.size() + 1;
  NumericVector prob(k);
  prob[0] = R::pnorm((alfa[0] - innerProduct(beta,f))/std::sqrt(sigma2),0.0,1.0,1,0);
  double sum = prob[0];
  for(int s = 1; s <= k - 2; s++){
    prob[s] = R::pnorm((alfa[s] - innerProduct(beta,f))/std::sqrt(sigma2),0.0,1.0,1,0)-R::pnorm((alfa[s-1] - innerProduct(beta,f))/std::sqrt(sigma2),0.0,1.0,1,0);
    sum += prob[s];
  }
  prob[k-1] = 1 - sum;
  return prob;
}


// [[Rcpp::export]]
NumericVector probLogit(NumericVector alfa, NumericVector beta, NumericVector f, double sigma2){
  int k = alfa.size() + 1;
  NumericVector prob(k);
  prob[0] = R::plogis(alfa[0],innerProduct(beta,f),std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0);
  double sum = prob[0];
  for(int s = 1; s <= k - 2; s++){
    prob[s] = R::plogis(alfa[s],innerProduct(beta,f),std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0)-R::plogis(alfa[s-1],innerProduct(beta,f),std::sqrt(sigma2 * 3/(pow(M_PI,2))),1,0);
    sum += prob[s];
  }
  prob[k-1] = 1 - sum;
  return prob;
}

// Só funciona para p = 1
// [[Rcpp::export]]
NumericVector probNominal(NumericMatrix beta, double f){
  int k = beta.ncol() + 1;
  NumericVector prob(k);
  NumericVector aux(k-1);
  aux[0] = exp(beta(0,0)*f);
  aux[1] = exp(beta(0,1)*f);
  aux[2] = exp(beta(0,2)*f);
  double aux2 = sum(aux);
  for(int s = 0; s <= k - 2; s++){
    prob[s] = exp(beta(0,s)*f)/(1 + aux2);
  }
  prob[k-1] = 1/(1 + aux2);
  return prob;
}

// [[Rcpp::export]]
double logcondcompbetajProbit(NumericMatrix f,NumericVector beta,NumericVector alfa,double sigma2, NumericMatrix y,int j){
  double aux = 0;
  int n = f.ncol();
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i) - 1;
    aux += std::log(probProbit(alfa,beta,f(_,i),sigma2)[ind]);
  }
  double aux2 = 0;
  int t=0;
  while(t<j-1 && t < beta.size()-1){
    aux2 -= (1/2000000)*std::pow(beta[t],2);
    t += 1;
  }
  if(t==j-1){
    aux2 -= (1/2000000)*std::pow(beta[t],2);
  }
  return aux+aux2;
}



// [[Rcpp::export]]
double logcondcompfiProbit(NumericVector f, NumericMatrix beta, NumericMatrix alfa, NumericVector sigma2,NumericMatrix y,int i){
  double aux = 0;
  int j = sigma2.size();
  int p = f.size();
  for(int l = 0; l <= j-1; l++){
    int ind = y(l,i-1) - 1;
    aux += log(probProbit(alfa(l,_),beta(l,_),f,sigma2[l])[ind]);
  }
  double aux2 = 0;
  for(int t=0; t<=p-1;t++){
    aux2 -= std::pow(f[t],2)/2;
  }
  return(aux + aux2);
}


// [[Rcpp::export]]
double logcondcompsigma2jProbit (NumericMatrix f,NumericVector beta,NumericVector alfa,double sigma2,NumericMatrix y, int j){
  double aux = 0;
  int n = f.ncol();
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i)-1;
    aux += log(probProbit(alfa,beta,f(_,i),sigma2)[ind]);
  }
  double aux2 = aux + 0.999*std::log(sigma2) - 0.001/sigma2;
  return(aux2); 
}

// [[Rcpp::export]]
double logcondcompbetajLogit(NumericMatrix f,NumericVector beta,NumericVector alfa,double sigma2, NumericMatrix y,int j){
  double aux = 0;
  int n = f.ncol();
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i) - 1;
    aux += std::log(probLogit(alfa,beta,f(_,i),sigma2)[ind]);
  }
  double aux2 = 0;
  int t=0;
  while(t<j-1 && t < beta.size()-1){
    aux2 -= (1/2000000)*std::pow(beta[t],2);
    t += 1;
  }
  if(t==j-1){
    aux2 -= (1/2000000)*std::pow(beta[t],2);
  }
  return aux+aux2;
}



// [[Rcpp::export]]
double logcondcompfiLogit(NumericVector f, NumericMatrix beta, NumericMatrix alfa, NumericVector sigma2,NumericMatrix y,int i){
  double aux = 0;
  int j = sigma2.size();
  int p = f.size();
  for(int l = 0; l <= j-1; l++){
    int ind = y(l,i-1) - 1;
    aux += log(probLogit(alfa(l,_),beta(l,_),f,sigma2[l])[ind]);
  }
  double aux2 = 0;
  for(int t=0; t<=p-1;t++){
    aux2 -= std::pow(f[t],2)/2;
  }
  return(aux + aux2);
}


// [[Rcpp::export]]
double logcondcompsigma2jLogit (NumericMatrix f,NumericVector beta,NumericVector alfa,double sigma2,NumericMatrix y, int j){
  double aux = 0;
  int n = f.ncol();
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i)-1;
    aux += log(probLogit(alfa,beta,f(_,i),sigma2)[ind]);
  }
  double aux2 = aux + 0.999*std::log(sigma2) - 0.001/sigma2;
  return(aux2); 
}


// [[Rcpp::export]]
double logcondcompbetajkNominal(NumericVector f,NumericMatrix beta, NumericMatrix y,int j, int k){
  double aux = 0;
  int n = f.size();
  for(int i = 0; i <= n-1; i++){
    int ind = y(j-1,i) - 1;
    aux += std::log(probNominal(beta,f[i])[ind]);
  }
  double aux2 = 0;
  int t=0;
  NumericVector aux3 = beta(_,k-1);
  while(t<j-1 && t < aux3.size()-1){
    aux2 -= (1/2000000)*std::pow(beta(t,k-1),2);
    t += 1;
  }
  if(t==j-1){
    aux2 -= (1/2000000)*std::pow(beta(t,k-1),2);
  }
  return aux+aux2;
}

// [[Rcpp::export]]
double logcondcompfiNominal(double f, NumericMatrix beta,NumericMatrix y,int i){
  double aux = 0;
  int j = beta.nrow();
  int k = beta.ncol();
  int p = 1;
  for(int l = 0; l <= j-1; l++){
    int ind = y(l,i-1) - 1;
    NumericMatrix aux3(1,k);
    aux3(0,_) = beta(l,_);
    aux += log(probNominal(aux3,f)[ind]);
  }
  double aux2 = 0;
  for(int t=0; t<=p-1;t++){
    aux2 -= std::pow(f,2)/2;
  }
  return(aux + aux2);
}


// // [[Rcpp::export]]
// List algoritmoProbit(int nit, int j, int p, int n, NumericMatrix alfa, NumericMatrix y,NumericMatrix f, NumericVector sigma2){
//   List beta(nit);
//   NumericMatrix inibeta(j,p);
//   inibeta(0,0) = 0.5;
//   beta[0] = inibeta;
//   // List sigma2(nit);
//   // sigma2[0] = rep(2,j);
//   // List f(nit);
//   // NumericMatrix inif(p,n);
//   // f[0] = inif;
//   double sdpropbeta = 0.05;
//   double sdpropbeta2 = 0.075;
//   // double sdpropf = 0.4;
//   // double sdpropsigma2 = 0.075;
//   NumericVector betaprop(p);
//   double A1 = 0.0;
//   double A2 = 0.0;
//   double A = 0.0;
//   // double A3 = 0.0;
//   // double A4 = 0.0;
//   // double B = 0.0;
//   // double A5 = 0.0;
//   // double A6 = 0.0;
//   // double C = 0.0;
//   double pa = 0.0;
//   double u = 0.0;
//   // double sigma2prop = 0.0;
// 
//   // NumericVector auxls(j);
//   NumericMatrix auxlb(j,p);
//   // NumericMatrix auxlf(p,n);
//   NumericMatrix aux4(j,p);
//   // NumericMatrix aux5(p,n);
//   NumericVector auxvec(2);
//   // NumericVector fprop(p);
// 
// 
//   for(int it=1; it <= nit-1;it++){
//     // auxls = sigma2[it-1];
//     aux4 = clone(as<NumericMatrix>(beta[it-1]));
//     // aux5 = clone(as<NumericMatrix>(f[it-1]));
//     for(int l=0; l <= j-1; l++){
//       for(int t=0; t<= p; t++){
//         betaprop[t] = R::rnorm(aux4(l,t),sdpropbeta)*indica(l>t) + RcppTN::rtn1(aux4(l,t),sdpropbeta2,0,R_PosInf)*indica(l==t);
//       }
//       // Rcout << "Valor de betaprop " << betaprop[0] << betaprop[1] << "\n";
//       // Rcout << "Valor de j" << l << "\n";
// 
//       A1 = logcondcompbetajCpp(f,betaprop,alfa(l,_),sigma2[l],y,l+1);
//       A2 = logcondcompbetajCpp(f,aux4(l,_),alfa(l,_),sigma2[l],y,l+1);
//       // Rcout << "Valor de A1 " << logcondcompbetajCpp(f,betaprop,alfa(l,_),sigma2[l],y,l+1) << "\n";
//       // Rcout << "Valor de A2 " << A2 << "\n";
// 
//       // A1 = logcondcompbetajCpp(f[it-1],betaprop,alfa(l,_),auxls[l],y,l);
//       // A2 = logcondcompbetajCpp(f[it-1],aux4(l,_),alfa(l,_),auxls[l],y,l);
//       A  = std::exp(A1 - A2);
//       // Rcout << "Valor de A " << A << "\n";
//       auxvec = {1.0,A};
//       pa = min(auxvec);
//       // Rcout << "Prob aceita beta " << pa << "\n";
//       u = R::runif(0.0,1.0);
//       if(u<pa){
//         auxlb(l,_) = betaprop;
//         aux4(l,_) = betaprop;
//       }else{
//         auxlb(l,_) = aux4(l,_);
//       }
//     }
//     beta[it] = clone(auxlb);
//     // for(int l=0; l <= j-1; l++){
//     //   sigma2prop = RcppTN::rtn1(auxls[l],sdpropsigma2,0,R_PosInf);
//     //   A3 = logcondcompsigma2jCpp(f[it-1],auxlb(l,_),alfa(l,_),sigma2prop,y,l);
//     //   A4 = logcondcompsigma2jCpp(f[it-1],auxlb(l,_),alfa(l,_),auxls[l],y,l);
//     //   B  = std::exp(A3 - A4);
//     //   Rcout << "Valor de B " << B << "\n";
//     //   auxvec = {1.0,B};
//     //   pa = min(auxvec);
//     //   Rcout << "Prob aceita sigma2 " << pa << "\n";
//     //   u = R::runif(0.0,1.0);
//     //   if(u<pa){
//     //     auxls[l] = sigma2prop;
//     //   }
//     // }
//     // sigma2[it] = clone(auxls);
//     // for(int i = 0; i <= n-1; i++){
//     //   fprop = aux5(_,i) + Rcpp::rnorm(p)*sdpropf;
//     //   A5 = logcondcompfiCpp(fprop,beta[it],alfa,sigma2[it],y,i);
//     //   A6 = logcondcompfiCpp(aux5(_,i),beta[it],alfa,sigma2[it],y,i);
//     //   C  = std::exp(A5 - A6);
//     //   // Rcout << "Valor de C " << C << "\n";
//     //   auxvec = {1.0,C};
//     //   pa = min(auxvec);
//     //   // Rcout << "Prob aceita f " << pa << "\n";
//     //   u = R::runif(0.0,1.0);
//     //   if(u<pa){
//     //     auxlf(_,i) = fprop;
//     //     aux5(_,i) = fprop;
//     //   }else{
//     //     auxlf(_,i) = aux5(_,i);
//     //   }
//     // }
//     // f[it] = clone(auxlf);
//     // Rcout << "Iteracao : " << it << "\n";
//   }
//   // List res = List::create(_["beta"] = beta , _["f"] = f, _["sigma2"]=sigma2);
//   // List res = List::create(_["beta"] = beta);
//   return beta;
// }


