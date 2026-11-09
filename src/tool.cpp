// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



using namespace arma;


// [[Rcpp::export]]
mat resampling(mat X, int n, int p, int B, int tau, int boot_type = 1) {
  // You can include R code blocks in C++ files processed with sourceCpp
  // (useful for testing and development). The R code will be automatically 
  // run after the compilation.
  // boot_type 1 means Rademacher, 2 means normal, 3 means mammen.
  const double sqrt5 = std::sqrt(5.0);
  const double mammen_low  = (1.0 - sqrt5) / 2.0;
  const double mammen_high = (1.0 + sqrt5) / 2.0;
  const double prob_low = (sqrt5 + 1.0) / (2.0 * sqrt5);
  double constant = sqrt(2)/double(n-1);
  double Gn_1e, sigma_n1e;
  mat Y, tmp, XXt;
  mat Hn_B(tau, B);
  vec Hne(tau, fill::zeros);
  vec rand_unif(n, fill::zeros);
  vec et(n);
  XXt = X * X.t();
  XXt.diag().zeros();
  double sigma_n1_const =
    constant * accu(pow(XXt,2));
  
  for(int i = 0; i < B; i++){
    if (boot_type == 1){
      rand_unif.randu();
      et.fill(-1.0);
      // vec et(n, fill::value(-1.0));
      et.elem(find(rand_unif > 0.5)) += 2;
    }
    else if (boot_type == 2){
      et.randn();
    }
    else if (boot_type == 3){
      rand_unif.randu();
      et.fill(mammen_high);
      et.elem(find(rand_unif <= prob_low)).fill(mammen_low);
    }
    //X.print();
    //et.print();
    // Y = et % X.each_col();
    // tmp = Y * Y.t();
    // tmp.diag().zeros();
    tmp = XXt % (et * et.t());
    if (boot_type == 1) {
      // Rademacher: et^2 is identically equal to 1, the denominator does not change with the bootstrap.
      sigma_n1e = sigma_n1_const;
    } else {
      // Normal and Mammen: et^2 is not identically equal to 1
      sigma_n1e = constant * accu(pow(tmp, 2));
    }
    Hne.zeros();
    for(int lag=1; lag<=tau; lag++){
      uvec indices1 = linspace<uvec>(n - lag, n-1, lag);
      uvec indices2 = linspace<uvec>(0, n - 1- lag, n - lag);
      uvec indices = join_cols(indices1, indices2);
      // if(i==0)
      //   indices.print();
      // Gn_1e = accu(tmp.submat(0 , 0, n-1-lag, n-1-lag ) % tmp.submat(lag , lag, n-1, n-1 ));
      Gn_1e = accu(tmp % tmp.submat(indices , indices));
      // Rcpp::Rcout << Gn_1e <<'\n'<< std::endl;
      if(lag > 1){
        Hne(lag-1) = Hne(lag-2) + Gn_1e/sigma_n1e;
      }
      else{
        Hne(lag-1) = Gn_1e / sigma_n1e;
      }
      
      // Hne(lag-1) = sum(Hne) + Gn_1e / sigma_n1e;
    }
    Hn_B.col(i) = Hne;
  }
  return sort(Hn_B, "ascend", 1);
}
// // [[Rcpp::export]]
// mat test(mat tmp, int n, int lag) {
//   uvec indices1 = linspace<uvec>(n - lag, n-1, lag);
//   uvec indices2 = linspace<uvec>(0, n - 1- lag, n - lag);
//   uvec indices = join_cols(indices1, indices2);
//   // if(i==0)
//   //   indices.print();
//   
//   return(tmp.submat(indices,indices));
// }




