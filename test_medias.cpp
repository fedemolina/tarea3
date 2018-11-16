#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <Rmath.h>


// [[Rcpp::export]]
double test_medias(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  double sumX=0.0, sumY=0.0, mediaX=0.0, mediaY=0.0, 
         varX=0.0, varY=0.0, sp=0.0, tstat=0.0,
         n = y.size(),m = x.size();
  
  int i=0;
  for(i = 0; i < m; ++i) {
    sumX = x[i] + sumX;
  }
  mediaX = sumX/m;
  
  for(i = 0; i < m; ++i){
    varX = (x[i] - mediaX)*(x[i]-mediaX) + varX;
  }
  varX = varX/(m-1);
  
  for(i = 0; i < m; ++i) {
    sumY = y[i] + sumY;
  }
  mediaY = sumY/n;
  
  for(i = 0; i < n; ++i){
    varY = (y[i] - mediaY)*(y[i] - mediaY) + varY;
  }
  varY=varY/(n-1);
  
  sp = sqrt((((m-1)*varX) + ((n-1)*varY)) / (m+n-2));
  tstat = (mediaX - mediaY) / (sp*sqrt(1/n + 1/m)); //Hasta acá perfecto. Dan exactamente lo mismo
  return 2*(1- R::pt(std::abs(tstat), (n+m-2), 1, 0)); // Acá se genera una diferencia por culpa del abs(). Uso el abs de std y listo.       
  // retorna el p-valor
}
