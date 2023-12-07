#include <Rcpp.h>
#include <iostream>
using namespace std;
using namespace Rcpp;

//' @title fibo
//' @description fibo sampler using Rcpp
//' @param n the number of samples
//' @return fibo sampler using Rcpp \code{n}
//' @examples
//' \dontrun{
//' a <- fibo(5)
//' }
//' @export
// [[Rcpp::export]]

int fibo(int n) {
    int a0 = 0;
    int a1 = 1;
    int an = a1;
    switch (n)
    {
    case 0:
        an = 0;
        break;
    case 1:
        an = 1;
        break;
    default:
        for (int i = 2; i <= n; i++)
        {
            a0 = a1;
            a1 = an;
            an = a0 + a1;
        }
        break;
    }

    return an;
}


//' @title rcppc is isprime
//' @description isprime sampler using Rcpp
//' @param x the number of samples
//' @return is or not rpime using Rcpp \code{n}
//' @examples
//' \dontrun{
//' b <- rcppc(5)
//' }
//' @export
// [[Rcpp::export]]

int rcppc(int x){
  if(x<2||(x & 1==0))
    return 0;
  for (int i=2;i<=sqrt(x);i++)
  {
    if(x%i==0)
      return 0;
  }
  return 1;
}
