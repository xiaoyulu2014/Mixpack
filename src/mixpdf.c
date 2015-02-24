#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <R.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


void sum_array(double *a,const int *num_elements, double *sum);
void mixpdf(const int *K, const int *nl, double *x, double *pi, double *mu, double *Sigma, double *pdf);


void sum_array(double *a, const int *num_elements, double *sum){
    
   for (int i=0; i < *num_elements; ++i){
   *sum += *(a+i);
  }
}
 
void mixpdf(const int *K, const int *nl, double *x, double *pi, double *mu, double *Sigma, double*pdf){
  
  double *sum;
 
   for (int i=0; i< *K; ++i) {
    *(pdf+i) = *(pi+i) * gsl_ran_bivariate_gaussian_pdf (*(mu + i*(*nl) +1), *(mu + i*(*nl) +2), 
                          *(Sigma+i*(*nl)*(*nl)+1*(*nl)+1),*(Sigma+i*(*nl)*(*nl)+2*(*nl)+2),*(Sigma+i*(*nl)*(*nl)+1*(*nl)+2));
  }
  
  sum_array(pdf, K, sum);

  for (int i=0; i< *K; ++i){
    *(pdf+i) = (double) *(pdf+i)/ (double) *sum;
  }
 
}

/*




int main(){
 arrays: double *b,  b[i]= *(b+i)
 matrices: double *B, B[i][j] = *(B + i*m +j)   (B nxm)
 3-d array: B[i][j][p] = *(B+i*m*k+j*k+p)   (B n*m*k)

double x[n][nl], double pi[K],double mu[K],double Sigma[K]
  double out[K];  
  const unsigned int n=100;
const unsigned int nl=2;
const unsigned int K=3;

  
  for (int i=1;i<K+1;i++) {
    out[i] = pi[i] * gsl_ran_bivariate_gaussian_pdf (mu[i][1], mu[i][1], Sigma[i][1][1], Sigma[i][2][2], Sigma[i][1][2]); 
    printf("prob is %d\n",out[i]);
  }
  
  for (int i=1; i<K+1; i++){
    out[i] = out[i]/sum_array(out,K);
    return(out[i]);
  }
  
  return(0);
}
*/
 
  