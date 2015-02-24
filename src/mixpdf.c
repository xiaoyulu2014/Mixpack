#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 /*#include <R.h>*/
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


void sum_array(double *a,const int *num_elements, double *sum);
void mixpdf(const int *K, const int *nl,double *pi, double *x, double *mu, double *Sigma, double *pdf);
void gsl_matrix_set_zero (gsl_matrix * m);
gsl_matrix * gsl_matrix_alloc (size_t n1, size_t n2);
void gsl_matrix_set (gsl_matrix * m, size_t i, size_t j, double x);
double gsl_linalg_LU_det (gsl_matrix * LU, int signum);  /*determinant*/
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);

void sum_array(double *a, const int *num_elements, double *sum){
    
   for (int i=0; i < *num_elements; ++i){
   *sum += *(a+i);
  }
}

void mixpdf(const int *K, const int *nl, double *pi, double *x, double *mu, double *Sigma, double *pdf){
 
   for (int i=0; i< *K; i++) {
     
     gsl_vector *V = gsl_vector_alloc(*nl);
     const gsl_vector *X = gsl_vector_alloc(*nl);
     gsl_matrix *M = gsl_matrix_alloc (*nl, *nl);
     
     
     for (int k=0; k< *nl; k++){
       double tmp0,tmp1;
       tmp1 = *(x+k);
       tmp0 = *(mu + i*(*nl) + k);
       gsl_vector_set(V,k,tmp0);
       gsl_vector_set(X,k,tmp1);

       for (int j=0; j< *nl; j++){
          double tmp;
          tmp = *(Sigma + i*(*nl)*(*nl) + k*(*nl) + j);
          gsl_matrix_set (M, k, j, tmp);
       }
     }
     int n;
     n = *nl;
     
     *(pdf+i) = *(pi+i) * dmvnorm( n, X, V, M);
   }
    
    /* *(pdf+i) = *(pi+i) * gsl_ran_bivariate_gaussian_pdf (*x-*(mu + i*(*nl)), *(x+1)-*(mu + i*(*nl)+1), sqrt(*(Sigma+i*(*nl)*(*nl))),
                     sqrt(*(Sigma+i*(*nl)*(*nl)+3)),(double) *(Sigma+i*(*nl)*(*nl)+1)/(double) sqrt(( *(Sigma+i*(*nl)*(*nl)) * *(Sigma+i*(*nl)*(*nl)+3)))) ;
                                } */
    
  
  double sum1=0;
  double *sum;
  sum = &sum1;
  sum_array(pdf, K, sum);

  for (int i=0; i< *K; i++){
    *(pdf+i) = (double) *(pdf+i)/ (double) *sum;
  }
 
}

double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
/* multivariate normal density function    */
/*
*  n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*/
int s;
double ax,ay;
gsl_vector *ym, *xm;
gsl_matrix *work = gsl_matrix_alloc(n,n), 
           *winv = gsl_matrix_alloc(n,n);
gsl_permutation *p = gsl_permutation_alloc(n);

gsl_matrix_memcpy( work, var );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(n);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, mean );
ym = gsl_vector_alloc(n);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);
ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

return ay;
}


/*
int main(){
    
     
    
    const int a = 1;
    const int *K;
    K = &a;
    double b = 1;
    double *pi;    
    pi = &b;
    const int c = 1;
    const int *nl;
    nl = &c;
    double out1 = 0;
    double *out;
    out = &out1;
    

    mixpdf(K,nl,pi,out);
    printf("the result is %f\n", *out);

    return(0);
}









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
 
  
