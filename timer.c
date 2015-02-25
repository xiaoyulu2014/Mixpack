#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <omp.h> 
#include <time.h>


void sum_array(double *a,const int *num_elements, double *sum);
void mixpdf(const int *nr, const int *K, const int *nl, double *pi, double *x, double *mu, double *Sigma,double *assignment);
void gsl_matrix_set_zero (gsl_matrix * m);
gsl_matrix * gsl_matrix_alloc (size_t n1, size_t n2);
void gsl_matrix_set (gsl_matrix * m, size_t i, size_t j, double x);
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);
size_t gsl_vector_max_index (const gsl_vector * v);
void gsl_vector_free (gsl_vector * v);
void gsl_matrix_set_identity (gsl_matrix * m);


int main(){
 
    const int nr1 = 100000, K1 = 10, nl1 = 10, *nr, *K, *nl;
    nr = &nr1;
    K = &K1;
    nl = &nl1;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

    double *pi;   
    pi = malloc((*nl)*(*K)*sizeof(double));
    for (int i=0; i < *K; i++){
      *(pi+i) = (double) 1.0/((double) (*K));
    }
        
    double *x;
    x = malloc((*nl)*(*nr)*sizeof(double));
    for (int i=0;i< (*nr)*(*nl); i++){
      *(x+i) = gsl_rng_uniform (r);
    }
    
    double *mu;
    mu = malloc((*nl)*(*K)*sizeof(double));
    for (int i=0; i< (*K)*(*nl); i++){
      *(mu+i) = gsl_rng_uniform (r);
    }

    double *Sigma;
    Sigma = malloc((*nl)*(*nl)*(*K)*sizeof(double));
    for (int i=0; i< (*K)*(*nl)*(*nl); i++){
      *(Sigma + i) = 0;
    } 
    for (int k=0; k < *K; k++){
      for (int i=0; i< (*nl); i++){
        *(Sigma + k*(*nl)*(*nl)+i*(*nl)+ i) = 1;
      }
    }
    
    double *assignment;   
    clock_t t;
    t = clock();
    mixpdf(nr,K,nl,pi,x,mu,Sigma,assignment); 
    t = clock() - t;
    printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    return 0;  
}

void mixpdf(const int *nr, const int *K, const int *nl, double *pi, double *x, double *mu, double *Sigma,double *assignment){
  
  assignment = malloc((*nr)*sizeof(double));
    #pragma omp parallel for 

  for (int m=0; m < *nr; m++){
  
    gsl_vector *X = gsl_vector_alloc(*nl);
    gsl_vector *pdf = gsl_vector_alloc(*K);


    for (int i=0; i< *nl; i++){
      double tmp1;
      tmp1 = *(x + m*(*nl) + i);
      gsl_vector_set(X,i,tmp1);    
    }  
  

   for (int i=0; i< *K; i++) {
     
     gsl_vector *V = gsl_vector_alloc(*nl);
     gsl_matrix *M = gsl_matrix_alloc (*nl, *nl);
     
     
     for (int k=0; k < *nl; k++){
       double tmp0;
       tmp0 = *(mu + i*(*nl) + k);
       gsl_vector_set(V,k,tmp0);
     }
     
     for (int k=0; k< *nl; k++){
      for (int j=0; j< *nl; j++){
          double tmp;
          tmp = *(Sigma + i*(*nl)*(*nl) + k*(*nl) + j);
          gsl_matrix_set (M, k, j, tmp);
       }
     } 
     int n;
     n = *nl;
     
     gsl_vector_set(pdf,i,*(pi+i) * dmvnorm( n, X, V, M)); 
     gsl_vector_free (V);
     gsl_matrix_free (M); 
     }  
     
     gsl_vector_free (X);    
    *(assignment+m) = gsl_vector_max_index (pdf);   
     gsl_vector_free (pdf);  
    
  }  
  
}

double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){

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
