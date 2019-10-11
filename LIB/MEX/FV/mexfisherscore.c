#include "mex.h"
#include <float.h>
#include <stdint.h>

#define BIT_LENGTH UINT32_C(64)

double binopdf(uint_fast64_t * x, double * m, uint_fast32_t k, uint_fast32_t D){
  uint_fast32_t i,j,offset;
  uint_fast64_t y;
  double p = 1;
  uint_fast32_t s_offset = k*D*BIT_LENGTH;
  for(i = 0; i < D; i++){
    y = x[i];
    /*printf("%llu ",y);*/
    offset = BIT_LENGTH*i + s_offset;
    for (j = 0; j < BIT_LENGTH; j++) {
      if((y >> j) & true){
        p *= m[j + offset];
      }
      else{
        p *= 1-m[j + offset];
      }
    }
  }
  /*printf("\n");*/
  return p;
}

void mexmultibinomialdev(uint_fast64_t * x, double * m1, double * m2, double * dev, double g, uint_fast32_t D, uint_fast32_t off)
{
  uint_fast64_t y;
  uint_fast32_t i,j,offset;
  
  for(i = 0; i < D; i++){
    y = x[i];
    offset = BIT_LENGTH * i;
    for (j = 0; j < BIT_LENGTH; j++) {
      if((y >> j) & true){
        dev[j + offset] = g * m1[j + offset + off];
      }
      else{
        dev[j + offset] = g * m2[j + offset + off];
      }
    }
  }
}

void mexmultibinomial(uint_fast64_t * x, double * w, double * m , double * gamma, uint_fast32_t D, uint_fast32_t N){
  uint_fast32_t k;
  double s = 0;
  for(k = 0; k < N; k++){
    gamma[k] = w[k] * binopdf(x,m,k,D);
    s += gamma[k]; 
  }
  for(k = 0; k < N; k++){
    gamma[k] = gamma[k]/s;
  }
}

void mexfisherscore(uint_fast64_t * X64, double * w, double * m, double * m1, double * m2, double * G, uint_fast32_t D, uint_fast32_t N,uint_fast32_t T){
  uint_fast32_t t,d,i,k,offset;
  double y[N];
  uint_fast32_t M = BIT_LENGTH * D;
  double dev[M];
  uint_fast64_t z[D];
      
  for(t = 0; t < T; t++){
    for(d = 0;d < D; d++){
      z[d] = X64[t + d * T];
    }
    mexmultibinomial(z, w, m, y, D, N);
    for(i = 0; i < N; i++){
      offset = i * M;    
      mexmultibinomialdev(z, m1, m2, dev, y[i], D, offset);
      
      for(k = 0;k < M; k++){
        G[k + offset] += dev[k];
      }
    }
  }
  
  for(k = 0;k < M * N; k++){
    G[k] = G[k]/T;
  }
  
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){
  uint_fast64_t * X64 = (uint_fast64_t *) mxGetPr(prhs[0]);
  double * w = mxGetPr(prhs[1]); /*w weights*/
  double * m = mxGetPr(prhs[2]); /*m means*/
  double * m1 = mxGetPr(prhs[3]); /*1/m means*/
  double * m2 = mxGetPr(prhs[4]); /*-1/(1-m) means*/
  uint_fast32_t D = (uint_fast32_t)mxGetN(prhs[0]); /*D unit64 dims*/
  uint_fast32_t T = (uint_fast32_t)mxGetM(prhs[0]); /*T unit64 features*/
  uint_fast32_t N = (uint_fast32_t)mxGetM(prhs[1]); /*N number of Gaussians*/
  plhs[0] = mxCreateDoubleMatrix(1,BIT_LENGTH*D*N,mxREAL);
  double * G = mxGetPr(plhs[0]);
  
  mexfisherscore(X64, w, m, m1,m2,G,D,N,T);
}

