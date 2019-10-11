#include "mex.h"
#include <math.h> 
#include <time.h>

int logical_ones(uint64_T x){
  int j;
  int n = 0;
  for (j = 0; j < 64; j++) {
    n += ((x >> j) & 1);
  }
  return n;
}

void randlbl(double * lbl, int m, int K, int seed){
  int k,x;
  srand (seed);
  for(k = 0; k < m; k++){
    x = rand() %K;
    lbl[k] = (double)x;
  }
}

void accumulator(uint64_T x, double * acc, int offset){
  int j;
  for (j = 0; j < 64; j++) {
    acc[j + offset] += (double)((x >> j) & 1);
  }
}

void clearacc(double * acc, int bits){
  int j;
  for(j = 0; j< bits; j++){
    acc[j] = 0;
  }
}

void cparray(double * x, double * ctd, int K){
  int k;
  for(k = 0; k < K; k++){
    x[k] = ctd[k];
  }
}

int isequal(double * x, double * ctd, int K){
  int k,test = 1;
  for(k = 0; k < K; k++){
    if((int)x[k] != (int)ctd[k]){
      test = 0;
      return test;
    }
  }
  return test;
}

void acc2uint64(uint64_T * z, uint64_T * map,double * acc, double sz_cl, int p){
  int k,j;
  uint64_T x;
  for(k = 0; k < p; k++){
    x = 0;
    for(j = 0; j < 64; j++){
      if(acc[j + 64 * k]/sz_cl > 0.5){
        x += map[j];    
      }
    }
    z[k] = x;    
  }
}


void mean(uint64_T * x, uint64_T * ctd, uint64_T * map, double * lbl, int m, int p, int K){
  int k,j,i;
  double acc[64*p];
  double sz_cl;
  uint64_T z[p];
  for(k = 0; k < K; k++){
    sz_cl = 0;
    clearacc(acc, 64 * p);
    for(j = 0; j < m; j++){      
      if(lbl[j] == (double)k){
        sz_cl += 1;        
        for(i = 0; i < p; i++){
          accumulator(x[j + i * m], acc, 64 * i);
        }
      }
    }
    acc2uint64(z, map,acc,sz_cl,p);    
    for(j = 0; j < p; j++){
      ctd[k + j * K] = z[j];
    }
  }
}

void initmap(uint64_T * map){
  int j;
  uint64_T x = 1;
  map[0] = 1;
  for(j = 1; j < 64; j++){
    x *= 2;
    map[j] = x;
    /*printf("%llu\n",x);*/
  }
}


void mexbinkmeans(uint64_T * x, double * lbl, double * sumd, uint64_T * ctd, mwSize N, mwSize D, int K, int seed){
  int k,j,i,max_,n_ones, it;
  uint64_T map[64];
  
  initmap(map);
  
  randlbl(lbl,N,K,seed);
  mean(x,ctd,map,lbl,N,D,K);
  
  
  double c;
  double cplbl[K];
  
  for(it = 0; it < 20; it++){
  
    for(j = 0; j < N; j++){
      max_ = 10000;
      for(k = 0; k < K; k++){
         n_ones = 0;
         for(i = 0;i < D; i++){
            n_ones += logical_ones(x[j + i * N] ^ ctd[k + i * K]);
         }
        if(n_ones < max_){
          max_ = n_ones;
          c = (double)k;
         }     
      }
      lbl[j] = c;
      sumd[j] = max_;      
    }
    
    
    if(isequal(lbl,cplbl,K) == 1){
      /*printf("convergency reached: %d\n",it);*/ 
      return;
    }
    
    /*new centroid*/
    mean(x,ctd,map,lbl,N,D,K);
    
    cparray(cplbl,lbl,K);
  }
  
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  uint64_T * x   = (uint64_T *) mxGetPr(prhs[0]); /*Feature Vector*/
  int K = (int)mxGetScalar(prhs[1]); /*K clusters*/
  
  int seed = time(NULL);
  
  if(nrhs > 2){
    seed = (int)mxGetScalar(prhs[2]);
  }
  
  size_t N = mxGetM(prhs[0]);/* Number of Features*/
  size_t D = mxGetN(prhs[0]);/* Dimensions*/

  plhs[0] = mxCreateDoubleMatrix((mwSize)N,1,mxREAL);
  double * lbl = mxGetPr(plhs[0]);
  
  plhs[1] = mxCreateNumericMatrix((mwSize)K, (mwSize)D,mxUINT64_CLASS,mxREAL);
  uint64_T * ctd = (uint64_T *)mxGetPr(plhs[1]);
  
  plhs[2] = mxCreateDoubleMatrix((mwSize)N,1,mxREAL);
  double * sumd = mxGetPr(plhs[2]);
  
  mexbinkmeans(x, lbl, sumd, ctd, (mwSize)N, (mwSize)D, K, seed);
}


