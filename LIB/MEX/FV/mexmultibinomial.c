#include "mex.h"
#include <math.h> 

#define W_LENGTH UINT8_C(64)

double binopdf(uint64_T * x, double * m, int k, int D){
  int i,j,offset;
  uint64_T y;
  double p = 1;
  for(i = 0; i < D; i++){
    y = x[i];
    offset = 64*i + k*D*64;
    for (j = 0; j < 64; j++) {
      if((y >> j) & true){
        p *= m[j + offset];
      }
      else
      {
        p *= 1-m[j + offset];
      }
    }
  }
  /*printf("\n");*/
  return p;
}

void mexmultibinomial(uint64_T * x, double * w, double * m , double * gamma, mwSize D, mwSize N){
  int k;
  double s = 0;
  for(k = 0; k < N; k++){
    gamma[k] = w[k] * binopdf(x,m,k,D);
    s += gamma[k]; 
  }

  for(k = 0; k < N; k++){
    gamma[k] = gamma[k]/s;
  }
    
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  uint64_T * x   = (uint64_T *) mxGetPr(prhs[0]);
  double * w = mxGetPr(prhs[1]); /*w weights*/
  double * m = mxGetPr(prhs[2]); /*m means*/
  size_t D = mxGetN(prhs[0]); /*D unit64 dims*/
  size_t N = mxGetM(prhs[1]); /*N number of Gaussians*/
  plhs[0] = mxCreateDoubleMatrix((mwSize)N,1,mxREAL);
  double * gamma = mxGetPr(plhs[0]);
  mexmultibinomial(x, w, m, gamma, (mwSize)D, (mwSize)N);
}

