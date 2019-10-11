#include "mex.h"
#include <math.h> 

void mexmeanbitwise64(uint64_T * x, double * gamma, double * mean, mwSize D, mwSize T){
  int k,i,j,offset;
  uint64_T y;
  
  for(k = 0; k < T; k++){
    for(i = 0; i < D;i++){
      y = x[k + i*T];
      offset = i * 64;
      for (j = 0; j < 64; j++){
        if((y >> j) & 1){
          mean[j + offset] += gamma[k]; 
        }
      }
    }
  }
    
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  uint64_T * x   = (uint64_T *) mxGetPr(prhs[0]);
  double * gamma = mxGetPr(prhs[1]); /*gamma projections*/
  size_t D = mxGetN(prhs[0]); /*D unit64 dims*/
  size_t T = mxGetM(prhs[0]); /*T number of Features*/
  plhs[0] = mxCreateDoubleMatrix((mwSize)(64 * D),1,mxREAL);
  double * mean = mxGetPr(plhs[0]);
  mexmeanbitwise64(x, gamma, mean, (mwSize)D, (mwSize)T);
}

