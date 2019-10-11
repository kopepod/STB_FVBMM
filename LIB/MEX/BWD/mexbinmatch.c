#include "mex.h"
#include <math.h> 

int logical_ones(uint64_T x){
  int j;
  int n = 0;
  for (j = 0; j < 64; j++) {
    n += ((x >> j) & 1);
  }
  return n;
}

void mexbinmatch(uint64_T * x, uint64_T * y, double * lbl, mwSize n, mwSize m,mwSize p){
  int k,j,i,max_,n_ones,c;
  for(j = 0; j < m; j++){
    max_ = 1000;
    for(k = 0; k < n; k++){
       n_ones = 0;
       for(i = 0;i < p; i++){
          n_ones += logical_ones(x[j + i * m] ^ y[k + i * n]);
       }
      if(n_ones < max_){
        max_ = n_ones;
        c = k;
       }
    }
    lbl[j] = c+1;
    lbl[j+m] = max_;
  }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  uint64_T * x   = (uint64_T *) mxGetPr(prhs[0]);/*features*/
  uint64_T * y   = (uint64_T *) mxGetPr(prhs[1]);/*BOV*/
  size_t m = mxGetM(prhs[0]);
  size_t n = mxGetM(prhs[1]);
  size_t p = mxGetN(prhs[1]);
    
  plhs[0] = mxCreateDoubleMatrix((mwSize)m,2,mxREAL);
  double * lbl = mxGetPr(plhs[0]);
  mexbinmatch(x, y, lbl, (mwSize)n, (mwSize)m, (mwSize)p);
}

