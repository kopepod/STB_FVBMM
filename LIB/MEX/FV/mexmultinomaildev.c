#include "mex.h"
#include <math.h> 
#include <float.h>

void mexmultibinomialdev(uint64_T * x, double * m , double * dev, mwSize D){
  int i,j,offset;
  double a,b;
  uint64_T y;

  for(i = 0; i < D; i++){
    y = x[i];
    offset = 64 * i;
    for (j = 0; j < 64; j++) {
      if((y >> j) & true){
        a = 1;
        b = m[j + offset] + DBL_MIN;
      }
      else
      {
        a = -1;
        b = 1-m[j + offset] - DBL_MIN;
      }
      dev[j + offset] = a/b;
    }  
  }

}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  uint64_T * x   = (uint64_T *) mxGetPr(prhs[0]);
  double * m = mxGetPr(prhs[1]); /*m means*/
  size_t D = mxGetN(prhs[0]); /*D unit64 dims*/
  plhs[0] = mxCreateDoubleMatrix((mwSize)(64 * D),1,mxREAL);
  double * dev = mxGetPr(plhs[0]);
  mexmultibinomialdev(x, m, dev, (mwSize)D);
}

