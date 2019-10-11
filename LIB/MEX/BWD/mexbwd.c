#include "mex.h"
#include <math.h>

void initmap(uint64_T * map){
  int j;
  uint64_T x = 1;
  map[0] = 1;
  for(j = 1; j < 64; j++){
    x *= 2;
    map[j] = x;
  }
}

void mexbwd(double * sv, bool * A, bool * B,mwSize N, mwSize p, uint64_T * map,uint64_T* code){
  int j,k,n;
  double sa,sb,ma,mb;
  n = 0;
  for(k = 0; k < p; k++){
    sa = 0;sb = 0;ma = 0; mb =0;
    for(j = 0; j < N; j++){      
      if(A[n]){
        sa += sv[j];
        ma += 1;
      }
      if(B[n]){
        sb += sv[j];
        mb += 1;
      }
      n += 1;
    }
    if(sa/ma > sb/mb){
      code[0] += map[63-k];
    }
  }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  double * sv = mxGetPr(prhs[0]);
  bool * A = mxGetLogicals(prhs[1]);
  bool * B = mxGetLogicals(prhs[2]);
  size_t N = mxGetM(prhs[0]); /*Number of pixels intensities*/
  size_t p = mxGetN(prhs[1]); /*Number of receptive pairs*/
  
  uint64_T map[64];
  initmap(map);
      
  plhs[0] = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
  uint64_T * code = (uint64_T *)mxGetPr(plhs[0]);
    
  mexbwd(sv,A,B,(mwSize) N,(mwSize) p,map, code);
    
}
