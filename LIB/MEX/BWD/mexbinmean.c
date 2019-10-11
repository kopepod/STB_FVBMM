#include "mex.h"

void addLogicalOnes(double * acc, uint64_T x, int offset){
  int j;
  for (j = 0; j < 64; j++) {
    acc[j+offset] += ((x >> j) & 1);
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

void mean(uint64_T * X, uint64_T * z, int N, int D){

  int n,d,j;
  double acc[64*D];
	uint64_T map[64];
	uint64_T x;
  
  initmap(map);
  
	for(d = 0; d < D; d++){
	  for (j = 0; j < 64; j++) {
	  	acc[j + 64*d] = 0;
		}
	}
      
  for(n = 0; n < N; n++){
    for(d = 0; d < D; d++){
    	addLogicalOnes(acc, X[n+d*N],64*d);
		}
  }
  
  for(d = 0; d < D; d++){
  	x = 0;
	  for (j = 0; j < 64; j++) {
	  	/*printf("%d %f\n",j + 64*d,acc[j + 64*d]);*/
	  	if(acc[j + 64*d]/N > 0.5){
	  		x += map[j];
	  	}
		}
		z[d] = x;
	}
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){

  uint64_T * X  = (uint64_T *) mxGetPr(prhs[0]);
  int N = mxGetM(prhs[0]);
  int D = mxGetN(prhs[0]);
  
	plhs[0] = mxCreateNumericMatrix(1, (mwSize)D,mxUINT64_CLASS,mxREAL);
  uint64_T * z = (uint64_T *)mxGetPr(plhs[0]);

	mean(X,z,N,D);
  
}


