#include "mex.h"
#include <math.h> 

int logical_ones(uint64_T x, double * V, double q, int n){
  int j;
  for (j = 0; j < 64; j++) {
    if ((x >> j) & 1){
    	V[n] += q;
    }
    n++;
  }
  return n;
}

void binaryvlad(uint64_T * x, uint64_T * ctd, double * qk, double * V, int N, int K, int D){
  int k,j,i,n,m;
  double q,s;
  for(j = 0; j < N; j++){
  	n = 0;
    for(k = 0; k < K; k++){
			q = qk[k];
			for(i = 0; i < D; i++){
				n = logical_ones(x[j + i * N] ^ ctd[k + i * K], V, q, n);
			}
    }
  }
  
  /*L2 intranorm*/

/*
  m = 0;
  
	for(k = 0; k < K; k++){
		s = 0;
		n = m;
		for(i = 0; i < 64 * D; i++){
			s += V[n] * V[n];
			n ++;
		}
		s = sqrt(s);
		n = m;
		for(i = 0; i < 64 * D; i++){
			V[n] /= s;
			n ++;
		}
		m = n;
	}
	*/
	/*L2 full*/
	s = 0;
	for(k = 0; k < K * 64 * D; k++){
		s += V[k];
	}
	/*s = sqrt(s);*/
	for(k = 0; k < K * 64 * D; k++){
		V[k] /= s;
	}
  
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){

  uint64_T * x   = (uint64_T *) mxGetPr(prhs[0]);/*features*/
  uint64_T * ctd = (uint64_T *) mxGetPr(prhs[1]);/*BOV*/
  double * qk = mxGetPr(prhs[2]);/*q cluster weights*/
  int N = (int)mxGetM(prhs[0]);/*N features to map*/
  int K = (int)mxGetM(prhs[1]);/*K clusters*/
  int D = (int)mxGetN(prhs[1]);/*D dimensions*/

  plhs[0] = mxCreateDoubleMatrix(1,K * 64 * D,mxREAL);
  double * V = mxGetPr(plhs[0]);/*VLAD vector*/
  binaryvlad(x, ctd, qk, V, N, K, D);
  
}

