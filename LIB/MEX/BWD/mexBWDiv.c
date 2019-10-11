#include "mex.h"
#include <stdint.h>
#include <math.h> 

#define X UINT16_C(32)
#define Y UINT16_C(32)
#define T UINT16_C(14)

void estimateorientation(double * v){

  uint_fast16_t x,y,t,idx0,idx1,n;
  uint_fast8_t orientation,j;
   
  double ip,dmin,d,a;
  double cx = 0;
  double cy = 0;
  double m_ = 0;
  
  double binx[4] = {1,0,-1,0};
  double biny[4] = {0,1,0,-1};
  
  for(x = 0; x < X; x++){
    for(y = 0; y < Y; y++){
      ip = 0;
      for(t = 0; t < 7; t++){
        idx0 = y + x*Y + t*X*Y;
        ip += v[idx0];
      }
      m_ += ip;
      cx += (x+0.5) * ip;
      cy += (y+0.5) * ip;
    }
  }
  
  cx = cx/m_-16;
  cy = cy/m_-16;
    
  m_ = sqrt(cx * cx + cy * cy);
  
  cx = cx/m_;
  cy = cy/m_;
  
  d = 0;
  orientation = 0;
  dmin = 1000;
  
  for(j = 0; j < 4; j++){
    d = fabs(binx[j]-cx) + fabs(biny[j]-cy);
    if(d < dmin){
      dmin = d;
      orientation = j;
    }
  }
  n = 0;
	switch(orientation){
		case 0:
			for(t = 0; t < T; t++){
				for(x = 0; x < X/2; x++){
			    for(y = 0; y < Y; y++){
				    idx0 = y + x*Y + t*X*Y;
				    idx1 = y + (X-x-1)*Y + t*X*Y;
			    	m_ = v[idx0];
			    	v[idx0] = v[idx1];
			    	v[idx1] = m_;
    		  }
    		}
  		}
			return;
		case 1:
			for(t = 0; t < T; t++){
				for(x = 0; x < X; x++){
			    for(y = 0; y < Y/2; y++){
				    idx0 = y + x*Y + t*X*Y;
				    idx1 = (Y-y-1) + x*Y + t*X*Y;
			    	m_ = v[idx0];
			    	v[idx0] = v[idx1];
			    	v[idx1] = m_;
    		  }
    		}
  		}
  		return;
		default:
			return;
	}
	
}

void initmap(uint64_T * map){
  uint_fast8_t j;
  uint64_T x = 1;
  map[63] = 1;
  for(j = 1; j < 64; j++){
    x *= 2;
    map[63-j] = x;
  }
}

void integralvideo(double * vx, double * vy){
/* integral video : iv*/

  uint_fast16_t x,y,t,idx0,idx1;
  double accx = 0;
  double accy = 0;
  
  for(t = 0; t < T; t++){
		for(x = 0; x < X; x++){
			accx = 0;
			accy = 0;
			for(y = 0; y < Y; y++){
  			idx0 = y + x*Y + t*X*Y;
   			accx = accx+vx[idx0];
   			accy = accy+vy[idx0];
   			if(x > 0){
   				idx1 = y + (x-1)*Y + t*X*Y;
    			vx[idx0] = accx + vx[idx1];
    			vy[idx0] = accy + vy[idx1];
    		}
    		else{
   				vx[idx0] = accx;
   				vy[idx0] = accy;
   			}
   		}
		}

 		if(t > 0){
			for(x = 0; x < X; x++){
				for(y = 0; y < Y; y++){
	 		 	 	idx0 = y + x*Y + t*X*Y;
	  		 	idx1 = y + x*Y + (t-1)*X*Y;
	  		 	accx = vx[idx0];
	  		 	accy = vy[idx0];
  		   	vx[idx0] = accx + vx[idx1];
  		   	vy[idx0] = accy + vy[idx1];
				}
			}
		}
	}
}

void mexBWDiv(
	double * vx,
	double * vy,
	uint_fast16_t * regions,
	uint_fast16_t N,
	uint64_T * map,
	uint64_T * code){
	
	uint64_T xcode = 0;
	uint64_T ycode = 0;
  uint_fast16_t k,n,idx,x,y,t;
  
  double sax,sbx,say,sby;
  n = 0;
	
	for(k = 0; k < 64; k++){
		sax = 0;say = 0;sbx = 0;sby = 0;
    while(regions[n] != 0){
     y = regions[n+N];
     x = regions[n+2*N];
     t = regions[n+3*N];
     idx = y + x*Y + t*X*Y;
     if(regions[n] == 1){
      sax += vx[idx];
      say += vy[idx];
     }
     else{
      sax -= vx[idx];
      say -= vy[idx];
     }
     n++;
    }
    n++;
	  
    while(regions[n] != 0){
     y = regions[n+N];
     x = regions[n+2*N];
     t = regions[n+3*N];
     idx = y + x*Y + t*X*Y;
     if(regions[n] == 1){
      sbx += vx[idx];
      sby += vy[idx];
     }
     else{
      sbx -= vx[idx];
      sby -= vy[idx];
     }
     n++;
    }
    n++;
    
    if(sax > sbx){
	    xcode += map[k];
    }
    
		if(say > sby){
	    ycode += map[k];
    }	
	}
	
  code[0] = xcode;
  code[1] = ycode;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  double * vx = mxGetPr(prhs[0]);
  double * vy = mxGetPr(prhs[1]);
  uint_fast64_t * regions = (uint_fast64_t *)mxGetPr(prhs[2]);
  uint_fast16_t N = (uint_fast16_t)mxGetM(prhs[2]); /*Number of regions*/
  uint64_T map[64];
  initmap(map);
  
	estimateorientation(vx);
	estimateorientation(vy);
	
	integralvideo(vx,vy);
	      
  plhs[0] = mxCreateNumericMatrix(1,2,mxUINT64_CLASS,mxREAL);
  uint64_T * code = (uint64_T *)mxGetPr(plhs[0]);
    
	mexBWDiv(vx,vy,regions,N,map,code);
    
}
