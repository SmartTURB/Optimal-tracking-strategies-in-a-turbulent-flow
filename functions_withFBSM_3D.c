#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include "common_vars.h"

//________________________________________________________________________________________//
/*  Function to allocate a vector instead of a matrix, e.g. A_history[3][3][NITER+1] -> *A_history  */
/* A_history[matrix_index(i,j,t)] gives the equivalent of A_history[i][j][t]. */
int matrix_index(int i, int j, int t){
  return 9*t + 3*i + j; 
}

//________________________________________________________________________________________//
/*   RND generators   box-muller                */
double gauss(void){
  double x,y,r,fac=0;
  while(fac==0){
    x=2.0*drand48()-1.0;
    y=2.0*drand48()-1.0;
    r=x*x+y*y;
    if(r<1.0 && r!=0.0){
      fac = sqrt(-2.0*log(r)/r);
    }
  }
  return x*fac;
}

//________________________________________________________________________________________//
/*     Distances between Lagrangian particles  */
/*     actually we have fixed x[0]=y[0]=z[0]=0 beacuse of the assumption of linearity */
double dist(double *x,double *y, double *z, int i){
  double scra;
  scra=(x[i]-x[0])*(x[i]-x[0])+(y[i]-y[0])*(y[i]-y[0]) + (z[i]-z[0])*(z[i]-z[0]);
  return sqrt(scra);
}

//________________________________________________________________________________________//
/* Function to compute gradients - useful to calculate surfing and perturbative control*/
void gradients(double xx[],double yy[], double zz[], double dux_dx[], double dux_dy[], double dux_dz[], double duy_dx[],double duy_dy[], double duy_dz[], double duz_dx[],double duz_dy[], double duz_dz[],int step){
  int  ip;
  for(ip=0;ip<NPART;ip++){
    dux_dx[ip]=AA_history[matrix_index(0,0,step)];
    dux_dy[ip]=AA_history[matrix_index(0,1,step)];
    dux_dz[ip]=AA_history[matrix_index(0,2,step)];
    duy_dx[ip]=AA_history[matrix_index(1,0,step)];
    duy_dy[ip]=AA_history[matrix_index(1,1,step)];
    duy_dz[ip]=AA_history[matrix_index(1,2,step)];
    duz_dx[ip]=AA_history[matrix_index(2,0,step)];
    duz_dy[ip]=AA_history[matrix_index(2,1,step)];
    duz_dz[ip]=AA_history[matrix_index(2,2,step)];
  }
}


//________________________________________________________________________________________// 
/*  Compute velocities at given position and time   */
double uxx, uyy, uzz; // velocity field
void derivs_ip(double xx,double yy, double zz,  double tensor_A[9]){
  uxx = tensor_A[0]*xx + tensor_A[1]*yy + tensor_A[2]*zz; 
  uyy = tensor_A[3]*xx + tensor_A[4]*yy + tensor_A[5]*zz;
  uzz = tensor_A[6]*xx + tensor_A[7]*yy + tensor_A[8]*zz;
}


//________________________________________________________________________________________//
/* Runge-Kutta to advance particles -  MAIN integration */ 
void rk4_advance_particles(double *x,double *y, double *z, int step, int c[NPART]){
  int ip;
  double k1x[NPART], k1y[NPART], k1z[NPART],k2x[NPART], k2y[NPART], k2z[NPART], k3x[NPART], k3y[NPART], k3z[NPART], k4x[NPART], k4y[NPART], k4z[NPART];
  double scra;
  double tensor_A[9];
  
  tensor_A[0] = AA_history[matrix_index(0,0,step)];
  tensor_A[1] = AA_history[matrix_index(0,1,step)];
  tensor_A[2] = AA_history[matrix_index(0,2,step)];
  tensor_A[3] = AA_history[matrix_index(1,0,step)];
  tensor_A[4] = AA_history[matrix_index(1,1,step)];
  tensor_A[5] = AA_history[matrix_index(1,2,step)];
  tensor_A[6] = AA_history[matrix_index(2,0,step)];
  tensor_A[7] = AA_history[matrix_index(2,1,step)];
  tensor_A[8] = AA_history[matrix_index(2,2,step)];
	
  //step1// 
  for(ip=0;ip<NPART;ip++){      
    derivs_ip(x[ip],y[ip],z[ip], tensor_A); //output uxx  uyy e uzz	  	  
    k1x[ip] = dt*(uxx + vs[ip]*px[ip]);
    k1y[ip] = dt*(uyy + vs[ip]*py[ip]);
    k1z[ip] = dt*(uzz + vs[ip]*pz[ip]);
  }
	
  //step2//
  tensor_A[0] = (AA_history[matrix_index(0,0,step)] + AA_history[matrix_index(0,0,step+1)] )/ 2.;
  tensor_A[1] = (AA_history[matrix_index(0,1,step)] + AA_history[matrix_index(0,1,step+1)] )/ 2.;
  tensor_A[2] = (AA_history[matrix_index(0,2,step)] + AA_history[matrix_index(0,2,step+1)] )/ 2.;
  tensor_A[3] = (AA_history[matrix_index(1,0,step)] + AA_history[matrix_index(1,0,step+1)] )/ 2.;
  tensor_A[4] = (AA_history[matrix_index(1,1,step)] + AA_history[matrix_index(1,1,step+1)] )/ 2.;
  tensor_A[5] = (AA_history[matrix_index(1,2,step)] + AA_history[matrix_index(1,2,step+1)] )/ 2.;
  tensor_A[6] = (AA_history[matrix_index(2,0,step)] + AA_history[matrix_index(2,0,step+1)] )/ 2.;
  tensor_A[7] = (AA_history[matrix_index(2,1,step)] + AA_history[matrix_index(2,1,step+1)] )/ 2.;
  tensor_A[8] = (AA_history[matrix_index(2,2,step)] + AA_history[matrix_index(2,2,step+1)] )/ 2.;
	
  for(ip=0;ip<NPART;ip++){	  
    derivs_ip(x[ip]+k1x[ip]/2.,y[ip]+k1y[ip]/2., z[ip] + k1z[ip]/2., tensor_A);	 
    k2x[ip] = dt*(uxx + vs[ip]*px[ip]);
    k2y[ip] = dt*(uyy + vs[ip]*py[ip]);
    k2z[ip] = dt*(uzz + vs[ip]*pz[ip]);
  }
	
  //step3//
  for(ip=0;ip<NPART;ip++){	
    derivs_ip(x[ip]+k2x[ip]/2,y[ip]+k2y[ip]/2., z[ip] + k2z[ip]/2., tensor_A);
    k3x[ip] = dt*(uxx + vs[ip]*px[ip]);
    k3y[ip] = dt*(uyy + vs[ip]*py[ip]);
    k3z[ip] = dt*(uzz + vs[ip]*pz[ip]);
  }
  //step4//
  tensor_A[0] = AA_history[matrix_index(0,0,step+1)];
  tensor_A[1] = AA_history[matrix_index(0,1,step+1)];
  tensor_A[2] = AA_history[matrix_index(0,2,step+1)];
  tensor_A[3] = AA_history[matrix_index(1,0,step+1)];
  tensor_A[4] = AA_history[matrix_index(1,1,step+1)];
  tensor_A[5] = AA_history[matrix_index(1,2,step+1)];
  tensor_A[6] = AA_history[matrix_index(2,0,step+1)];
  tensor_A[7] = AA_history[matrix_index(2,1,step+1)];
  tensor_A[8] = AA_history[matrix_index(2,2,step+1)];
  for(ip=0;ip<NPART;ip++){	  
    derivs_ip(x[ip]+k3x[ip],y[ip]+k3y[ip],z[ip]+k3z[ip], tensor_A);	  
    k4x[ip] = dt*(uxx + vs[ip]*px[ip]);
    k4y[ip] = dt*(uyy + vs[ip]*py[ip]);
    k4z[ip] = dt*(uzz + vs[ip]*pz[ip]);
  }	
		
  for(ip=INDEX_TRACER; ip<NPART; ip++){                       
    scra = dist(x,y,z,ip);
    if(c[ip]==NO_CAPTURE){ //check if the agent has capture
      if(scra<=capture_distance){                                                                                                         
	c[ip] = OK_CAPTURE;
      }
    }

    //Advance particle only if capture is not achieved, otherwise the episode is considered succeed and the dynamics is stopped.
    if(c[ip]==NO_CAPTURE){ 
      x[ip] += k1x[ip]/6 + k2x[ip]/3 + k3x[ip]/3 + k4x[ip]/6;                                                                         
      y[ip] += k1y[ip]/6 + k2y[ip]/3 + k3y[ip]/3 + k4y[ip]/6;    
      z[ip] += k1z[ip]/6 + k2z[ip]/3 + k3z[ip]/3 + k4z[ip]/6;	                                                          
    }                                                                                                                                 
  }	  
}

//________________________________________________________________________________________//
/* Runge-Kutta to advance the optimal particle in the FBSM algorithm  */
void rk4_forward_R(double Rv_x[], double Rv_y[], double Rv_z[], double pv_x[], double pv_y[], double pv_z[], int ip){
  double k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y,k4z;
  int t;
  double scra, filter;
  double tensor_A[9];
  
  for (t=0;t<NITER;t++){    
    tensor_A[0] = AA_history[matrix_index(0,0,t)];
    tensor_A[1] = AA_history[matrix_index(0,1,t)];
    tensor_A[2] = AA_history[matrix_index(0,2,t)];
    tensor_A[3] = AA_history[matrix_index(1,0,t)];
    tensor_A[4] = AA_history[matrix_index(1,1,t)];
    tensor_A[5] = AA_history[matrix_index(1,2,t)];
    tensor_A[6] = AA_history[matrix_index(2,0,t)];
    tensor_A[7] = AA_history[matrix_index(2,1,t)];
    tensor_A[8] = AA_history[matrix_index(2,2,t)];
    
    //step1//
    derivs_ip(Rv_x[t],Rv_y[t], Rv_z[t], tensor_A); //output ux uy uz
    scra = sqrt(Rv_x[t]*Rv_x[t] + Rv_y[t]*Rv_y[t] + Rv_z[t]*Rv_z[t]);
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k1x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k1y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k1z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    //step2//
    tensor_A[0] = (AA_history[matrix_index(0,0,t)] + AA_history[matrix_index(0,0,t+1)] )/ 2.;
    tensor_A[1] = (AA_history[matrix_index(0,1,t)] + AA_history[matrix_index(0,1,t+1)] )/ 2.;
    tensor_A[2] = (AA_history[matrix_index(0,2,t)] + AA_history[matrix_index(0,2,t+1)] )/ 2.;
    tensor_A[3] = (AA_history[matrix_index(1,0,t)] + AA_history[matrix_index(1,0,t+1)] )/ 2.;
    tensor_A[4] = (AA_history[matrix_index(1,1,t)] + AA_history[matrix_index(1,1,t+1)] )/ 2.;
    tensor_A[5] = (AA_history[matrix_index(1,2,t)] + AA_history[matrix_index(1,2,t+1)] )/ 2.;
    tensor_A[6] = (AA_history[matrix_index(2,0,t)] + AA_history[matrix_index(2,0,t+1)] )/ 2.;
    tensor_A[7] = (AA_history[matrix_index(2,1,t)] + AA_history[matrix_index(2,1,t+1)] )/ 2.;
    tensor_A[8] = (AA_history[matrix_index(2,2,t)] + AA_history[matrix_index(2,2,t+1)] )/ 2.;
		
    derivs_ip(Rv_x[t]+k1x/2.,Rv_y[t]+k1y/2., Rv_z[t]+k1z/2., tensor_A);
    scra = sqrt((Rv_x[t]+ k1x/2)*(Rv_x[t]+ k1x/2)+ (Rv_y[t]+ k1y/2)*(Rv_y[t]+ k1y/2) + (Rv_z[t]+ k1z/2)*(Rv_z[t]+ k1z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k2x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k2y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k2z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    //step3//
    derivs_ip(Rv_x[t]+k2x/2.,Rv_y[t]+k2y/2., Rv_z[t]+k2z/2., tensor_A);
    scra = sqrt((Rv_x[t]+ k2x/2)*(Rv_x[t]+ k2x/2)+ (Rv_y[t]+ k2y/2)*(Rv_y[t]+ k2y/2) +  (Rv_z[t]+ k2z/2)*(Rv_z[t]+ k2z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k3x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k3y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k3z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    //step4//
    tensor_A[0] = AA_history[matrix_index(0,0,t+1)];
    tensor_A[1] = AA_history[matrix_index(0,1,t+1)];
    tensor_A[2] = AA_history[matrix_index(0,2,t+1)];
    tensor_A[3] = AA_history[matrix_index(1,0,t+1)];
    tensor_A[4] = AA_history[matrix_index(1,1,t+1)];
    tensor_A[5] = AA_history[matrix_index(1,2,t+1)];
    tensor_A[6] = AA_history[matrix_index(2,0,t+1)];
    tensor_A[7] = AA_history[matrix_index(2,1,t+1)];
    tensor_A[8] = AA_history[matrix_index(2,2,t+1)];
		
    derivs_ip(Rv_x[t]+k3x,Rv_y[t]+k3y,Rv_z[t]+k3z, tensor_A);
    scra = sqrt((Rv_x[t]+ k3x)*(Rv_x[t]+ k3x)+ (Rv_y[t]+ k3y)*(Rv_y[t]+ k3y) + (Rv_z[t]+ k3z)*(Rv_z[t]+ k3z));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k4x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k4y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k4z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    //integration
    Rv_x[t+1] = Rv_x[t] + (k1x/6 + k2x/3 + k3x/3 + k4x/6);
    Rv_y[t+1] = Rv_y[t] + (k1y/6 + k2y/3 + k3y/3 + k4y/6);
    Rv_z[t+1] = Rv_z[t] + (k1z/6 + k2z/3 + k3z/3 + k4z/6);
		
  }
} 


// Initialization FBSM with PURE PURSUIT
void rk4_forward_R_PP_initialization(double Rv_x[], double Rv_y[], double Rv_z[],double pv_x[], double pv_y[], double pv_z[], int ip){
  double k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y, k4z;
  int t;
  double scra, filter;
  double tensor_A[9];
	
  for (t=0;t<NITER+1;t++){ 
    pv_x[t] = -Rv_x[t];
    pv_y[t] = -Rv_y[t];
    pv_z[t] = -Rv_z[t];
    scra = sqrt(pv_x[t]*pv_x[t] + pv_y[t] * pv_y[t] + pv_z[t]*pv_z[t]);
    pv_x[t]/=scra;
    pv_y[t]/=scra;
    pv_z[t]/=scra;
		
    tensor_A[0] = AA_history[matrix_index(0,0,t)];
    tensor_A[1] = AA_history[matrix_index(0,1,t)];
    tensor_A[2] = AA_history[matrix_index(0,2,t)];
    tensor_A[3] = AA_history[matrix_index(1,0,t)];
    tensor_A[4] = AA_history[matrix_index(1,1,t)];
    tensor_A[5] = AA_history[matrix_index(1,2,t)];
    tensor_A[6] = AA_history[matrix_index(2,0,t)];
    tensor_A[7] = AA_history[matrix_index(2,1,t)];
    tensor_A[8] = AA_history[matrix_index(2,2,t)];
		
    //step1//
    derivs_ip(Rv_x[t],Rv_y[t], Rv_z[t], tensor_A);//output ux e uy
    scra = sqrt(Rv_x[t]*Rv_x[t] + Rv_y[t]*Rv_y[t] + Rv_z[t]*Rv_z[t]);
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k1x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k1y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k1z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    //step2//
    tensor_A[0] = (AA_history[matrix_index(0,0,t)] + AA_history[matrix_index(0,0,t+1)] )/ 2.;
    tensor_A[1] = (AA_history[matrix_index(0,1,t)] + AA_history[matrix_index(0,1,t+1)] )/ 2.;
    tensor_A[2] = (AA_history[matrix_index(0,2,t)] + AA_history[matrix_index(0,2,t+1)] )/ 2.;
    tensor_A[3] = (AA_history[matrix_index(1,0,t)] + AA_history[matrix_index(1,0,t+1)] )/ 2.;
    tensor_A[4] = (AA_history[matrix_index(1,1,t)] + AA_history[matrix_index(1,1,t+1)] )/ 2.;
    tensor_A[5] = (AA_history[matrix_index(1,2,t)] + AA_history[matrix_index(1,2,t+1)] )/ 2.;
    tensor_A[6] = (AA_history[matrix_index(2,0,t)] + AA_history[matrix_index(2,0,t+1)] )/ 2.;
    tensor_A[7] = (AA_history[matrix_index(2,1,t)] + AA_history[matrix_index(2,1,t+1)] )/ 2.;
    tensor_A[8] = (AA_history[matrix_index(2,2,t)] + AA_history[matrix_index(2,2,t+1)] )/ 2.;

    derivs_ip(Rv_x[t]+k1x/2,Rv_y[t]+k1y/2,Rv_z[t]+k1z/2, tensor_A);
    scra = sqrt((Rv_x[t]+ k1x/2)*(Rv_x[t]+ k1x/2)+ (Rv_y[t]+ k1y/2)*(Rv_y[t]+ k1y/2) + (Rv_z[t]+ k1z/2)*(Rv_z[t]+ k1z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k2x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k2y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k2z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    //step3//
    derivs_ip(Rv_x[t]+k2x/2,Rv_y[t]+k2y/2,Rv_z[t]+k2z/2, tensor_A);
    scra = sqrt((Rv_x[t]+ k2x/2)*(Rv_x[t]+ k2x/2)+ (Rv_y[t]+ k2y/2)*(Rv_y[t]+ k2y/2) +  (Rv_z[t]+ k2z/2)*(Rv_z[t]+ k2z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k3x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k3y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k3z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    //step4//
    tensor_A[0] = AA_history[matrix_index(0,0,t+1)];
    tensor_A[1] = AA_history[matrix_index(0,1,t+1)];
    tensor_A[2] = AA_history[matrix_index(0,2,t+1)];
    tensor_A[3] = AA_history[matrix_index(1,0,t+1)];
    tensor_A[4] = AA_history[matrix_index(1,1,t+1)];
    tensor_A[5] = AA_history[matrix_index(1,2,t+1)];
    tensor_A[6] = AA_history[matrix_index(2,0,t+1)];
    tensor_A[7] = AA_history[matrix_index(2,1,t+1)];
    tensor_A[8] = AA_history[matrix_index(2,2,t+1)];
    derivs_ip(Rv_x[t]+k3x,Rv_y[t]+k3y, Rv_z[t]+k3z, tensor_A);
    scra = sqrt((Rv_x[t]+ k3x)*(Rv_x[t]+ k3x)+ (Rv_y[t]+ k3y)*(Rv_y[t]+ k3y) +(Rv_z[t]+ k3z)*(Rv_z[t]+ k3z));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k4x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k4y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k4z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
		
    // integration
    Rv_x[t+1] = Rv_x[t] + (k1x/6 + k2x/3 + k3x/3 + k4x/6);
    Rv_y[t+1] = Rv_y[t] + (k1y/6 + k2y/3 + k3y/3 + k4y/6);
    Rv_z[t+1] = Rv_z[t] + (k1z/6 + k2z/3 + k3z/3 + k4z/6);

  }
}

// function to initialize the FBSM with the surfing control strategy
void control_surfing_forinitialization(double Rx[],double Ry[], double Rz[], int t){
  int i,j;
  double dux_dx[NPART],dux_dy[NPART], dux_dz[NPART],duy_dx[NPART],duy_dy[NPART], duy_dz[NPART], duz_dx[NPART], duz_dy[NPART], duz_dz[NPART]; 
  double R0[3]; 
  double scra;
  double vect_grad[9]; //vector that containts all the elements of the gradients matrix, i.e., grad_matrix[i,j] = vect_grad[i*3+j] 
  double vect_forexpgrad[9] = {0};
	
  gradients(Rx,Ry,Rz,dux_dx,dux_dy,dux_dz,duy_dx,duy_dy, duy_dz, duz_dx, duz_dy, duz_dz, t); 

  scra=sqrt(Rx[t]*Rx[t]+Ry[t]*Ry[t]+Rz[t]*Rz[t]);	
  R0[0] = Rx[t]/scra; 
  R0[1] = Ry[t]/scra;
  R0[2] = Rz[t]/scra;	 

  grad_matrix[0][0] = tau_surfing*dux_dx[INDEX_SURFING];
  grad_matrix[0][1] = tau_surfing*dux_dy[INDEX_SURFING];
  grad_matrix[0][2] = tau_surfing*dux_dz[INDEX_SURFING];
  grad_matrix[1][0] = tau_surfing*duy_dx[INDEX_SURFING];
  grad_matrix[1][1] = tau_surfing*duy_dy[INDEX_SURFING];
  grad_matrix[1][2] = tau_surfing*duy_dz[INDEX_SURFING];
  grad_matrix[2][0] = tau_surfing*duz_dx[INDEX_SURFING];
  grad_matrix[2][1] = tau_surfing*duz_dy[INDEX_SURFING];
  grad_matrix[2][2] = tau_surfing*duz_dz[INDEX_SURFING];
  vect_grad[0] = grad_matrix[0][0];
  vect_grad[1] = grad_matrix[0][1];
  vect_grad[2] = grad_matrix[0][2];
  vect_grad[3] = grad_matrix[1][0];
  vect_grad[4] = grad_matrix[1][1];
  vect_grad[5] = grad_matrix[1][2];
  vect_grad[6] = grad_matrix[2][0];
  vect_grad[7] = grad_matrix[2][1];
  vect_grad[8] = grad_matrix[2][2];
	
  gsl_matrix_view m = gsl_matrix_view_array(vect_grad,3,3);
  gsl_matrix_view em = gsl_matrix_view_array(vect_forexpgrad,3,3);
  gsl_linalg_exponential_ss(&m.matrix,&em.matrix, .001);

  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_grad[i][j] = gsl_matrix_get(&em.matrix, i,j);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_gradT[i][j] = Exp_grad[j][i];

  px[INDEX_SURFING] = -Exp_gradT[0][0]*R0[0] - Exp_gradT[0][1]*R0[1] - Exp_gradT[0][2]*R0[2];
  py[INDEX_SURFING] = -Exp_gradT[1][0]*R0[0] - Exp_gradT[1][1]*R0[1] - Exp_gradT[1][2]*R0[2]; 
  pz[INDEX_SURFING] = -Exp_gradT[2][0]*R0[0] - Exp_gradT[2][1]*R0[1] - Exp_gradT[2][2]*R0[2]; 
	
  scra=sqrt(px[INDEX_SURFING]*px[INDEX_SURFING]+py[INDEX_SURFING]*py[INDEX_SURFING] + pz[INDEX_SURFING]*pz[INDEX_SURFING]);   
  px[INDEX_SURFING]/=scra;  
  py[INDEX_SURFING]/=scra;  
  pz[INDEX_SURFING]/=scra;
					
} //end function control surfing for initialization



// Initialization FBSM with SURFING CONTROL
void rk4_forward_R_SC_initialization(double Rv_x[], double Rv_y[], double Rv_z[],double pv_x[], double pv_y[], double pv_z[], int ip){
  double k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y, k4z;
  int t;
  double scra, filter;
  double tensor_A[9];
	
  for (t=0;t<NITER+1;t++){
    control_surfing_forinitialization(Rv_x,Rv_y,Rv_z,t);
    pv_x[t] = px[INDEX_SURFING];
    pv_y[t] = py[INDEX_SURFING];
    pv_z[t] = pz[INDEX_SURFING];
    tensor_A[0] = AA_history[matrix_index(0,0,t)];
    tensor_A[1] = AA_history[matrix_index(0,1,t)];
    tensor_A[2] = AA_history[matrix_index(0,2,t)];
    tensor_A[3] = AA_history[matrix_index(1,0,t)];
    tensor_A[4] = AA_history[matrix_index(1,1,t)];
    tensor_A[5] = AA_history[matrix_index(1,2,t)];
    tensor_A[6] = AA_history[matrix_index(2,0,t)];
    tensor_A[7] = AA_history[matrix_index(2,1,t)];
    tensor_A[8] = AA_history[matrix_index(2,2,t)];
	  
    //step1//
    derivs_ip(Rv_x[t],Rv_y[t], Rv_z[t], tensor_A);//in output ux e uy
    scra = sqrt(Rv_x[t]*Rv_x[t] + Rv_y[t]*Rv_y[t] + Rv_z[t]*Rv_z[t]);
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k1x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k1y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k1z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //step2//
    tensor_A[0] = (AA_history[matrix_index(0,0,t)] + AA_history[matrix_index(0,0,t+1)] )/ 2.;
    tensor_A[1] = (AA_history[matrix_index(0,1,t)] + AA_history[matrix_index(0,1,t+1)] )/ 2.;
    tensor_A[2] = (AA_history[matrix_index(0,2,t)] + AA_history[matrix_index(0,2,t+1)] )/ 2.;
    tensor_A[3] = (AA_history[matrix_index(1,0,t)] + AA_history[matrix_index(1,0,t+1)] )/ 2.;
    tensor_A[4] = (AA_history[matrix_index(1,1,t)] + AA_history[matrix_index(1,1,t+1)] )/ 2.;
    tensor_A[5] = (AA_history[matrix_index(1,2,t)] + AA_history[matrix_index(1,2,t+1)] )/ 2.;
    tensor_A[6] = (AA_history[matrix_index(2,0,t)] + AA_history[matrix_index(2,0,t+1)] )/ 2.;
    tensor_A[7] = (AA_history[matrix_index(2,1,t)] + AA_history[matrix_index(2,1,t+1)] )/ 2.;
    tensor_A[8] = (AA_history[matrix_index(2,2,t)] + AA_history[matrix_index(2,2,t+1)] )/ 2.;
	  
    derivs_ip(Rv_x[t]+k1x/2,Rv_y[t]+k1y/2,Rv_z[t]+k1z/2, tensor_A);
    scra = sqrt((Rv_x[t]+ k1x/2)*(Rv_x[t]+ k1x/2)+ (Rv_y[t]+ k1y/2)*(Rv_y[t]+ k1y/2) + (Rv_z[t]+ k1z/2)*(Rv_z[t]+ k1z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k2x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k2y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k2z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //step3//
    derivs_ip(Rv_x[t]+k2x/2,Rv_y[t]+k2y/2,Rv_z[t]+k2z/2, tensor_A);
    scra = sqrt((Rv_x[t]+ k2x/2)*(Rv_x[t]+ k2x/2)+ (Rv_y[t]+ k2y/2)*(Rv_y[t]+ k2y/2) +  (Rv_z[t]+ k2z/2)*(Rv_z[t]+ k2z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k3x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k3y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k3z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //step4//
    tensor_A[0] = AA_history[matrix_index(0,0,t+1)];
    tensor_A[1] = AA_history[matrix_index(0,1,t+1)];
    tensor_A[2] = AA_history[matrix_index(0,2,t+1)];
    tensor_A[3] = AA_history[matrix_index(1,0,t+1)];
    tensor_A[4] = AA_history[matrix_index(1,1,t+1)];
    tensor_A[5] = AA_history[matrix_index(1,2,t+1)];
    tensor_A[6] = AA_history[matrix_index(2,0,t+1)];
    tensor_A[7] = AA_history[matrix_index(2,1,t+1)];
    tensor_A[8] = AA_history[matrix_index(2,2,t+1)];
    derivs_ip(Rv_x[t]+k3x,Rv_y[t]+k3y, Rv_z[t]+k3z, tensor_A);
    scra = sqrt((Rv_x[t]+ k3x)*(Rv_x[t]+ k3x)+ (Rv_y[t]+ k3y)*(Rv_y[t]+ k3y) +(Rv_z[t]+ k3z)*(Rv_z[t]+ k3z));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k4x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k4y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k4z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //integration
    Rv_x[t+1] = Rv_x[t] + (k1x/6 + k2x/3 + k3x/3 + k4x/6);
    Rv_y[t+1] = Rv_y[t] + (k1y/6 + k2y/3 + k3y/3 + k4y/6);
    Rv_z[t+1] = Rv_z[t] + (k1z/6 + k2z/3 + k3z/3 + k4z/6);
	  
  }
}

// function to initialize the FBSM with the perturbative optimal control 
void control_perturbative_forinitialization(double Rx[],double Ry[], double Rz[],int t){
  int i,j;
  double dux_dx[NPART],dux_dy[NPART], dux_dz[NPART],duy_dx[NPART],duy_dy[NPART], duy_dz[NPART], duz_dx[NPART], duz_dy[NPART], duz_dz[NPART];
  double R0[NPART][3];
  double RTAU[3];
  double scra;
  double vect_grad[9]; //vector that containts all the elements of the gradients matrix, i.e., grad_matrix[i,j] = vect_grad[i*3+j]
  double vect_forexpgrad[9] = {0};	
	
  gradients(Rx,Ry,Rz,dux_dx,dux_dy,dux_dz,duy_dx,duy_dy, duy_dz, duz_dx, duz_dy, duz_dz, t); 

  scra=sqrt(Rx[t]*Rx[t]+Ry[t]*Ry[t]+Rz[t]*Rz[t]);	
  R0[INDEX_PERTURBATIVE][0] = Rx[t]/scra; 
  R0[INDEX_PERTURBATIVE][1] = Ry[t]/scra;
  R0[INDEX_PERTURBATIVE][2] = Rz[t]/scra;	 
   
  grad_matrix[0][0] = tau_perturbative*dux_dx[INDEX_PERTURBATIVE];
  grad_matrix[0][1] = tau_perturbative*dux_dy[INDEX_PERTURBATIVE];
  grad_matrix[0][2] = tau_perturbative*dux_dz[INDEX_PERTURBATIVE];
  grad_matrix[1][0] = tau_perturbative*duy_dx[INDEX_PERTURBATIVE];
  grad_matrix[1][1] = tau_perturbative*duy_dy[INDEX_PERTURBATIVE];
  grad_matrix[1][2] = tau_perturbative*duy_dz[INDEX_PERTURBATIVE];
  grad_matrix[2][0] = tau_perturbative*duz_dx[INDEX_PERTURBATIVE];
  grad_matrix[2][1] = tau_perturbative*duz_dy[INDEX_PERTURBATIVE];
  grad_matrix[2][2] = tau_perturbative*duz_dz[INDEX_PERTURBATIVE];
  vect_grad[0] = grad_matrix[0][0];
  vect_grad[1] = grad_matrix[0][1];
  vect_grad[2] = grad_matrix[0][2];
  vect_grad[3] = grad_matrix[1][0];
  vect_grad[4] = grad_matrix[1][1];
  vect_grad[5] = grad_matrix[1][2];
  vect_grad[6] = grad_matrix[2][0];
  vect_grad[7] = grad_matrix[2][1];
  vect_grad[8] = grad_matrix[2][2];
	
  gsl_matrix_view m = gsl_matrix_view_array(vect_grad,3,3);
  gsl_matrix_view em = gsl_matrix_view_array(vect_forexpgrad,3,3);
  gsl_linalg_exponential_ss(&m.matrix,&em.matrix, .001); 
  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_grad[i][j] = gsl_matrix_get(&em.matrix, i,j);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_gradT[i][j]=Exp_grad[j][i];
		
  RTAU[0] = Exp_grad[0][0]*R0[INDEX_PERTURBATIVE][0]+Exp_grad[0][1]*R0[INDEX_PERTURBATIVE][1] + Exp_grad[0][2]*R0[INDEX_PERTURBATIVE][2];
  RTAU[1] = Exp_grad[1][0]*R0[INDEX_PERTURBATIVE][0]+Exp_grad[1][1]*R0[INDEX_PERTURBATIVE][1] + Exp_grad[1][2]*R0[INDEX_PERTURBATIVE][2];
  RTAU[2] = Exp_grad[2][0]*R0[INDEX_PERTURBATIVE][0]+Exp_grad[2][1]*R0[INDEX_PERTURBATIVE][1] + Exp_grad[2][2]*R0[INDEX_PERTURBATIVE][2];

  px[INDEX_PERTURBATIVE]=-Exp_gradT[0][0]*RTAU[0]-Exp_gradT[0][1]*RTAU[1] -Exp_gradT[0][2]*RTAU[2]; 
  py[INDEX_PERTURBATIVE]=-Exp_gradT[1][0]*RTAU[0]-Exp_gradT[1][1]*RTAU[1] -Exp_gradT[1][2]*RTAU[2]; 
  pz[INDEX_PERTURBATIVE]=-Exp_gradT[2][0]*RTAU[0]-Exp_gradT[2][1]*RTAU[1] -Exp_gradT[2][2]*RTAU[2]; 
	
  scra=sqrt(px[INDEX_PERTURBATIVE]*px[INDEX_PERTURBATIVE]+py[INDEX_PERTURBATIVE]*py[INDEX_PERTURBATIVE] + pz[INDEX_PERTURBATIVE]*pz[INDEX_PERTURBATIVE]);   
  px[INDEX_PERTURBATIVE]/=scra;  
  py[INDEX_PERTURBATIVE]/=scra;  
  pz[INDEX_PERTURBATIVE]/=scra;
		
} //end function control PERTURBATIVE for initialization


// Initialization FBSM with PERTURBATIVE OPTIMAL CONTROL
void rk4_forward_R_PO_initialization(double Rv_x[], double Rv_y[], double Rv_z[],double pv_x[], double pv_y[], double pv_z[], int ip){
  double k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y, k4z;
  int t;
  double scra, filter;
  double tensor_A[9];
	
  for (t=0;t<NITER+1;t++){
    control_perturbative_forinitialization(Rv_x,Rv_y,Rv_z,t);
    pv_x[t] = px[INDEX_PERTURBATIVE];
    pv_y[t] = py[INDEX_PERTURBATIVE];
    pv_z[t] = pz[INDEX_PERTURBATIVE];
    tensor_A[0] = AA_history[matrix_index(0,0,t)];
    tensor_A[1] = AA_history[matrix_index(0,1,t)];
    tensor_A[2] = AA_history[matrix_index(0,2,t)];
    tensor_A[3] = AA_history[matrix_index(1,0,t)];
    tensor_A[4] = AA_history[matrix_index(1,1,t)];
    tensor_A[5] = AA_history[matrix_index(1,2,t)];
    tensor_A[6] = AA_history[matrix_index(2,0,t)];
    tensor_A[7] = AA_history[matrix_index(2,1,t)];
    tensor_A[8] = AA_history[matrix_index(2,2,t)];
	  
    //step1//
    derivs_ip(Rv_x[t],Rv_y[t], Rv_z[t], tensor_A);// output ux e uy
    scra = sqrt(Rv_x[t]*Rv_x[t] + Rv_y[t]*Rv_y[t] + Rv_z[t]*Rv_z[t]);
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k1x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k1y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k1z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //step2//
    tensor_A[0] = (AA_history[matrix_index(0,0,t)] + AA_history[matrix_index(0,0,t+1)] )/ 2.;
    tensor_A[1] = (AA_history[matrix_index(0,1,t)] + AA_history[matrix_index(0,1,t+1)] )/ 2.;
    tensor_A[2] = (AA_history[matrix_index(0,2,t)] + AA_history[matrix_index(0,2,t+1)] )/ 2.;
    tensor_A[3] = (AA_history[matrix_index(1,0,t)] + AA_history[matrix_index(1,0,t+1)] )/ 2.;
    tensor_A[4] = (AA_history[matrix_index(1,1,t)] + AA_history[matrix_index(1,1,t+1)] )/ 2.;
    tensor_A[5] = (AA_history[matrix_index(1,2,t)] + AA_history[matrix_index(1,2,t+1)] )/ 2.;
    tensor_A[6] = (AA_history[matrix_index(2,0,t)] + AA_history[matrix_index(2,0,t+1)] )/ 2.;
    tensor_A[7] = (AA_history[matrix_index(2,1,t)] + AA_history[matrix_index(2,1,t+1)] )/ 2.;
    tensor_A[8] = (AA_history[matrix_index(2,2,t)] + AA_history[matrix_index(2,2,t+1)] )/ 2.;
	  
    derivs_ip(Rv_x[t]+k1x/2,Rv_y[t]+k1y/2,Rv_z[t]+k1z/2, tensor_A);
    scra = sqrt((Rv_x[t]+ k1x/2)*(Rv_x[t]+ k1x/2)+ (Rv_y[t]+ k1y/2)*(Rv_y[t]+ k1y/2) + (Rv_z[t]+ k1z/2)*(Rv_z[t]+ k1z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k2x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k2y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k2z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //step3//
    derivs_ip(Rv_x[t]+k2x/2,Rv_y[t]+k2y/2,Rv_z[t]+k2z/2, tensor_A);
    scra = sqrt((Rv_x[t]+ k2x/2)*(Rv_x[t]+ k2x/2)+ (Rv_y[t]+ k2y/2)*(Rv_y[t]+ k2y/2) +  (Rv_z[t]+ k2z/2)*(Rv_z[t]+ k2z/2));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k3x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k3y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k3z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //step4//
    tensor_A[0] = AA_history[matrix_index(0,0,t+1)];
    tensor_A[1] = AA_history[matrix_index(0,1,t+1)];
    tensor_A[2] = AA_history[matrix_index(0,2,t+1)];
    tensor_A[3] = AA_history[matrix_index(1,0,t+1)];
    tensor_A[4] = AA_history[matrix_index(1,1,t+1)];
    tensor_A[5] = AA_history[matrix_index(1,2,t+1)];
    tensor_A[6] = AA_history[matrix_index(2,0,t+1)];
    tensor_A[7] = AA_history[matrix_index(2,1,t+1)];
    tensor_A[8] = AA_history[matrix_index(2,2,t+1)];
    derivs_ip(Rv_x[t]+k3x,Rv_y[t]+k3y, Rv_z[t]+k3z, tensor_A);
    scra = sqrt((Rv_x[t]+ k3x)*(Rv_x[t]+ k3x)+ (Rv_y[t]+ k3y)*(Rv_y[t]+ k3y) +(Rv_z[t]+ k3z)*(Rv_z[t]+ k3z));
    filter = (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    k4x = dt*(filter*(uxx + vs[ip]*pv_x[t]));
    k4y = dt*(filter*(uyy + vs[ip]*pv_y[t]));
    k4z = dt*(filter*(uzz + vs[ip]*pv_z[t]));
	  
    //integration
    Rv_x[t+1] = Rv_x[t] + (k1x/6 + k2x/3 + k3x/3 + k4x/6);
    Rv_y[t+1] = Rv_y[t] + (k1y/6 + k2y/3 + k3y/3 + k4y/6);
    Rv_z[t+1] = Rv_z[t] + (k1z/6 + k2z/3 + k3z/3 + k4z/6);
	  
  }
}


//________________________________________________________________________________________//
/* Runge-Kutta backward integration for the lagrangian multiplier in the FBSM algorithm */

void rk4_backward_phi(double phi_x[], double phi_y[], double phi_z[],double Rv_x[], double Rv_y[], double Rv_z[], double pv_x[], double pv_y[], double pv_z[], int iep, int iteration,int ip){
  double k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y,k4z;
  int t,j;
  double scra, filter, deriv_filter, phiAR, phivp;
  double A00,A01,A02,A10, A11, A12, A20, A21, A22;
  
  for (t=0;t<NITER;t++){
    j = NITER-t; 	    
    A00 = AA_history[matrix_index(0,0,j)];
    A01 = AA_history[matrix_index(0,1,j)];
    A02 = AA_history[matrix_index(0,2,j)];
    A10 = AA_history[matrix_index(1,0,j)];
    A11 = AA_history[matrix_index(1,1,j)];
    A12 = AA_history[matrix_index(1,2,j)];
    A20 = AA_history[matrix_index(2,0,j)];
    A21 = AA_history[matrix_index(2,1,j)];
    A22 = AA_history[matrix_index(2,2,j)];
		
    //step1//
    scra = sqrt(Rv_x[j]*Rv_x[j] + Rv_y[j]*Rv_y[j] + Rv_z[j]*Rv_z[j]);
    filter =  (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    deriv_filter = alpha/(cosh(alpha * (scra - capture_distance)/capture_distance) * cosh(alpha * (scra - capture_distance)/capture_distance)  * capture_distance * scra * 2);
    phiAR = phi_x[j] * (A00* Rv_x[j] + A01*Rv_y[j] + A02*Rv_z[j]) + phi_y[j] * ( A10 * Rv_x[j] + A11*Rv_y[j] + A12*Rv_z[j]) + phi_z[j] * ( A20 * Rv_x[j] + A21*Rv_y[j] + A22*Rv_z[j]);
    phivp = vs[ip] * (phi_x[j] * pv_x[j] + phi_y[j] * pv_y[j] + phi_z[j] * pv_z[j]);		
    k1x = - (filter*(A00 * phi_x[j] + A10 * phi_y[j] + A20 * phi_z[j]) + deriv_filter * Rv_x[j] * (phiAR + phivp + capture_weight)  ) * dt;
    k1y = - (filter*(A01 * phi_x[j] + A11 * phi_y[j] + A21 * phi_z[j]) + deriv_filter * Rv_y[j] * (phiAR + phivp + capture_weight)  ) * dt;
    k1z = - (filter*(A02 * phi_x[j] + A12 * phi_y[j] + A22 * phi_z[j]) + deriv_filter * Rv_z[j] * (phiAR + phivp + capture_weight)  ) * dt;
		
    //step2//
    A00 = (AA_history[matrix_index(0,0,j)] + AA_history[matrix_index(0,0,j-1)] )/ 2 ;  
    A01 = (AA_history[matrix_index(0,1,j)] + AA_history[matrix_index(0,1,j-1)] )/ 2 ;  
    A02 = (AA_history[matrix_index(0,2,j)] + AA_history[matrix_index(0,2,j-1)] )/ 2 ;  
    A10 = (AA_history[matrix_index(1,0,j)] + AA_history[matrix_index(1,0,j-1)] )/ 2 ;  
    A11 = (AA_history[matrix_index(1,1,j)] + AA_history[matrix_index(1,1,j-1)] )/ 2 ;  
    A12 = (AA_history[matrix_index(1,2,j)] + AA_history[matrix_index(1,2,j-1)] )/ 2 ;  
    A20 = (AA_history[matrix_index(2,0,j)] + AA_history[matrix_index(2,0,j-1)] )/ 2 ;
    A21 = (AA_history[matrix_index(2,1,j)] + AA_history[matrix_index(2,1,j-1)] )/ 2 ;    
    A22 = (AA_history[matrix_index(2,2,j)] + AA_history[matrix_index(2,2,j-1)] )/ 2 ; 
    scra = sqrt(0.5*(Rv_x[j] + Rv_x[j-1])*0.5*(Rv_x[j] + Rv_x[j-1]) + 0.5*(Rv_y[j] + Rv_y[j-1])*0.5*(Rv_y[j] + Rv_y[j-1]) + 0.5*(Rv_z[j] + Rv_z[j-1])*0.5*(Rv_z[j] + Rv_z[j-1]));
    filter =  (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    deriv_filter = alpha/(cosh(alpha * (scra - capture_distance)/capture_distance) * cosh(alpha * (scra - capture_distance)/capture_distance)  * capture_distance * scra * 2);
    phiAR = (phi_x[j] - k1x/2)* (A00 *  0.5*(Rv_x[j] + Rv_x[j-1])+ A01*0.5*(Rv_y[j] + Rv_y[j-1]) +  A02*0.5*(Rv_z[j] + Rv_z[j-1])) + (phi_y[j] - k1y/2) *( A10 * 0.5*(Rv_x[j] + Rv_x[j-1]) + A11*0.5*(Rv_y[j] + Rv_y[j-1]) + A12*0.5*(Rv_z[j] + Rv_z[j-1])) + (phi_z[j] -k1z/2) * (A20* 0.5*(Rv_x[j] + Rv_x[j-1])/2 + A21 * (Rv_y[j]+Rv_y[j-1]) + A22*(Rv_z[j] + Rv_z[j-1]));
    phivp = vs[ip] * ((phi_x[j] - k1x/2) * 0.5*(pv_x[j]+ pv_x[j-1]) + (phi_y[j] - k1y/2)* 0.5*(pv_y[j] + pv_y[j-1]) + (phi_z[j] - k1z/2)* 0.5*(pv_z[j] + pv_z[j-1]));
    k2x = - (filter*(A00 * (phi_x[j] - k1x/2) + A10 * (phi_y[j] - k1y/2) + A20 * (phi_z[j] - k1z/2))  + deriv_filter * 0.5*(Rv_x[j] +Rv_x[j-1]) * (phiAR  +phivp +  capture_weight))*dt;
    k2y = - (filter*(A01 * (phi_x[j] - k1x/2) + A11 * (phi_y[j] - k1y/2) + A21 * (phi_z[j] - k1z/2))  + deriv_filter * 0.5*(Rv_y[j] +Rv_y[j-1]) * (phiAR  + phivp + capture_weight))*dt;
    k2z = - (filter*(A02 * (phi_x[j] - k1x/2) + A12 * (phi_y[j] - k1y/2) + A22 * (phi_z[j] - k1z/2))  + deriv_filter * 0.5*(Rv_z[j] +Rv_z[j-1]) * (phiAR  + phivp + capture_weight))*dt;
		
    //step3//
    scra = sqrt(0.5*(Rv_x[j] + Rv_x[j-1])*0.5*(Rv_x[j] + Rv_x[j-1]) + 0.5*(Rv_y[j] + Rv_y[j-1])*0.5*(Rv_y[j] + Rv_y[j-1]) + 0.5*(Rv_z[j] + Rv_z[j-1])*0.5*(Rv_z[j] + Rv_z[j-1]));
    filter =  (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    deriv_filter = alpha/(cosh(alpha * (scra - capture_distance)/capture_distance) * cosh(alpha * (scra - capture_distance)/capture_distance)  * capture_distance * scra * 2);
    phiAR = (phi_x[j] - k2x/2)* (A00 *  0.5*(Rv_x[j] + Rv_x[j-1])+ A01*0.5*(Rv_y[j] + Rv_y[j-1]) +  A02*0.5*(Rv_z[j] + Rv_z[j-1])) + (phi_y[j] - k2y/2) *( A10 * 0.5*(Rv_x[j] + Rv_x[j-1]) + A11*0.5*(Rv_y[j] + Rv_y[j-1]) + A12*0.5*(Rv_z[j] + Rv_z[j-1])) + (phi_z[j] -k2z/2) * (A20* 0.5*(Rv_x[j] + Rv_x[j-1])/2 + A21 * (Rv_y[j]+Rv_y[j-1]) + A22*(Rv_z[j] + Rv_z[j-1]));
    phivp = vs[ip] * ((phi_x[j] - k2x/2) * 0.5*(pv_x[j]+ pv_x[j-1]) + (phi_y[j] - k2y/2)* 0.5*(pv_y[j] + pv_y[j-1]) + (phi_z[j] - k2z/2)* 0.5*(pv_z[j] + pv_z[j-1]));
    k3x = - (filter*(A00 * (phi_x[j] - k2x/2) + A10 * (phi_y[j] - k2y/2) + A20 * (phi_z[j] - k2z/2))  + deriv_filter * 0.5*(Rv_x[j] +Rv_x[j-1]) * (phiAR  +phivp +  capture_weight))*dt;
    k3y = - (filter*(A01 * (phi_x[j] - k2x/2) + A11 * (phi_y[j] - k2y/2) + A21 * (phi_z[j] - k2z/2))  + deriv_filter * 0.5*(Rv_y[j] +Rv_y[j-1]) * (phiAR  + phivp + capture_weight))*dt;
    k3z = - (filter*(A02 * (phi_x[j] - k2x/2) + A12 * (phi_y[j] - k2y/2) + A22 * (phi_z[j] - k2z/2))  + deriv_filter * 0.5*(Rv_z[j] +Rv_z[j-1]) * (phiAR  + phivp + capture_weight))*dt;
		
		
    //step4//
    A00 = AA_history[matrix_index(0,0,j-1)];
    A01 = AA_history[matrix_index(0,1,j-1)];
    A02 = AA_history[matrix_index(0,2,j-1)];
    A10 = AA_history[matrix_index(1,0,j-1)];
    A11 = AA_history[matrix_index(1,1,j-1)];
    A12 = AA_history[matrix_index(1,2,j-1)];
    A20 = AA_history[matrix_index(2,0,j-1)];
    A21 = AA_history[matrix_index(2,1,j-1)];
    A22 = AA_history[matrix_index(2,2,j-1)];
    scra = sqrt(Rv_x[j-1]*Rv_x[j-1] + Rv_y[j-1]*Rv_y[j-1] + Rv_z[j-1]*Rv_z[j-1]);
    filter =  (1 + tanh(alpha * (scra-capture_distance)/capture_distance ))/2;
    deriv_filter = alpha/(cosh(alpha * (scra - capture_distance)/capture_distance) * cosh(alpha * (scra - capture_distance)/capture_distance)  * capture_distance * scra * 2);
    phiAR = (phi_x[j] -k3x) * (A00 * Rv_x[j-1] + A01*Rv_y[j-1] + A02*Rv_x[j-1]) + (phi_y[j]-k3x) * ( A10 * Rv_x[j-1] + A11*Rv_y[j-1] + A12*Rv_z[j-1]) + (phi_z[j]-k3z) * (A20 * Rv_x[j-1] + A21 * Rv_y[j-1] + A22 * Rv_z[j-1]);
    phivp = vs[ip] * ((phi_x[j]-k3x) * pv_x[j-1] + (phi_y[j]-k3y) * pv_y[j-1] + (phi_z[j]-k3z) * pv_z[j-1]);
    k4x = - (filter*(A00 * (phi_x[j]-k3x) + A10* (phi_y[j]-k3y) + A20*(phi_z[j]-k3z)) + deriv_filter * Rv_x[j-1] * (phiAR  + phivp + capture_weight)  )*dt;
    k4y = - (filter*(A01 * (phi_x[j]-k3x) + A11 * (phi_y[j]-k3y) +A21*(phi_z[j]-k3z)  )  + deriv_filter * Rv_y[j-1] * (phiAR + phivp + capture_weight)  )  *dt;
    k4z = - (filter*(A02 * (phi_x[j]-k3x) + A12 * (phi_y[j]-k3y) +A22*(phi_z[j]-k3z)  )  + deriv_filter * Rv_z[j-1] * (phiAR + phivp + capture_weight)  )  *dt;
		
    //backward integration
    phi_x[j-1] = phi_x[j] -  (k1x/6 + k2x/3 + k3x/3 + k4x/6);
    phi_y[j-1] = phi_y[j] -  (k1y/6 + k2y/3 + k3y/3 + k4y/6);
    phi_z[j-1] = phi_z[j] -  (k1z/6 + k2z/3 + k3z/3 + k4z/6);
  }
}


/* ********************************************* CONTROL FUNCTIONS ********************************************* */

/*---------------------------------------------------------------------------------------*/
/* OPTIMAL CONTROL - FORWARD-BACKWARD SWEEM METHOD (FBSM) ALGORITHM */
/*---------------------------------------------------------------------------------------*/
void control_ForwardBackwardSweepMethod(double*x, double*y, double *z, int iep){
  int i,t, iteration, k, max_step_checkconv;
  double test,captPP,captSC,captPO,capturetime_iteration;
  double temp[9];
  double scra;
  double min;
  double Rv_x[NITER+1], Rv_y[NITER+1], Rv_z[NITER+1];
  double oldRv_x[NITER+1], oldRv_y[NITER+1], oldRv_z[NITER+1];
  double pv_x[NITER+1], pv_y[NITER+1], pv_z[NITER+1];
  double pv1_x[NITER+1], pv1_y[NITER+1], pv1_z[NITER+1];
  double oldpv_x[NITER+1], oldpv_y[NITER+1], oldpv_z[NITER+1];
  double phi_x[NITER+1], phi_y[NITER+1], phi_z[NITER+1];
  double oldphi_x[NITER+1], oldphi_y[NITER+1], oldphi_z[NITER+1];  
	       	
  iteration = 0; // algorithm iteration
  test = -1; // check te convergence; initialization of the test variable to be < 0 
		
  //***** step 1. initialization *****//
  for (t=0;t<NITER+1;t++){
    pv_x[t] = 0;
    pv_y[t] = 0;
    pv_z[t] = 0;
    Rv_x[t] = Rv_y[t] = Rv_z[t] = 0; 
    phi_x[t] = phi_y[t] = phi_z[t] =  0;
  }

  // Check the best heuristic strategy to initialize the algorithm (1st possibility Pure Pursuit)
  Rv_x[0] = (x[INDEX_OPTCTRL] - x[0]); // Fix the initial condition for the state variable
  Rv_y[0] = (y[INDEX_OPTCTRL] - y[0]);
  Rv_z[0] = (z[INDEX_OPTCTRL] - z[0]);  
  rk4_forward_R_PP_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);   
  capturetime_iteration=0;  
  for(t=0;t<NITER+1;t++){
    scra=sqrt(Rv_x[t]*Rv_x[t]+Rv_y[t]*Rv_y[t]+Rv_z[t]*Rv_z[t]);
    if (capturetime_iteration==0 && scra<=capture_distance)
      capturetime_iteration=t*dt;
  }
  if(capturetime_iteration==0)
    capturetime_iteration=NITER*dt;
  captPP=capturetime_iteration;
  scra=sqrt(Rv_x[NITER]*Rv_x[NITER]+Rv_y[NITER]*Rv_y[NITER]+Rv_z[NITER]*Rv_z[NITER]);
  J_PP=scra*scra + capture_weight*capturetime_iteration;
  // Check the best heuristic strategy to initialize the algorithm (2nd possibility Surfing Control)
  Rv_x[0] = (x[INDEX_OPTCTRL] - x[0]); // Fix the initial condition for the state variable
  Rv_y[0] = (y[INDEX_OPTCTRL] - y[0]);
  Rv_z[0] = (z[INDEX_OPTCTRL] - z[0]);	
  rk4_forward_R_SC_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);	
  capturetime_iteration=0;
  for(t=0;t<NITER+1;t++){
    scra=sqrt(Rv_x[t]*Rv_x[t]+Rv_y[t]*Rv_y[t]+Rv_z[t]*Rv_z[t]);
    if (capturetime_iteration==0 && scra<=capture_distance)
      capturetime_iteration=t*dt;
  }
  if(capturetime_iteration==0)
    capturetime_iteration=NITER*dt;
  captSC=capturetime_iteration;
  scra=sqrt(Rv_x[NITER]*Rv_x[NITER]+Rv_y[NITER]*Rv_y[NITER]+Rv_z[NITER]*Rv_z[NITER]);
  J_SC=scra*scra + capture_weight*capturetime_iteration;
  // Check the best heuristic strategy to initialize the algorithm (3rd possibility Perturbative Optimal control)
  Rv_x[0] = (x[INDEX_OPTCTRL] - x[0]); // Fix the initial condition for the state variable
  Rv_y[0] = (y[INDEX_OPTCTRL] - y[0]);
  Rv_z[0] = (z[INDEX_OPTCTRL] - z[0]);	
  rk4_forward_R_PO_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL); 
  capturetime_iteration=0;
  for(t=0;t<NITER+1;t++){
    scra=sqrt(Rv_x[t]*Rv_x[t]+Rv_y[t]*Rv_y[t]+Rv_z[t]*Rv_z[t]);
    if (capturetime_iteration==0 && scra<=capture_distance)
      capturetime_iteration=t*dt;
  }
  if(capturetime_iteration==0)
    capturetime_iteration=NITER*dt;
  captPO=capturetime_iteration;
  scra=sqrt(Rv_x[NITER]*Rv_x[NITER]+Rv_y[NITER]*Rv_y[NITER]+Rv_z[NITER]*Rv_z[NITER]);	
  J_PO=scra*scra + capture_weight*capturetime_iteration;

  // Check the best heuristic strategy to initialize the algorithm
  if(J_PP<J_SC && J_PP<J_PO){
    minJ_OC=J_PP;	  
    rk4_forward_R_PP_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);
  }	
  else if(J_SC<J_PP && J_SC<J_PO){
    minJ_OC=J_SC;	     
    rk4_forward_R_SC_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);
  }
  else if(J_PO<J_SC && J_PO<J_PP){
    minJ_OC=J_PO;
    rk4_forward_R_PO_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);
  }


  //***** ITERATION OVER THE FBSM ****//
  while (test <0){ 
    iteration += 1;
    for (t=0;t<NITER+1;t++){
      oldpv_x[t] = pv_x[t];
      oldpv_y[t] = pv_y[t];
      oldpv_z[t] = pv_z[t];
      oldRv_x[t] = Rv_x[t];
      oldRv_y[t] = Rv_y[t];
      oldRv_z[t] = Rv_z[t];
      oldphi_x[t] = phi_x[t];
      oldphi_y[t] = phi_y[t];
      oldphi_z[t] = phi_z[t];
    }
	  
    //***** step 2. forward evolution of the state variable R *****//
    rk4_forward_R(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL); 
	  
    //***** step 3. backward evolution of the co-state variable phi *****//
    phi_x[NITER] = Rv_x[NITER];  
    phi_y[NITER] = Rv_y[NITER];
    phi_z[NITER] = Rv_z[NITER];
    rk4_backward_phi(phi_x, phi_y, phi_z, Rv_x, Rv_y, Rv_z,  pv_x, pv_y, pv_z,  iep, iteration,INDEX_OPTCTRL); 
	  
    //***** step 4. update of the control variables *****//
    for (t=0;t<NITER+1;t++){	    
      pv1_x[t] = -phi_x[t];
      pv1_y[t] = -phi_y[t];
      pv1_z[t] = -phi_z[t];	    
      scra = sqrt(pv1_x[t] * pv1_x[t] + pv1_y[t] * pv1_y[t] + pv1_z[t] * pv1_z[t]);
      if(scra>0){ 
	pv1_x[t]/=scra;
	pv1_y[t]/=scra;
	pv1_z[t]/=scra;
      }
      if(scra==0){
	pv1_x[t] = 0;
	pv1_y[t] = 0;
	pv1_z[t] = 0;
      }	    	    
      /* weighted average with a learning rate gamma_lr */
      pv_x[t] = gamma_lr*pv1_x[t] + (1-gamma_lr)*oldpv_x[t]; 
      pv_y[t] = gamma_lr*pv1_y[t] + (1-gamma_lr)*oldpv_y[t]; 
      pv_z[t] = gamma_lr*pv1_z[t] + (1-gamma_lr)*oldpv_z[t]; 
      scra = sqrt(pv_x[t] * pv_x[t] + pv_y[t] * pv_y[t] + pv_z[t] * pv_z[t]); 
      if(scra>0){ 
	pv_x[t]/=scra;
	pv_y[t]/=scra;
	pv_z[t]/=scra;
      }
      if(scra==0){
	pv_x[t] = 0;
	pv_y[t]=  0;
	pv_z[t]=  0;
      } 	    
    } //end step 4
	  
    //***** step 5. check convergence *****//
    for (i=0; i<9; i++){
      temp[i] = 0;
    }	  
    for (t=0;t<NITER+1;t++){     
      temp[0] +=delta*fabs(pv_x[t]) - fabs(oldpv_x[t] - pv_x[t]); 
      temp[1] +=delta*fabs(pv_y[t]) - fabs(oldpv_y[t] - pv_y[t]); 
      temp[2] +=delta*fabs(pv_z[t]) - fabs(oldpv_z[t] - pv_z[t]); 
      temp[3] +=delta*fabs(Rv_x[t]) - fabs(oldRv_x[t] - Rv_x[t]);
      temp[4] +=delta*fabs(Rv_y[t]) - fabs(oldRv_y[t] - Rv_y[t]);
      temp[5] +=delta*fabs(Rv_z[t]) - fabs(oldRv_z[t] - Rv_z[t]);  
      temp[6] +=delta*fabs(phi_x[t]) - fabs(oldphi_x[t] - phi_x[t]);
      temp[7] +=delta*fabs(phi_y[t]) - fabs(oldphi_y[t] - phi_y[t]); 
      temp[8] +=delta*fabs(phi_z[t]) - fabs(oldphi_z[t] - phi_z[t]); 	    
    }
	  
    min = temp[0];
    k = 0;
    for (i=0; i<9; i++){
      if (min>temp[i]){
	min=temp[i];
	k = i;
      }
    }
    test = min;	  
    convergence = OK_CONVERGENCE;

    // Save the performance index along the iteration
    capturetime_iteration=0;    
    for(t=0;t<NITER+1;t++){
      scra=sqrt(Rv_x[t]*Rv_x[t]+Rv_y[t]*Rv_y[t]+Rv_z[t]*Rv_z[t]);
      if (capturetime_iteration==0 && scra<=capture_distance)
	capturetime_iteration=t*dt;
    }
    if(capturetime_iteration==0)
      capturetime_iteration=NITER*dt;
    scra=sqrt(Rv_x[NITER]*Rv_x[NITER]+Rv_y[NITER]*Rv_y[NITER]+Rv_z[NITER]*Rv_z[NITER]);
    fprintf(performance_index, "%d %d %g %g %g \n", iep, iteration, scra*scra, capturetime_iteration, scra*scra + capture_weight*capturetime_iteration);

    J_OC=scra*scra + capture_weight*capturetime_iteration;
    if (J_OC<minJ_OC){
      minJ_OC=J_OC;     
    }
	  
    if(iteration>=max_iteration){
      test=1;
      convergence=NO_CONVERGENCE;
    }
	  
  } // end while  
	
	
  //***** step 6. assign the control *****//
  for (t=0;t<NITER+1;t++){
    px_ott[t] = pv_x[t];
    py_ott[t] = pv_y[t];
    pz_ott[t] = pv_z[t];
    rk4_forward_R(Rv_x, Rv_y, Rv_z,pv_x, pv_y, pv_z, INDEX_OPTCTRL);
    Rx_ott[t] = Rv_x[t];
    Ry_ott[t] = Rv_y[t];
    Rz_ott[t] = Rv_z[t];
	  
    if(sqrt(Rx_ott[t]*Rx_ott[t] + Ry_ott[t]*Ry_ott[t]+ Rz_ott[t]*Rz_ott[t]) <= capture_distance){
      if(capturetime_OC_FBSM == 0)
	capturetime_OC_FBSM = (double) t*dt;	    
    }
	  
  }
	
  // output variable
  capture = NO_CAPTURE;       
  if (sqrt(Rv_x[NITER]*Rv_x[NITER] + Rv_y[NITER]*Rv_y[NITER] + Rv_z[NITER]*Rv_z[NITER]) <= capture_distance)
    capture = OK_CAPTURE;	
  fprintf(performance_index, "\n\n");
  
} //end function control_ForwardBackwardSweepMethod


void control_OC_maxiteration(double*x, double*y, double *z, int iep, int max_NumIteration){
  int i,t, iteration, k, max_step_checkconv; 
  double temp[9];
  double scra;
  double min;
  double capturetime_iteration;
  double Rv_x[NITER+1], Rv_y[NITER+1], Rv_z[NITER+1];
  double oldRv_x[NITER+1], oldRv_y[NITER+1], oldRv_z[NITER+1];
  double pv_x[NITER+1], pv_y[NITER+1], pv_z[NITER+1];
  double pv1_x[NITER+1], pv1_y[NITER+1], pv1_z[NITER+1];
  double oldpv_x[NITER+1], oldpv_y[NITER+1], oldpv_z[NITER+1];
  double phi_x[NITER+1], phi_y[NITER+1], phi_z[NITER+1];
  double oldphi_x[NITER+1], oldphi_y[NITER+1], oldphi_z[NITER+1];
   
  //***** step 1. initialization *****//
  for (t=0;t<NITER+1;t++){
    pv_x[t] = 0;
    pv_y[t] = 0;
    pv_z[t] = 0;
    Rv_x[t] = Rv_y[t] = Rv_z[t] = 0; 
    phi_x[t] = phi_y[t] = phi_z[t] =  0;
  }	
  Rv_x[0] = (x[INDEX_OPTCTRL] - x[0]); // Fix the initial condition for the state variable R
  Rv_y[0] = (y[INDEX_OPTCTRL] - y[0]);
  Rv_z[0] = (z[INDEX_OPTCTRL] - z[0]);	
	
  // ----------------------- INITIALIZATION ------------------------ //
  rk4_forward_R_PP_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL); 	
  capturetime_iteration=0;
  for(t=0;t<NITER+1;t++){
    scra=sqrt(Rv_x[t]*Rv_x[t]+Rv_y[t]*Rv_y[t]+Rv_z[t]*Rv_z[t]);
    if (capturetime_iteration==0 && scra<=capture_distance)
      capturetime_iteration=t*dt;
  }
  if(capturetime_iteration==0)
    capturetime_iteration=NITER*dt;
  scra=sqrt(Rv_x[NITER]*Rv_x[NITER]+Rv_y[NITER]*Rv_y[NITER]+Rv_z[NITER]*Rv_z[NITER]);        
  J_PP=scra*scra + capture_weight*capturetime_iteration;

  rk4_forward_R_SC_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL); 
  capturetime_iteration=0;
  for(t=0;t<NITER+1;t++){
    scra=sqrt(Rv_x[t]*Rv_x[t]+Rv_y[t]*Rv_y[t]+Rv_z[t]*Rv_z[t]);
    if (capturetime_iteration==0 && scra<=capture_distance)
      capturetime_iteration=t*dt;
  }
  if(capturetime_iteration==0)
    capturetime_iteration=NITER*dt;
  scra=sqrt(Rv_x[NITER]*Rv_x[NITER]+Rv_y[NITER]*Rv_y[NITER]+Rv_z[NITER]*Rv_z[NITER]);
  J_SC=scra*scra + capture_weight*capturetime_iteration;

  rk4_forward_R_PO_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL); 
  capturetime_iteration=0;
  for(t=0;t<NITER+1;t++){
    scra=sqrt(Rv_x[t]*Rv_x[t]+Rv_y[t]*Rv_y[t]+Rv_z[t]*Rv_z[t]);
    if (capturetime_iteration==0 && scra<=capture_distance)
      capturetime_iteration=t*dt;
  }
  if(capturetime_iteration==0)
    capturetime_iteration=NITER*dt;
  scra=sqrt(Rv_x[NITER]*Rv_x[NITER]+Rv_y[NITER]*Rv_y[NITER]+Rv_z[NITER]*Rv_z[NITER]);
  J_PO=scra*scra + capture_weight*capturetime_iteration;

  if(J_PP<=J_SC && J_PP<=J_PO){
    minJ_OC=J_PP;
    rk4_forward_R_PP_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);
  }
  else if(J_SC<=J_PP && J_SC<=J_PO){
    minJ_OC=J_SC;
    rk4_forward_R_SC_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);
  }
  else if(J_PO<=J_SC && J_PO<=J_PP){
    minJ_OC=J_PO;
    rk4_forward_R_PO_initialization(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);
  }
        
  //***** ITERATION OVER THE FBSM  ****//
  for (iteration=1;iteration<max_NumIteration;iteration++){	  
    for (t=0;t<NITER+1;t++){
      oldpv_x[t] = pv_x[t];
      oldpv_y[t] = pv_y[t];
      oldpv_z[t] = pv_z[t];
      oldRv_x[t] = Rv_x[t];
      oldRv_y[t] = Rv_y[t];
      oldRv_z[t] = Rv_z[t];
      oldphi_x[t] = phi_x[t];
      oldphi_y[t] = phi_y[t];
      oldphi_z[t] = phi_z[t];
    }
    
    //***** step 2. forward evolution of the state variable R *****//
    rk4_forward_R(Rv_x, Rv_y, Rv_z, pv_x, pv_y, pv_z, INDEX_OPTCTRL);
    
    //***** step 3. backward evolution of the co-state phi *****//
    phi_x[NITER] = Rv_x[NITER];  
    phi_y[NITER] = Rv_y[NITER];
    phi_z[NITER] = Rv_z[NITER];
    rk4_backward_phi(phi_x, phi_y, phi_z, Rv_x, Rv_y, Rv_z,  pv_x, pv_y, pv_z,  iep, iteration,INDEX_OPTCTRL); 
	  
    //***** step 4. update the control variables *****//
    for (t=0;t<NITER+1;t++){	    
      pv1_x[t] = -phi_x[t];
      pv1_y[t] = -phi_y[t];
      pv1_z[t] = -phi_z[t];      
      scra = sqrt(pv1_x[t] * pv1_x[t] + pv1_y[t] * pv1_y[t] + pv1_z[t] * pv1_z[t]);
      if(scra>0){ 
	pv1_x[t]/=scra;
	pv1_y[t]/=scra;
	pv1_z[t]/=scra;
      }
      if(scra==0){
	pv1_x[t] = 0;
	pv1_y[t] = 0;
	pv1_z[t] = 0;
      }
	    	         
      pv_x[t] = gamma_lr*pv1_x[t] + (1-gamma_lr)*oldpv_x[t]; 
      pv_y[t] = gamma_lr*pv1_y[t] + (1-gamma_lr)*oldpv_y[t]; 
      pv_z[t] = gamma_lr*pv1_z[t] + (1-gamma_lr)*oldpv_z[t]; 
      scra = sqrt(pv_x[t] * pv_x[t] + pv_y[t] * pv_y[t] + pv_z[t] * pv_z[t]); 
      if(scra>0){ 
	pv_x[t]/=scra;
	pv_y[t]/=scra;
	pv_z[t]/=scra;
      }
      if(scra==0){
	pv_x[t] = 0;
	pv_y[t]=  0;
	pv_z[t]=  0;
      } 
	    
    } //end step 4
	  	  
    if(iteration>=max_NumIteration){
      fprintf(stderr, "error");
      exit(0);	   
    }
	  
  } // end cycle iteration  
	
	
  //***** step 6. assign the control *****//
  for (t=0;t<NITER+1;t++){
    px_ott[t] = pv_x[t];
    py_ott[t] = pv_y[t];
    pz_ott[t] = pv_z[t];
    rk4_forward_R(Rv_x, Rv_y, Rv_z,pv_x, pv_y, pv_z, INDEX_OPTCTRL);
    Rx_ott[t] = Rv_x[t];
    Ry_ott[t] = Rv_y[t];
    Rz_ott[t] = Rv_z[t];
    
    if(sqrt(Rx_ott[t]*Rx_ott[t] + Ry_ott[t]*Ry_ott[t]+ Rz_ott[t]*Rz_ott[t]) <= capture_distance){
      if(capturetime_OC_FBSM == 0)
	capturetime_OC_FBSM = (double) t*dt;	    
    }
    
  }
  capture = NO_CAPTURE;
  
  if (sqrt(Rv_x[NITER]*Rv_x[NITER] + Rv_y[NITER]*Rv_y[NITER] + Rv_z[NITER]*Rv_z[NITER]) <= capture_distance)
      capture = OK_CAPTURE;
    
} //end function FBSM max iteration fixed


/*----------------------------------------------------------------------------------------*/
/* SURFING CONTROL */
/*----------------------------------------------------------------------------------------*/

void control_surfing(double *x,double *y, double *z, int t){
  int i,j;
  double dux_dx[NPART],dux_dy[NPART], dux_dz[NPART],duy_dx[NPART],duy_dy[NPART], duy_dz[NPART], duz_dx[NPART], duz_dy[NPART], duz_dz[NPART]; 
  double R0[NPART][3]; 
  double scra;
  double vect_grad[9]; //vector that containts all the elements of the gradients matrix, i.e., grad_matrix[i,j] = vect_grad[i*3+j] 
  double vect_forexpgrad[9] = {0};
	
  gradients(x,y,z,dux_dx,dux_dy,dux_dz,duy_dx,duy_dy, duy_dz, duz_dx, duz_dy, duz_dz, t); 

  scra = dist(x,y,z,INDEX_SURFING);
  R0[INDEX_SURFING][0] = (x[INDEX_SURFING] - x[0])/scra; 
  R0[INDEX_SURFING][1] = (y[INDEX_SURFING] - y[0])/scra;
  R0[INDEX_SURFING][2] = (z[INDEX_SURFING] - z[0])/scra;	 

  grad_matrix[0][0] = tau_surfing*dux_dx[INDEX_SURFING];
  grad_matrix[0][1] = tau_surfing*dux_dy[INDEX_SURFING];
  grad_matrix[0][2] = tau_surfing*dux_dz[INDEX_SURFING];
  grad_matrix[1][0] = tau_surfing*duy_dx[INDEX_SURFING];
  grad_matrix[1][1] = tau_surfing*duy_dy[INDEX_SURFING];
  grad_matrix[1][2] = tau_surfing*duy_dz[INDEX_SURFING];
  grad_matrix[2][0] = tau_surfing*duz_dx[INDEX_SURFING];
  grad_matrix[2][1] = tau_surfing*duz_dy[INDEX_SURFING];
  grad_matrix[2][2] = tau_surfing*duz_dz[INDEX_SURFING];
  vect_grad[0] = grad_matrix[0][0];
  vect_grad[1] = grad_matrix[0][1];
  vect_grad[2] = grad_matrix[0][2];
  vect_grad[3] = grad_matrix[1][0];
  vect_grad[4] = grad_matrix[1][1];
  vect_grad[5] = grad_matrix[1][2];
  vect_grad[6] = grad_matrix[2][0];
  vect_grad[7] = grad_matrix[2][1];
  vect_grad[8] = grad_matrix[2][2];
	
  gsl_matrix_view m = gsl_matrix_view_array(vect_grad,3,3);
  gsl_matrix_view em = gsl_matrix_view_array(vect_forexpgrad,3,3);
  gsl_linalg_exponential_ss(&m.matrix,&em.matrix, .001);

  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_grad[i][j] = gsl_matrix_get(&em.matrix, i,j);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_gradT[i][j] = Exp_grad[j][i];

  px[INDEX_SURFING] = -Exp_gradT[0][0]*R0[INDEX_SURFING][0] - Exp_gradT[0][1]*R0[INDEX_SURFING][1] - Exp_gradT[0][2]*R0[INDEX_SURFING][2];
  py[INDEX_SURFING] = -Exp_gradT[1][0]*R0[INDEX_SURFING][0] - Exp_gradT[1][1]*R0[INDEX_SURFING][1] - Exp_gradT[1][2]*R0[INDEX_SURFING][2]; 
  pz[INDEX_SURFING] = -Exp_gradT[2][0]*R0[INDEX_SURFING][0] - Exp_gradT[2][1]*R0[INDEX_SURFING][1] - Exp_gradT[2][2]*R0[INDEX_SURFING][2]; 
	
  scra=sqrt(px[INDEX_SURFING]*px[INDEX_SURFING]+py[INDEX_SURFING]*py[INDEX_SURFING] + pz[INDEX_SURFING]*pz[INDEX_SURFING]);   
  px[INDEX_SURFING]/=scra;  
  py[INDEX_SURFING]/=scra;  
  pz[INDEX_SURFING]/=scra;
					
} //end function control surfing


/*----------------------------------------------------------------------------------------*/
/* PURE PURSUIT CONTROL */
/*----------------------------------------------------------------------------------------*/
void control_purepursuit(double *x,double *y, double *z){
  double R0[3]; 
  double scra = dist(x,y,z,INDEX_PUREPURSUIT);

  R0[0] = (x[INDEX_PUREPURSUIT] - x[0])/scra; 
  R0[1] = (y[INDEX_PUREPURSUIT] - y[0])/scra;	
  R0[2] = (z[INDEX_PUREPURSUIT] - z[0])/scra;	 

  px[INDEX_PUREPURSUIT] = -R0[0];
  py[INDEX_PUREPURSUIT] = -R0[1];
  pz[INDEX_PUREPURSUIT] = -R0[2];
}


/*-----------------------------------------------------------------------------------------*/
/* PERTURBATIVE OPTIMAL CONTROL */
/*-----------------------------------------------------------------------------------------*/

void control_perturbative(double *x,double *y, double *z, int t){
  int i,j;
  double dux_dx[NPART],dux_dy[NPART], dux_dz[NPART],duy_dx[NPART],duy_dy[NPART], duy_dz[NPART], duz_dx[NPART], duz_dy[NPART], duz_dz[NPART];
  double R0[NPART][3];
  double RTAU[3];
  double scra;
  double vect_grad[9]; //vector that containts all the elements of the gradients matrix, i.e., grad_matrix[i,j] = vect_grad[i*3+j] 
  double vect_forexpgrad[9] = {0};	
	
  gradients(x,y,z,dux_dx,dux_dy,dux_dz,duy_dx,duy_dy, duy_dz, duz_dx, duz_dy, duz_dz, t); 

  scra = dist(x,y,z,INDEX_PERTURBATIVE);
  R0[INDEX_PERTURBATIVE][0] = (x[INDEX_PERTURBATIVE] - x[0])/scra; 
  R0[INDEX_PERTURBATIVE][1] = (y[INDEX_PERTURBATIVE] - y[0])/scra;
  R0[INDEX_PERTURBATIVE][2] = (z[INDEX_PERTURBATIVE] - z[0])/scra; 

  grad_matrix[0][0] = tau_perturbative*dux_dx[INDEX_PERTURBATIVE];
  grad_matrix[0][1] = tau_perturbative*dux_dy[INDEX_PERTURBATIVE];
  grad_matrix[0][2] = tau_perturbative*dux_dz[INDEX_PERTURBATIVE];
  grad_matrix[1][0] = tau_perturbative*duy_dx[INDEX_PERTURBATIVE];
  grad_matrix[1][1] = tau_perturbative*duy_dy[INDEX_PERTURBATIVE];
  grad_matrix[1][2] = tau_perturbative*duy_dz[INDEX_PERTURBATIVE];
  grad_matrix[2][0] = tau_perturbative*duz_dx[INDEX_PERTURBATIVE];
  grad_matrix[2][1] = tau_perturbative*duz_dy[INDEX_PERTURBATIVE];
  grad_matrix[2][2] = tau_perturbative*duz_dz[INDEX_PERTURBATIVE];
  vect_grad[0] = grad_matrix[0][0];
  vect_grad[1] = grad_matrix[0][1];
  vect_grad[2] = grad_matrix[0][2];
  vect_grad[3] = grad_matrix[1][0];
  vect_grad[4] = grad_matrix[1][1];
  vect_grad[5] = grad_matrix[1][2];
  vect_grad[6] = grad_matrix[2][0];
  vect_grad[7] = grad_matrix[2][1];
  vect_grad[8] = grad_matrix[2][2];
	
  gsl_matrix_view m = gsl_matrix_view_array(vect_grad,3,3);
  gsl_matrix_view em = gsl_matrix_view_array(vect_forexpgrad,3,3);
  gsl_linalg_exponential_ss(&m.matrix,&em.matrix, .001); 
  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_grad[i][j] = gsl_matrix_get(&em.matrix, i,j);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      Exp_gradT[i][j]=Exp_grad[j][i];
		
  RTAU[0] = Exp_grad[0][0]*R0[INDEX_PERTURBATIVE][0]+Exp_grad[0][1]*R0[INDEX_PERTURBATIVE][1] + Exp_grad[0][2]*R0[INDEX_PERTURBATIVE][2];
  RTAU[1] = Exp_grad[1][0]*R0[INDEX_PERTURBATIVE][0]+Exp_grad[1][1]*R0[INDEX_PERTURBATIVE][1] + Exp_grad[1][2]*R0[INDEX_PERTURBATIVE][2];
  RTAU[2] = Exp_grad[2][0]*R0[INDEX_PERTURBATIVE][0]+Exp_grad[2][1]*R0[INDEX_PERTURBATIVE][1] + Exp_grad[2][2]*R0[INDEX_PERTURBATIVE][2];

  px[INDEX_PERTURBATIVE]=-Exp_gradT[0][0]*RTAU[0]-Exp_gradT[0][1]*RTAU[1] -Exp_gradT[0][2]*RTAU[2]; 
  py[INDEX_PERTURBATIVE]=-Exp_gradT[1][0]*RTAU[0]-Exp_gradT[1][1]*RTAU[1] -Exp_gradT[1][2]*RTAU[2]; 
  pz[INDEX_PERTURBATIVE]=-Exp_gradT[2][0]*RTAU[0]-Exp_gradT[2][1]*RTAU[1] -Exp_gradT[2][2]*RTAU[2]; 
	
  scra=sqrt(px[INDEX_PERTURBATIVE]*px[INDEX_PERTURBATIVE]+py[INDEX_PERTURBATIVE]*py[INDEX_PERTURBATIVE] + pz[INDEX_PERTURBATIVE]*pz[INDEX_PERTURBATIVE]);   
  px[INDEX_PERTURBATIVE]/=scra;  
  py[INDEX_PERTURBATIVE]/=scra;  
  pz[INDEX_PERTURBATIVE]/=scra;
		
} //end function control PERTURBATIVE
