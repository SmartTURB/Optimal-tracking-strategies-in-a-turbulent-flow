///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// THIS CODE SOLVES THE OPTIMAL CONTROL PROBLEM OF CATCHING A DRIFTING TARGET IN A 3D TURBULENT FLOW /////////////
/// We assume that distance between the target and the pursuer is always below the Kolmogorv scale, eta=0.0043 ////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////                       
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include <time.h>

#include "common_vars.h"
double *full = NULL;

// consts
const double dt=0.00225; // time step (it is defined as the dt used to extrapolate gradients from Lagrangian trajectories)
const double alpha = 10; // stiffness used to integrate the dynamics in the FBSM algorithm 
const double capture_weight = 5.8*(0.004/100.)*(0.004/100.)*100; // parameter in the performance index - determines the weight of the cost integral 
const double tau_surfing = 0.014;  // free parameter
const double tau_perturbative = 0.03;
double time_stepdt;
int count;
int par_iep; // index of episode;
double dx, dy, dz; // to initialize particles 
int numTimes;

#include "functions_withFBSM_3D.h"
#include "functions_h5_read.h"
#include "mpi.h"

int me, NP;

/************************************************/
/*                MAIN                          */
/************************************************/
int main(int argc, char *argv[]) {
  int step, seme, i, ip, t, count_conv_iep, count_NOconv_iep, count_conv_iep_nocapture, count_conv_iep_capture, check_final_linear_regime;
  double Tc_PP, Tc_SC, Tc_PO, Tc_OC;
  double dd0, vel_ag;
  double scra;
  char nome[1024];
  double **dist_t_nocapture, **logdist_t_nocapture, **dist2_t_nocapture, *logdist_t_capture; 
  FILE  *capture_times, *final_linear_regime_distances, *final_distances, *dist_nocapture, *logdist_nocapture,*dist2_nocapture, *logdist_capture, *conv_andInitialCondition, *final_performance_index ;
  
  /* MPI Init here */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  
  double x[NPART], y[NPART], z[NPART]; // particles coordinates
  int c[NPART]; // boolean variables that determine if the particles reach the capture or not.
  // if c[i-particle]==NO_CAPTURE integrate the dynamics of the i-particle; otherwise the problem is considered solved and the dynamics of the i-particle is not integrated anymore.  
  
  if (argc !=3) {
    fprintf(stderr, "usage %s <vel_ag> <dd0> \n", argv[0]);
    fprintf(stderr, "         <vel_ag> is velocity of agents \n");
    fprintf(stderr, "         <dd0> is initial distance \n");
    return 1;
  }
  
  vel_ag=atof(argv[1]);
  dd0=atof(argv[2]);  
  capture_distance = dd0/100.;
  
  /* OPEN FILES */
  sprintf(nome,"file_output/proc%d_Convergence_and_InitialConditions.dat", me);
  conv_andInitialCondition = fopen(nome,"w");  
  sprintf(nome,"file_output/proc%d_capture_times.dat",me);
  capture_times = fopen(nome,"w");  
  sprintf(nome,"file_output/proc%d_performance_index.dat",me);
  performance_index = fopen(nome,"w");  
  sprintf(nome,"file_output/proc%d_final_performance_index.dat",me);
  final_performance_index = fopen(nome,"w");    
  sprintf(nome,"file_output/proc%d_final_distances.dat",me);
  final_distances = fopen(nome,"w");
  sprintf(nome,"file_output/proc%d_final_linear_regime_distances.dat",me);
  final_linear_regime_distances = fopen(nome,"w");
  
  /* *********************** Initialization ************************* */
  seme=6667;
  if (AMIROOT) fprintf(stderr, "#seme=%d\n",seme);
  srand48(seme);
  
  /* ALLOCATION */
  x_target = (double *) malloc((NITER+2)*sizeof(double));
  y_target = (double *) malloc((NITER+2)*sizeof(double));
  z_target = (double *) malloc((NITER+2)*sizeof(double));
  
  logdist_t_capture = (double *) malloc((NITER+1)*sizeof(double));
  dist_t_nocapture = (double **) malloc(NPART*sizeof (double *));
  logdist_t_nocapture =(double **) malloc(NPART*sizeof (double *));
  dist2_t_nocapture =(double **) malloc(NPART*sizeof (double *));
  AA_history = (double *) malloc(9*(NITER+2)*sizeof(double));
  
  for(i=0;i<NPART;i++){
    dist_t_nocapture[i] = (double *) malloc((NITER+1)*sizeof (double));
    logdist_t_nocapture[i] = (double *) malloc((NITER+1)*sizeof (double));
    dist2_t_nocapture[i] = (double *) malloc((NITER+1)*sizeof (double));
  }
  
  vs[INDEX_TARGET] = vs[INDEX_TRACER] = 0.0; //the first two particles are tracers 
  for (ip=INDEX_PUREPURSUIT; ip<NPART; ip++){
    vs[ip] = vel_ag;
  }
  
  for (step=0; step<NITER+1; step++){
    logdist_t_capture[step]=0.0;
    for (ip = 0; ip<NPART; ip++){
      dist_t_nocapture[ip][step]    = 0.0;
      logdist_t_nocapture[ip][step] = 0.0;
      dist2_t_nocapture[ip][step]   = 0.0;
    }
  }
  
  count_conv_iep = 0;
  count_NOconv_iep = 0;
  count_conv_iep_nocapture = 0;
  count_conv_iep_capture = 0;
  
  /* ***************** Write Info **************** */
  if (AMIROOT) {
    fprintf(stderr, "#Swimming in linear flow with gradients extrapolated from Lagrangian trajectories in 3d turbulent flows \n" );
    fprintf(stderr, "dt=%g",dt);
    fprintf(stderr, "#Number of PART =%d\n", (int) NPART);
    fprintf(stderr, "Initial distance =%g\n", dd0);
    fprintf(stderr, "Agent velocity =%g\n", vel_ag);
  }
  
  /* --------------------------------------------------------------------------------------------------------- */
  /* ****************************************** loop episodes ************************************************ */
  /* --------------------------------------------------------------------------------------------------------- */
  
  time_declare(totTime);
  time_start(totTime);
  
  for (par_iep = 1; par_iep<=Neps/NP*me; par_iep++) {
    dx=drand48()-0.5;
    dy=drand48()-0.5;
    dz=drand48()-0.5;
  }
  
  for (par_iep = 1; par_iep<=Neps/NP; par_iep++){
    int iep = par_iep + me*Neps/NP;    
    time_declare(op);    
    time_start(op);
    
    for(ip=0;ip<NPART;ip++)
      c[ip] = NO_CAPTURE;
    
    Tc_PP = Tc_SC = Tc_PO = Tc_OC = capturetime_OC_FBSM = 0;
    
    if(iep%5==0){
      fflush(conv_andInitialCondition);
      fflush(capture_times);
      fflush(final_distances);
      fflush(final_linear_regime_distances);
      if (AMIROOT) {
	fprintf(stderr,"STUDIED %d EPISODES \n", iep);
	fprintf(stderr, "number of converged episodes = %d \n", count_conv_iep);
	fprintf(stderr, "number of no-converged episodes  = %d \n", iep-1-count_conv_iep);
      }			
    }
    
    /* initialize particles position */
    x[0] = y[0] = z[0] = 0; // I work in a linear regime - one particle (the Target) in the target reference is always in the origin.
    dx=drand48()-0.5;
    dy=drand48()-0.5;
    dz=drand48()-0.5;
    for(ip=1;ip<NPART;ip++){
      x[ip] = x[0] + dd0 * dx/sqrt(dx*dx + dy*dy + dz*dz);
      y[ip] = y[0] + dd0 * dy/sqrt(dx*dx + dy*dy + dz*dz);
      z[ip] = z[0] + dd0 * dz/sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    /* read gradients from a file and put in AA_history */ 
    numTimes = traj_read("full_traj_tracers.h5", iep-1, 1);
    if (numTimes<1) { // Didn't read a traj
      fprintf(stderr,"numTimes=%d\n", numTimes);
      fprintf(stderr, "\n !!!!!! Traj database was not read. Aborting... \n");
      return -1;
    }    
    int isOk = traj_dump(iep-1, 1, NITER+1); 
    
    /* calculate Optimal Control troughth FBSM algorithm */
    gamma_lr=0.0005; //learning rate 
    delta=0.0001; //threshold for convergence
    max_iteration=20000; //max iteration for FBSM
    control_ForwardBackwardSweepMethod(x,y,z,iep);
    if(convergence==NO_CONVERGENCE){
      gamma_lr=0.00005; // change the hyperparameters
      delta=0.00001;
      max_iteration=200000;
      capturetime_OC_FBSM =0;
      control_ForwardBackwardSweepMethod(x,y,z,iep);
    }
    
    /* Update stats over episodes */
    if(convergence==NO_CONVERGENCE){
      capture = UNDEF_CAPTURE; // if OC does not converge -- > capture = UNDEF_CAPTURE
      count_NOconv_iep++;
    }
    if(convergence==OK_CONVERGENCE) {
      count_conv_iep++;
      if(capture==OK_CAPTURE) count_conv_iep_capture++;
      if(capture==NO_CAPTURE) count_conv_iep_nocapture++;
    }
		
    /* save convergence/capture/initial_conditions for each episode */
    fprintf(conv_andInitialCondition, "%d %d %d %g %g %g \n", iep, convergence, capture,
	    x[INDEX_OPTCTRL]-x[INDEX_TARGET], y[INDEX_OPTCTRL]-y[INDEX_TARGET], z[INDEX_OPTCTRL]-z[INDEX_TARGET]);


    /* ---------------------------------------------------------------------------------------- */
    /* -------------------------------- CASE NO CONVERGENCE ------------------------------------ */
    /* ---------------------------------------------------------------------------------------- */
    if (convergence==NO_CONVERGENCE){ // if OC does not converge I integrate heuristic policies to save the final performance index
      time_stepdt=0;
      
      for (step=0; step<NITER; step++){

	px[INDEX_TARGET] = py[INDEX_TARGET] = px[INDEX_TRACER] = py[INDEX_TRACER] = 0; // tracers
	control_purepursuit(x,y,z);  
	control_surfing(x,y,z,step);
	control_perturbative(x,y,z,step); 		   
	// Pure Pursuit
	if (c[INDEX_PUREPURSUIT] == OK_CAPTURE){ 
	  px[INDEX_PUREPURSUIT] = py[INDEX_PUREPURSUIT] = 0; 
	  if(Tc_PP==0){ 
	    Tc_PP = (step-1)*dt;
	  }
	}
	// Surfing Control
	if (c[INDEX_SURFING] == OK_CAPTURE){
	  px[INDEX_SURFING] = py[INDEX_SURFING] = 0;
	  if(Tc_SC==0){
	    Tc_SC = (step-1)*dt;
	  }
	}
	// Perturbative Optimal control
	if (c[INDEX_PERTURBATIVE] == OK_CAPTURE){
	  px[INDEX_PERTURBATIVE] = py[INDEX_PERTURBATIVE] = 0;
	  if(Tc_PO==0){
	    Tc_PO = (step-1)*dt;
	  }		     
	}
	
	/* integrate particles position */ 
	rk4_advance_particles(x,y,z, step, c);
	time_stepdt = time_stepdt + dt;       
      } 

      scra=dist(x,y,z,INDEX_PUREPURSUIT);
      J_PP=scra*scra+capture_weight*Tc_PP;
      scra=dist(x,y,z,INDEX_SURFING);
      J_SC = scra*scra + capture_weight*Tc_SC;
      scra=dist(x,y,z,INDEX_PERTURBATIVE);
      J_PO = scra*scra + capture_weight*Tc_PO;
      
      fprintf(final_performance_index, "%d %d %d %g %g %g %g %g \n", iep, convergence, capture, J_PP, J_SC, J_PO, J_OC, minJ_OC);
      
    } 
    /* ---------------------------------------------------------------------------------------- */
    /* -------------------------------- END NO CONVERGENCE ------------------------------------ */
    /* ---------------------------------------------------------------------------------------- */


    /* ---------------------------------------------------------------------------------------- */
    /* -------------------------------- CASE OK CONVERGENCE ------------------------------------ */
    /* ---------------------------------------------------------------------------------------- */
    if (convergence == OK_CONVERGENCE) { // Investigate the episode only if the FBSM algorithm for the OC is convergent
      // Note that with further fine-tuning of the hyper-parameters involved in the iterative process, one might still improve the convergence up to 100%.
      Tc_PP = Tc_SC = Tc_PO = Tc_OC = 0;
      check_final_linear_regime=0;
      
      time_stepdt=0;
      for (step=0; step<NITER; step++) {
	
	px[INDEX_TARGET] = py[INDEX_TARGET] = pz[INDEX_TARGET] = px[INDEX_TRACER] = py[INDEX_TRACER] = pz[INDEX_TRACER] = 0; // tracers			
	control_purepursuit(x,y,z);  
	control_surfing(x,y,z,step);
	control_perturbative(x,y,z,step);
	px[INDEX_OPTCTRL] = px_ott[step];
	py[INDEX_OPTCTRL] = py_ott[step];
	pz[INDEX_OPTCTRL] = pz_ott[step];
	/* look at the capture time for each policy */
	// Pure Pursuit // 
	if (c[INDEX_PUREPURSUIT] == OK_CAPTURE){ 
	  px[INDEX_PUREPURSUIT] = py[INDEX_PUREPURSUIT] = pz[INDEX_PUREPURSUIT] = 0; 
	  if(Tc_PP==0){
	    Tc_PP = (step-1)*dt;
	  }
	}
	// Surfing Control
	if (c[INDEX_SURFING] == OK_CAPTURE){
	  px[INDEX_SURFING] = py[INDEX_SURFING] = pz[INDEX_SURFING] = 0;
	  if(Tc_SC==0){
	    Tc_SC = (step-1)*dt;
	  }
	}
	// Perturbative Optimal control
	if (c[INDEX_PERTURBATIVE] == OK_CAPTURE){
	  px[INDEX_PERTURBATIVE] = py[INDEX_PERTURBATIVE] = pz[INDEX_PERTURBATIVE] = 0;
	  if(Tc_PO==0){
	    Tc_PO = (step-1)*dt;
	  }
	}
	// Optimal Control
	if (c[INDEX_OPTCTRL] == OK_CAPTURE){
	  px[INDEX_OPTCTRL] = py[INDEX_OPTCTRL] = pz[INDEX_OPTCTRL] = 0;
	  if(Tc_OC==0){
	    Tc_OC = (step-1)*dt;
	  }
	}
		       					
	/* save statistical info when the agent optimally-controlled does not capture the target */
	if (capture == NO_CAPTURE){ 
	  for(ip=1; ip<NPART; ip++){
	    scra=dist(x,y,z,ip);
	    dist_t_nocapture[ip][step] +=  scra;
	    logdist_t_nocapture[ip][step] += log10(scra);
	    dist2_t_nocapture[ip][step] += scra*scra;
	  }			
	}
	if (capture == OK_CAPTURE){ /* save statistical info only for the tracer when OC captures */
	  scra=dist(x,y,z,INDEX_TRACER);
	  logdist_t_capture[step]+=log10(scra);
	}
			
	/* integrate particles position */ 
	rk4_advance_particles(x,y,z, step, c);
	time_stepdt = time_stepdt + dt;
			
	if (capture == NO_CAPTURE){
	  scra = dist(x,y,z,INDEX_OPTCTRL);
	  if (check_final_linear_regime==0 && scra >= 10*0.0043){ //10*0.0043 \sim 10 \eta (border of the linear regime)
	    fprintf(final_linear_regime_distances, "%d %g %g %g %g %g \n",iep, time_stepdt, dist(x,y,z,INDEX_PUREPURSUIT), dist(x,y,z,INDEX_SURFING), dist(x,y,z,INDEX_PERTURBATIVE), dist(x,y,z,INDEX_OPTCTRL));
	    check_final_linear_regime=1;
	  } 
	}
	
      } //end loop step     

      if (capture == NO_CAPTURE){		  
	if (check_final_linear_regime==0){
	  fprintf(final_linear_regime_distances, "%d %g %g %g %g %g \n",iep, time_stepdt, dist(x,y,z,INDEX_PUREPURSUIT), dist(x,y,z,INDEX_SURFING), dist(x,y,z,INDEX_PERTURBATIVE), dist(x,y,z,INDEX_OPTCTRL));
	  check_final_linear_regime=1;
	} 
      }
      
      /* save final state */
      /* save statistical properties along the no-capture episodes */
      if (capture == NO_CAPTURE){ 
	for(ip=1; ip<NPART; ip++){
	  scra=dist(x,y,z,ip);
	  dist_t_nocapture[ip][step] +=  scra;
	  logdist_t_nocapture[ip][step] += log10(scra);
	  dist2_t_nocapture[ip][step] += scra*scra;
	}

	sprintf(nome,"file_output/proc%d_mean_dist_nocaptureepisodes.dat",me);
	dist_nocapture = fopen(nome,"w");
	sprintf(nome,"file_output/proc%d_mean_logdist_nocaptureepisodes.dat",me);
	logdist_nocapture =  fopen(nome,"w");
	sprintf(nome,"file_output/proc%d_mean_dist2_nocaptureepisodes.dat",me);
	dist2_nocapture =  fopen(nome,"w");

	for (step=0; step<NITER+1; step++){
	  fprintf(dist_nocapture, "%g %g %g %g %g %g  \n", step*dt,  dist_t_nocapture[INDEX_TRACER][step],\
		  dist_t_nocapture[INDEX_PUREPURSUIT][step], dist_t_nocapture[INDEX_SURFING][step], dist_t_nocapture[INDEX_PERTURBATIVE][step], dist_t_nocapture[INDEX_OPTCTRL][step]);
	  fprintf(logdist_nocapture, "%g %g %g %g %g %g \n", step*dt,logdist_t_nocapture[INDEX_TRACER][step], \
		  logdist_t_nocapture[INDEX_PUREPURSUIT][step],logdist_t_nocapture[INDEX_SURFING][step], logdist_t_nocapture[INDEX_PERTURBATIVE][step],	\
		  logdist_t_nocapture[INDEX_OPTCTRL][step]);
	  fprintf(dist2_nocapture, "%g %g %g %g %g %g \n", step*dt,dist2_t_nocapture[INDEX_TRACER][step], \
		  dist2_t_nocapture[INDEX_PUREPURSUIT][step], dist2_t_nocapture[INDEX_SURFING][step], dist2_t_nocapture[INDEX_PERTURBATIVE][step], dist2_t_nocapture[INDEX_OPTCTRL][step]);
	}
	fclose(dist_nocapture);
	fclose(logdist_nocapture);
	fclose(dist2_nocapture);
		  
      }// end if NO_CAPTURE


      if (capture == OK_CAPTURE){	
	scra=dist(x,y,z,INDEX_TRACER);
	logdist_t_capture[NITER]+=log10(scra);
	sprintf(nome,"file_output/proc%d_mean_logdist_tracer_captureepisodes.dat",me);
	logdist_capture =  fopen(nome,"w");
	for(step=0; step<NITER+1; step++){
	  fprintf(logdist_capture, "%g %g \n", step*dt, logdist_t_capture[step]);
	}
	fclose(logdist_capture);
      }


      /* save capture time */                                                                                                                             
      if(Tc_PP == 0) Tc_PP = NITER*dt; 	// if the capture time is still 0 it means that the agent does not capture.
      if(Tc_PO == 0) Tc_PO = NITER*dt; 	// in this case put the time horizon as a capture time (-->FAILURE). 
      if(Tc_SC == 0) Tc_SC = NITER*dt;		
      if(Tc_OC == 0) Tc_OC = NITER*dt;
      if(capturetime_OC_FBSM==0) capturetime_OC_FBSM = NITER*dt;
      fprintf(capture_times, "%g %g %g %g %g \n", Tc_PP, Tc_SC, Tc_PO, Tc_OC, capturetime_OC_FBSM);

      /* save final distances (not necessary within the limit of 10\eta)*/
      fprintf(final_distances, "%g %g %g %g %g %d \n", dist(x,y,z,INDEX_TRACER), dist(x,y,z,INDEX_PUREPURSUIT), dist(x,y,z,INDEX_SURFING), dist(x,y,z,INDEX_PERTURBATIVE), dist(x,y,z,INDEX_OPTCTRL),iep);

      scra=dist(x,y,z,INDEX_PUREPURSUIT);
      J_PP = scra*scra + capture_weight*Tc_PP; 
      scra=dist(x,y,z,INDEX_SURFING);
      J_SC = scra*scra + capture_weight*Tc_SC; 
      scra=dist(x,y,z,INDEX_PERTURBATIVE);
      J_PO = scra*scra + capture_weight*Tc_PO;
      fprintf(final_performance_index, "%d %d %d %g %g %g %g %g \n", iep, convergence, capture, J_PP, J_SC, J_PO, J_OC, minJ_OC);

    } // end if OK_CONVERGENCE	
  }
  time_end(totTime, "Total Time for episodes]");

  /* --------------------------------------------------------------------------------------------------------- */
  /* *************************************** end loop episodes *********************************************** */
  /* --------------------------------------------------------------------------------------------------------- */
	
  int all_count_conv_iep, all_count_NOconv_iep, all_count_conv_iep_capture, all_count_conv_iep_nocapture;

  MPI_Reduce(&count_conv_iep, &all_count_conv_iep, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&count_NOconv_iep, &all_count_NOconv_iep, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_conv_iep_capture, &all_count_conv_iep_capture, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_conv_iep_nocapture, &all_count_conv_iep_nocapture, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (AMIROOT) {
    fprintf(stderr, "number of converged episodes = %d\n", all_count_conv_iep);
    fprintf(stderr, "number of no-converged episodes = %d\n", all_count_NOconv_iep);
    fprintf(stderr, "number of capture episodes = %d \n", all_count_conv_iep_capture);
    fprintf(stderr, "number of no-capture episodes = %d \n", all_count_conv_iep_nocapture);
  }

  MPI_Finalize ();
	
  return 0;
}
