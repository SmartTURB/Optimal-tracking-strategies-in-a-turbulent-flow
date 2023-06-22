/************************************************/
/*  Global Params                               */
/************************************************/
#define myprintf(...)  aamyprintf(__func__,__VA_ARGS__)

/* parameters */
#define NITER 500       /* number of Iterations */
#define NPART 6         /* number of particles. 0 = Lagrangian target, 1 = Tracer, 2 = Pure Pursuit, 3 = Surfing Control, 4 = Perturbative Optimal Control, 5 = Optimal Control (OC). */
#define Neps 200000     /* number of episodes */
#define N 3             /* dimension of the system */                             

#define NO_CONVERGENCE 0 /* convergence of the iterative algorithm */
#define OK_CONVERGENCE 1

#define NO_CAPTURE 0 /* capture or no capture with the OC */
#define OK_CAPTURE 1
#define UNDEF_CAPTURE 2 /* the capture is undefined for non-convergent episodes */

#define INDEX_TARGET 0  /* particle index */
#define INDEX_TRACER 1 
#define INDEX_PUREPURSUIT 2 
#define INDEX_SURFING 3 
#define INDEX_PERTURBATIVE 4
#define INDEX_OPTCTRL 5


/************************************************/
/*  Global Variables                            */
/************************************************/

// For all file
double *AA_history; // Gradients history used for each episode

// For h5read
double *full;  // to read gradients from a file
double *x_target, *y_target, *z_target; // lagrangian target trajectory

// For functions
double px_ott[NITER+1], py_ott[NITER+1], pz_ott[NITER+1]; // output of the FBSM --> Optimal Control 
double Rx_ott[NITER+1], Ry_ott[NITER+1], Rz_ott[NITER+1]; // output of the FBSM --> Optimal relative trajectory 
double vs[NPART]; // self-propulsion speed 
double px[NPART], py[NPART], pz[NPART]; // control variables
double grad_matrix[3][3], Exp_grad[3][3], Exp_gradT[3][3]; // matrices used to calculate heuristic strategies
double capture_distance;  // capture distance
int convergence, capture; // boolean variables

double gamma_lr,delta; /* hyper-parameters for the FBSM algorithm */
int max_iteration;
double J_PP, J_SC, J_PO, J_OC, minJ_OC; // performance index along Pure Pursuit (PP), Surfing Control (SC), Perturbative OC (PO) and Optimal Control (OC)
double capturetime_OC_FBSM; //optimal capture time calculated from FBSM algorithm (i.e. with the dynamic smoothed with tanh())

FILE *performance_index; /* to evaluate the performance index along the iteration of the FBSM algorithm */

// consts
const double dt; // time step (it is defined as the dt used to extrapolate gradients from Lagrangian trajectories)
const double alpha; // stiffness used to integrate the dynamics in the FBSM algorithm 
const double capture_weight; // parameter in the performance index - determines the weight of the cost integral 
const double tau_surfing;  // free parameter for the surfing policy 
const double tau_perturbative; // free parameter for the perturbative optimal control 

extern int me, NP;
#define AMIROOT (me==0)

#define time_declare(varprefix) double varprefix##_start;
#define time_start(varprefix) { varprefix##_start = MPI_Wtime (); }
#define time_end(varprefix, msg) if (AMIROOT) { double varprefix##_end = MPI_Wtime (); if (AMIROOT) fprintf(stderr, msg ", time:%f secs\n", varprefix##_end - varprefix##_start); }
