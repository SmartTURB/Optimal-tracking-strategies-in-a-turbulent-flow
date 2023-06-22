double dist(double *x,double *y, double *z, int i);

void control_ForwardBackwardSweepMethod(double*x, double*y, double *z, int iep);
void control_purepursuit(double *x,double *y, double *z);
void control_surfing(double *x,double *y, double *z, int t);
void control_perturbative(double *x,double *y, double *z, int t);

void rk4_advance_particles(double *x,double *y, double *z, int step, int c[NPART]);

int matrix_index(int i, int j, int t);
