#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdarg.h>
#include "hdf5.h"

#include "common_vars.h"
int matrix_index(int i, int j, int t);

typedef struct {
  double xo, yo, zo;
  double uxo, uyo, uzo;
  double axo, ayo, azo;
  unsigned int name;
} particle_lagr_fast_type;

#define DEBUG_LOG 0

int aamyprintf(const char *fname, int priority, const char *format, ...)
{

  if (AMIROOT) {

    fprintf(stderr,"[%s] ",fname);
    va_list args;
    va_start(args, format);
    
    vfprintf(stderr, format, args);
    
    va_end(args);
  }
}

void myAbort(int ierr) {
	exit(ierr);
}
	
hid_t MakeTypeOld() {
   hid_t H5_TYPE_IN_USE = H5T_NATIVE_DOUBLE;

   hid_t hdf5_pop_type = H5Tcreate (H5T_COMPOUND, sizeof(particle_lagr_fast_type));

   H5Tinsert(hdf5_pop_type, "xo", HOFFSET(particle_lagr_fast_type, xo), H5_TYPE_IN_USE);
   H5Tinsert(hdf5_pop_type, "yo", HOFFSET(particle_lagr_fast_type, yo), H5_TYPE_IN_USE);
   H5Tinsert(hdf5_pop_type, "zo", HOFFSET(particle_lagr_fast_type, zo), H5_TYPE_IN_USE);

   H5Tinsert(hdf5_pop_type, "uxo", HOFFSET(particle_lagr_fast_type, uxo), H5_TYPE_IN_USE);
   H5Tinsert(hdf5_pop_type, "uyo", HOFFSET(particle_lagr_fast_type, uyo), H5_TYPE_IN_USE);
   H5Tinsert(hdf5_pop_type, "uzo", HOFFSET(particle_lagr_fast_type, uzo), H5_TYPE_IN_USE);

   H5Tinsert(hdf5_pop_type, "axo", HOFFSET(particle_lagr_fast_type, axo), H5_TYPE_IN_USE);
   H5Tinsert(hdf5_pop_type, "ayo", HOFFSET(particle_lagr_fast_type, ayo), H5_TYPE_IN_USE);
   H5Tinsert(hdf5_pop_type, "azo", HOFFSET(particle_lagr_fast_type, azo), H5_TYPE_IN_USE);

   H5Tinsert (hdf5_pop_type, "name", HOFFSET (particle_lagr_fast_type, name), H5T_NATIVE_UINT);

  return hdf5_pop_type;
}



int traj_read(const char *fname, long startPart, long cntPart) {
  static int isFirst = 1;
      
   /* Open a new file using default properties. */
   hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
   int file_exist = file_id > 0;

   hid_t group = H5Gopen2(file_id, "/", H5P_DEFAULT);

   hid_t hdf5_pop_type = MakeTypeOld();

   int part_num = 0, numTimes = 0;

   /* Open the dataset. */
   if (file_exist) {
        hid_t dataset = H5Dopen(group, "traj3d", H5P_DEFAULT);
        if (dataset < 0) {
            myprintf(DEBUG_LOG, "Not able to open dataset \n");
            myAbort(3);
        }

        hid_t filespace = H5Dget_space(dataset);

        int rank = H5Sget_simple_extent_ndims(filespace);

        hsize_t dimsr[3];
        hid_t status = H5Sget_simple_extent_dims(filespace, dimsr, NULL);
        if (isFirst) {
            isFirst = 0;

            myprintf(DEBUG_LOG, "dataset) rank:%d dims1=%d dims2=%d dims3=%d \n", rank, dimsr[0], dimsr[1], dimsr[2]);
        }

        part_num = dimsr[0];
        if (part_num <= 0) {
            myprintf(DEBUG_LOG, "!!! particle number(%d) <=0. Bye.\n", part_num);
            myAbort(4);
        }

        numTimes = dimsr[1];
        if (numTimes <= 0) {
            myprintf(DEBUG_LOG, "!!! particle times(%d) <=0. Bye.\n", numTimes);
            myAbort(4);
        }

        /* Define hyperslab in the dataset. */
        hsize_t count[3];
        count[0] = cntPart;
        count[1] = numTimes;
        count[2] = 18;
        hid_t memspace = H5Screate_simple(3, count, NULL);


        hsize_t fileCount[3];
        fileCount[0] = cntPart;
        fileCount[1] = numTimes;
        fileCount[2] = 18;
        
        hsize_t fileSt[3];
        fileSt[0] = startPart;
        fileSt[1] = 0;
        fileSt[2] = 0;

        status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, fileSt, NULL, fileCount, NULL);
        //fprintf(stderr, "status fileSel = %d \n", status);

        //fprintf(stderr, "Reading Start particle name=%ld, num particles=%ld numTimes=%d \n", startPart,cntPart, numTimes);
        full = (double *) malloc(sizeof(double) * cntPart * numTimes * 18);

        /* Read the dataset in memory */
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, full);
        //fprintf(stderr, "status = %d \n", status);

        H5Tclose(hdf5_pop_type);

        H5Dclose(dataset);
        H5Sclose(memspace);
        H5Sclose(filespace);
   }

   H5Gclose(group);
   H5Fclose(file_id);

   return numTimes;
}

int traj_dump(long startPart, long cntPart, int TimeHorizon) {
   double *p = full;
   int j, ee;
   long i;
   double velacc_not_toread;
    
   for (i=0; i<cntPart; i++) {      
     for (j=0; j<TimeHorizon; j++) {        
       x_target[j] = p[0];
       y_target[j] = p[1];
       z_target[j] = p[2];

       // I skip the velocity and acceleration info (from index 3 to index 8)
	
       AA_history[matrix_index(0,0,j)] = p[9];
       AA_history[matrix_index(0,1,j)] = p[12];
       AA_history[matrix_index(0,2,j)] = p[15];
       AA_history[matrix_index(1,0,j)] = p[10];
       AA_history[matrix_index(1,1,j)] = p[13];
       AA_history[matrix_index(1,2,j)] = p[16];
       AA_history[matrix_index(2,0,j)] = p[11];
       AA_history[matrix_index(2,1,j)] = p[14];
       AA_history[matrix_index(2,2,j)] = p[17];
	
       p += 18;
     }
   }       
   return 1;
}
