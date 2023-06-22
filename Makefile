MPI-INCLUDES= -I/usr/include/openmpi-x86_64
HDF5-INCLUDES= -I/usr/include/hdf5/openmpi

LIB-HDF5= /usr/lib/x86_64-linux-gnu/libhdf5_openmpi.so


CC = mpicc -O3
LINK = mpicc -O3

lagrangianControl.x:	lagrangianControl.o functions_withFBSM_3D.o functions_h5_read.o
	$(LINK) lagrangianControl.o functions_withFBSM_3D.o functions_h5_read.o -o lagrangianControl.x $(LIB-HDF5) -lgsl -lm

clean:
	rm -f lagrangianControl.x lagrangianControl.o functions_withFBSM_3D.o functions_h5_read.o

functions_h5_read.o:	functions_h5_read.c common_vars.h
	$(CC) functions_h5_read.c -c $(HDF5-INCLUDES) $(MPI-INCLUDES)

functions_withFBSM_3D.o:	functions_withFBSM_3D.c common_vars.h
	$(CC) functions_withFBSM_3D.c -c $(HDF5-INCLUDES)

lagrangianControl.o:	lagrangianControl.c common_vars.h
	$(CC) lagrangianControl.c $(MPI-INCLUDES) $(HDF5-INCLUDES) -c

