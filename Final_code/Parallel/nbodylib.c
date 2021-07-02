#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#include "nbodylib.h"


void checkr(int r, char *txt)
{
	if (r != MPI_SUCCESS)
	{
		fprintf(stderr, "Error: %s\n", txt);
		exit(-1);
	}
}

// Indicates the rank of the processor

int proc()
{
	int r, rank;
	r = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	checkr(r, "proc");
	return (rank);
}

// Indicates the number of processors being used

int nproc()
{
	int r, size;
	r = MPI_Comm_size(MPI_COMM_WORLD, &size);
	checkr(r, "nproc");
	return (size);
}

// Worksplit

void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end)
{
	// Number of tasks
	int ntask = end - start + 1;
	// Number of tasks per processor
	int interval = ntask / nproc;
	// Tasks left
	int remainder = ntask % nproc;

	if (ntask < nproc)
	{
		printf("Less tasks than processors\n");
		exit(-1);
	}
	if (remainder != 0)
	{
		if (proc < remainder)
		{
			*mystart = start + proc * (interval + 1);
			*myend = *mystart + interval;
		}
		else
		{
			*mystart = start + remainder * (interval + 1) + (proc - remainder) * interval;
			*myend = *mystart + interval - 1;
		}
	}
	else
	{
		*mystart = start + proc * interval;
		*myend = *mystart + (interval - 1);
	}
}

void NBody(double t, double* var, int n, double* varp) {

    /*
	  Body 1: Perturbated
		  var[0] = rx   var[3] = vx
		  var[1] = ry   var[4] = vy
		  var[2] = rz   var[5] = vz
    */
    // varp = prime (derivative of var | var prime)

    int r; 		// for error checking
	int start = 0;
    MPI_Status st;

    if(proc() != 0) {
		start = planet2proc[proc()-1] + 1;
	  }

    double* y;
    y = (double*)malloc(3*NBODY*sizeof(double));

    if (proc() != 0) {

        double* pos; // vectors of positions
    	pos = (double*)malloc(3*(BODIES(proc()))*sizeof(double));

		for (int body = 0; body < (BODIES(proc())); body++)
		{
			pos[3*body] 	= VAR(body,0);
			pos[3*body + 1] = VAR(body, 1);
			pos[3*body + 2] = VAR(body, 2);
		}

       r = MPI_Ssend(pos, 3*(BODIES(proc())), MPI_DOUBLE, 0, proc() + nproc(), MPI_COMM_WORLD);
	   checkr(r, "send");

       free(pos);

       r = MPI_Recv(y, 3*NBODY, MPI_DOUBLE, 0, proc(), MPI_COMM_WORLD, &st);
	   checkr(r, "receive");
    }
    else { // If rank == 0

    	double* var_loc;
        
		for (int y_counter = 0; y_counter <= planet2proc[0]; y_counter++) {
			Y(y_counter,0) = VAR(y_counter,0);
			Y(y_counter,1) = VAR(y_counter,1);
			Y(y_counter,2) = VAR(y_counter,2);
		}

		for (int xproc = 1; xproc < nproc(); xproc++) {

			var_loc = (double*)malloc(3*(BODIES(xproc))*sizeof(double));

			r = MPI_Recv(var_loc, 3*(BODIES(xproc)), MPI_DOUBLE, xproc, xproc + nproc(), MPI_COMM_WORLD, &st);
			checkr(r, "receive");

			for (int p = 0; p < (3*(BODIES(xproc))); p++) {
				Y(planet2proc[xproc - 1] + 1,p) = var_loc[p];
			}

			free(var_loc);

			r = MPI_Ssend(y, 3*NBODY,MPI_DOUBLE, xproc, xproc, MPI_COMM_WORLD);
			checkr(r, "send");
		}
	}

    // Each planet with respect to other
    for (int loc_body = 0; loc_body<(n/6); loc_body++) {

        VARP(loc_body, 0) = 0;
        VARP(loc_body, 1) = 0;
        VARP(loc_body, 2) = 0;
        VARP(loc_body, 3) = 0;
        VARP(loc_body, 4) = 0;
        VARP(loc_body, 5) = 0;

        // For all relations with other particles
        for (int glob_body = 0; glob_body<NBODY; glob_body++) {

            // Do not include the particle with itself
            if (IGLOB(loc_body,start) == glob_body) continue;

            // Distance modulus
            double r = sqrt(DIST(IGLOB(loc_body,start),glob_body,0)*DIST(IGLOB(loc_body,start),glob_body,0) +
                            DIST(IGLOB(loc_body,start),glob_body,1)*DIST(IGLOB(loc_body,start),glob_body,1) +
                            DIST(IGLOB(loc_body,start),glob_body,2)*DIST(IGLOB(loc_body,start),glob_body,2));

            // Perturbated
            VARP(loc_body, 3) += (G*MASS[glob_body] / r / r / r) * DIST(IGLOB(loc_body,start),glob_body,0);  // ax
            VARP(loc_body, 4) += (G*MASS[glob_body] / r / r / r) * DIST(IGLOB(loc_body,start),glob_body,1);  // ay
            VARP(loc_body, 5) += (G*MASS[glob_body] / r / r / r) * DIST(IGLOB(loc_body,start),glob_body,2);  // az

        }

        VARP(loc_body, 0) = VAR(loc_body,3);                                  // vx
        VARP(loc_body, 1) = VAR(loc_body,4);                                  // vy
        VARP(loc_body, 2) = VAR(loc_body,5);                                  // vz

    }

	free(y);

}
