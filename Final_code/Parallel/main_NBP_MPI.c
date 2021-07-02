/* N-body problem
Date: 12/04/2021
Author/s:
  - Angel Pan Du
  - Alba Molina Cuadrado
  - Iván Sermanoukian Molina
  - Yi Qiang Ji Zhang
Subject: High Performance Computing for Aerospace Engineering
Professor: Manel Soria & Arnau Miro

// Problem statement
------------------------------------------------------------------------
N Body Problem with MPI
-------------------------------------------------------------------------
*/


#include <stdio.h>	           // Standard Input and Output Library
#include <stdlib.h>            // Library for defining functions to perform general functions (malloc())
#include <math.h>	             // Math library
#include <string.h>            // Library for manipulating arrays of characters

#include "mpi.h"               // MPI Library
#include "nbodylib.h"          // N-body functions and actions
#include "RK_MPI.h"            // Runge Kutta library

// Global variables
int* planet2proc;
int array_size;

// Global vector of mass
double MASS[NBODY] = {M0,M1,M2/*,M3,M4,M5,M6,M7,M8,M9*/};

// Main
int main(int argc, char* argv[]) {

    int r; // For error checking

    // Start MPI
    r = MPI_Init(&argc, &argv);
    checkr(r, "Initiate");
    MPI_Status st;

    // Check rank and size of the processors
    int rank, size;
    rank = proc();
    size = nproc();

    array_size = size;
    planet2proc = (int*)malloc(sizeof(int)*array_size);

    // If the processor is higher than the number of processors, exit
    if (rank >= size) {exit(-1);}

    // Declare worksplit variables
    int mystart,mylocstart;
    int myend, mylocend;
    int bodies;
    for (int xproc=0; xproc<size; xproc++) {
      worksplit(&mystart, &myend, xproc, size, 0, NBODY - 1);
      if (xproc == rank) {
        bodies = myend - mystart + 1;
        mylocstart = mystart;
        mylocend = myend;
      }
      planet2proc[xproc] = myend;
    }

    // Data calculations

    double V1 = Vx(M0,R1);  // Mercury velocity
    double T1 = Tx(M0,R1);  // Mercury period
    double V2 = Vx(M0,R2);  // Venus velocity
    double T2 = Tx(M0,R2);  // Venus period
    // double V3 = Vx(M0,R3);  // Earth velocity
    // double T3 = Tx(M0,R3);  // Earth period
    // double V4 = Vx(M0,R4);  // Mars velocity
    // double T4 = Tx(M0,R4);  // Mars period
    // double V5 = Vx(M0,R5);  // Jupiter velocity
    // double T5 = Tx(M0,R5);  // Jupiter period
    // double V6 = Vx(M0,R6);  // Saturn velocity
    // double T6 = Tx(M0,R6);  // Saturn period
    // double V7 = Vx(M0,R7);  // Uranus velocity
    // double T7 = Tx(M0,R7);  // Uranus period
    // double V8 = Vx(M0,R8);  // Neptune velocity
    // double T8 = Tx(M0,R8);  // Neptune period
    // double V9 = Vx(M0,R9);  // Pluto velocity
    // double T9 = Tx(M0,R9);  // Pluto period

    // Integration bounds and initial solution
    double tspan[2] = {0.,5*T2};

    // All initial conditions
    double *y0_loc;
	  y0_loc = (double *)malloc(bodies*6*sizeof(double));

    // Global initial conditions for all particles
    for(int planet_p = 6*mylocstart; planet_p < (6*(mylocend + 1)); planet_p++) {

        double y0_temp[N] = {0., 0., 0., 0., 0., 0.,    /* Sun */
                            R1, 0., 0., 0., V1, 0.,     /* Mercury */
                            0., R2, 0., -V2, 0., 0.,    /* Venus */
                            // 0., R3, 0., -V3, 0., 0., /* Earth */
                            // 0., R4, 0., -V4, 0., 0., /* Mars */
                            // 0., R5, 0., -V5, 0., 0., /* Jupiter */
                            // 0., R6, 0., -V6, 0., 0., /* Saturn */
                            // 0., R7, 0., -V7, 0., 0., /* Uranus */
                            // 0., R8, 0., -V8, 0., 0., /* Neptune */
                            // 0., R9, 0., -V9, 0., 0.  /* Pluto */
                        };
        y0_loc[planet_p-6*mylocstart] = y0_temp[planet_p];
    }

    // Runge-Kutta parameters
    RK_PARAM rkp = rkparams(tspan);
    rkp.h0 = 24.*3600.;
    rkp.eps = 1e-11;

    // Launch the integrator
    RK_OUT rko = odeRK("hiroshi912", NBody, tspan, y0_loc, 6*bodies, &rkp);

    // Finish
    printf("error = %.2e with %d steps\n", rko.err, rko.n);

	  char str[15+nproc()];

    for(int xproc = 0; xproc < nproc(); xproc++) {
      if(xproc == proc()) {
        sprintf(str, "Output%d.csv", xproc);
        writerkout(str, &rko, 6*bodies);
      }
    }

    freerkout(&rko);
    free(y0_loc);
    free(planet2proc);

    // Finalise MPI
    MPI_Finalize();

    // End of the program
    exit(0);
}
