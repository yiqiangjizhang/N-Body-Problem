/* Three-body problem
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
Three Body Problem in Serie
-------------------------------------------------------------------------
*/

#include <stdio.h>	           // Standard Input and Output Library
#include <stdlib.h>            // Library for defining functions to perform general functions (malloc())
#include <math.h>	           // Math library
#include <string.h>            // Library for manipulating arrays of characters

#include "RK.h"                // Runge Kutta library


#define PI 4.*atan(1)          // Pi definition

// Parameters Sun - Earth
#define G     6.674e-11*1.e-9  // [km^3/kg/s^2]
#define M0    1.9891e30        // [kg] Sun's mass
#define M1    3.30104e23       // [kg] Mercury's mass
#define M2    4.86732e24       // [kg] Venus' mass
#define M3    5.97219e24       // [kg] Earth's mass 
#define M4    6.41693e23       // [kg] Mars' mass
#define M5    1.89813e27       // [kg] Jupiter's mass
#define M6    5.68319e26       // [kg] Saturn's mass
#define M7    8.68103e25       // [kg] Uranus' mass
#define M8    1.02410e26       // [kg] Neptune's mass
#define M9    1.46e22          // [kg] Pluto's mass

#define R1    57909227         // [km] Mercury distance from Sun
#define R2    108209475        // [km] Venus distance from Sun
#define R3    149598262        // [km] Earth distance from Sun
#define R4    227943824 	   // [km] Mars distance from Sun
#define R5    778340821 	   // [km] Jupiter distance from Sun
#define R6    1426666422 	   // [km] Saturn distance from Sun
#define R7    2870658186 	   // [km] Uranus distance from Sun
#define R8    4498396441       // [km] Neptune distance from Sun
#define R9    5906.4e6 	       // [km] Pluto distance from Sun

#define NBODY 10               // Number of bodies
#define N     6*NBODY          // Number of initial conditions

#define Vx(Mx,Rx) sqrt(G * Mx / Rx)
#define Tx(Mx,Rx) 2*PI/sqrt(G * Mx / Rx / Rx / Rx)

// Macros
#define VAR(b,p) var[6*(b) + (p)]              // b: body, p: parameter, 
#define VARP(b,p) varp[6*(b) + (p)]           // b: body, p: parameter, 
#define DIST(b1,b2,p) (VAR(b2,p) - VAR(b1,p))  // b1: body 1, b2: body 2 p: (0,1,2) = (x,y,z)

// Prototypes
void NBody(double t, double* var, int n, double* varp);

// Main
int main(int argc, char* argv[]) {

    double V1 = Vx(M0,R1);  // Mercury velocity
    double T1 = Tx(M0,R1);  // Mercury period
    double V2 = Vx(M0,R2);  // Venus velocity
    double T2 = Tx(M0,R2);  // Venus period
    double V3 = Vx(M0,R3);  // Earth velocity
    double T3 = Tx(M0,R3);  // Earth period
    double V4 = Vx(M0,R4);  // Mars velocity
    double T4 = Tx(M0,R4);  // Mars period
    double V5 = Vx(M0,R5);  // Jupiter velocity
    double T5 = Tx(M0,R5);  // Jupiter period
    double V6 = Vx(M0,R6);  // Saturn velocity
    double T6 = Tx(M0,R6);  // Saturn period
    double V7 = Vx(M0,R7);  // Uranus velocity
    double T7 = Tx(M0,R7);  // Uranus period
    double V8 = Vx(M0,R8);  // Neptune velocity
    double T8 = Tx(M0,R8);  // Neptune period
    double V9 = Vx(M0,R9);  // Pluto velocity
    double T9 = Tx(M0,R9);  // Pluto period
   
    // Integration bounds and initial solution
    double tspan[2] = {0.,T9};
   
    double y0[N] = {0., 0., 0., 0., 0., 0.,  /* Sun */ 
                    R1, 0., 0., 0., V1, 0.,  /* Mercury */ 
                    0., R2, 0., -V2, 0., 0., /* Venus */
                    0., R3, 0., -V3, 0., 0., /* Earth */
                    0., R4, 0., -V4, 0., 0., /* Mars */
                    0., R5, 0., -V5, 0., 0., /* Jupiter */
                    0., R6, 0., -V6, 0., 0., /* Saturn */
                    0., R7, 0., -V7, 0., 0., /* Uranus */
                    0., R8, 0., -V8, 0., 0., /* Neptune */
                    0., R9, 0., -V9, 0., 0.  /* Pluto */
                    };
            
    // Runge-Kutta parameters
    RK_PARAM rkp = rkparams(tspan);
    rkp.h0 = 24.*3600.;
    rkp.eps = 1e-11;

    // Launch the integrator
    RK_OUT rko = odeRK("hiroshi912", NBody, tspan, y0, N, &rkp);

    // Finish
    printf("error = %.2e with %d steps\n", rko.err, rko.n);
    
    writerkout("output.csv", &rko, N);

    freerkout(&rko);


    // End of the program
    exit(0);
}

void NBody(double t, double* var, int n, double* varp) {

    /*
	  Body 1: Perturbated
		  var[0] = rx   var[3] = vx
		  var[1] = ry   var[4] = vy
		  var[2] = rz   var[5] = vz
	  Body 2: Perturber
		  var[6] = rx   var[9]  = vx
		  var[7] = ry   var[10] = vy
		  var[8] = rz   var[11] = vz
    */
    // varp = prime (derivative of var | var prime)

    // Vector of mass
    double M[NBODY] = {M0,M1,M2,M3,M4,M5,M6,M7,M8,M9};

    // Each planet with respect to other
    for (int i = 0; i<NBODY; i++) {   

        VARP(i, 0) = 0;  
        VARP(i, 1) = 0;
        VARP(i, 2) = 0;
        VARP(i, 3) = 0;
        VARP(i, 4) = 0;
        VARP(i, 5) = 0;   

        // For all relations with other particles
        for (int j = 0; j<NBODY; j++) {

            // Do not include the particle with itself            
            if (i == j) continue;

            // Distance modulus
            double r = sqrt(DIST(i,j,0)*DIST(i,j,0) + 
                            DIST(i,j,1)*DIST(i,j,1) +
                            DIST(i,j,2)*DIST(i,j,2));

            // Perturbated
            VARP(i, 3) += (G*M[j] / r / r / r) * DIST(i,j,0);  // ax
            VARP(i, 4) += (G*M[j] / r / r / r) * DIST(i,j,1);  // ay
            VARP(i, 5) += (G*M[j] / r / r / r) * DIST(i,j,2);  // az

        }

        VARP(i, 0) = VAR(i,3);                                  // vx 
        VARP(i, 1) = VAR(i,4);                                  // vy
        VARP(i, 2) = VAR(i,5);                                  // vz
           
    }

}

































