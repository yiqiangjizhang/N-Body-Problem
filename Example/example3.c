/*
	EXAMPLE 3
	
	Two-Body problem.

	Arnau Miro, 2018
	Last rev: 2020
*/

/* Two-body problem
Date: 18/03/2021
Author/s:
  - Angel Pan Du
  - Alba Molina Cuadrado
  - Iván Sermanoukian Molina
  - Yi Qiang Ji Zhang
Subject: High Performance Computing for Aerospace Engineering
Professor: Manel Soria & Arnau Miro

// Problem statement
------------------------------------------------------------------------
Two Body Problem
-------------------------------------------------------------------------
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string.h>

#include "RK_c.h"

#define PI 4. * atan(1)

// Parameters
#define G   6.674e-11*1.e-9  // km3/kg s2
#define M1  5.97219e24       // kg
#define M2  1.9891e30        // kg
#define R0  1.5e8            // km
#define N   12

// EXAMPLE 2
// #define G   6.674e-11*1.e-9  // km3/kg s2
// #define M1  5.97219e5       // kg
// #define M2  5.97219e5        // kg
// #define R0  2e3            // km
// #define N   12


// Define the function to integrate
void TwoBody(double t, double *var, int n, double *varp)
{
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
    double r = sqrt((var[6] - var[0]) * (var[6] - var[0]) +
                    (var[7] - var[1]) * (var[7] - var[1]) +
                    (var[8] - var[2]) * (var[8] - var[2]));

    // Perturbated
    varp[0] = var[3];
    varp[1] = var[4];
    varp[2] = var[5];
    varp[3] = (G * M2 / r / r / r) * (var[6] - var[0]);
    varp[4] = (G * M2 / r / r / r) * (var[7] - var[1]);
    varp[5] = (G * M2 / r / r / r) * (var[8] - var[2]);

    // Perturber
    varp[6]  = var[9];
    varp[7]  = var[10];
    varp[8]  = var[11];
    varp[9]  = (G * M1 / r / r / r) * (var[0] - var[6]);
    varp[10] = (G * M1 / r / r / r) * (var[1] - var[7]);
    varp[11] = (G * M1 / r / r / r) * (var[2] - var[8]);
}

int main()
{
    double V0 = sqrt(G * M2 / R0);
    double T  = 2*PI/sqrt(G * M2 / R0 / R0 / R0);
    
    // Integration bounds and initial solution
    double tspan[2] = {0., 2. * T};
    double y0[N] = {R0, 0., 0., 0., V0, 0., 0., 0., 0., 0., 0., 0.};
    double var[9];
    double varp[N];


    // Runge-Kutta parameters
    RK_PARAM rkp = rkparams(tspan);
    rkp.h0 = 24.*3600.;
    rkp.eps = 1e-6;

    // Launch the integrator
    RK_OUT rko = odeRK("hiroshi912", TwoBody, tspan, y0, N, &rkp);

    // Finish
    printf("error = %.2e with %d steps\n", rko.err, rko.n);

    writerkout("output.csv", &rko, N);

    freerkout(&rko);

    return 0;
}