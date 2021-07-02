/*
	RUNGE-KUTTA Integration
	Library to perform numerical integration using Runge-Kutta methods, implemented
	based on MATLAB's ode45.
	The function to call is ODERK, a generic Runge-Kutta variable step integrator.
	The inputs of rkm function are:
		> scheme: Runge-Kutta scheme to use. Options are:
			* Hiroshi 9(12) (hiroshi912)
		> odefun: a function so that
			dydx = odefun(x,y,n)
		where n is the number of y variables to integrate.
		> xspan: integration start and end point.
		> y0: initial conditions (must be size n).
		> n: number of initial conditions.
		An additional parameter can be passed to the integrator, which contains
		parameters for the integrator:
			> h0:       Initial step for the interval
			> hmin:      Minimum interval allowed
			> eps:       Tolerance to meet.
			> eventfcn:  Event function.
			> outputfcn: Output function.
	If the inputted tolerance is not met, the integrator will use hmin and will
	produce a successful step. Warning! Accuracy might be compromised.
	The function will return a structure containing the following information:
		> retval: Return value. Will be negative in case of errors or positive if successful.
			* 1: indicates a successful run.
			* 2: indicates that at a certain point hmin was used for the step.
		> n: Number of steps taken.
		> err: Maximum error achieved.
		> x: solution x values of size n.
		> y: solution y values of size n per each variable.
	Arnau Miro, Elena Terzic 2021
	Yi Qiang Ji, Alba Molina, Angel Pan, Iván Sermanoukian 2021
	Last rev: 2021
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "RK.h"

#define NNSTEP 1000
#define AIND(i,j)  rkm.A[rkm.nstep*(i) + (j)]
#define AINDP(i,j) rkm->A[rkm->nstep*(i) + (j)]
#define EQSTR(STR1, STR2) strcmp(STR1, STR2) == 0 


/*
	ODERK
	A generic Runge-Kutta variable step integrator.
	Inputs:
		> scheme: Runge-Kutta scheme to use
		> odefun: must be a function so that dydx = odefun(x,y,n)
		  and n is the number of derivatives to compute.
		> xspan: integration start and end point.
		> y0: initial conditions (must be size n).
		> n: number of initial conditions.
		> rkp: parameters to the Runge-Kutta integrator.
	Outputs:
		> RK_OUT structure containing the exit code and the solution.
*/

RK_OUT odeRK(const char *scheme, void (*odefun)(double, double *, int, double *),
			  double xspan[2], double *y0, const int n, const RK_PARAM *rkp)
{
	RK_OUT rko;

	const double hmin = fabs(xspan[1] - xspan[0])*rkp->minstep;

	// Select Runge-Kutta method
	RKMethod rkm = RKMethodSelection(str2RKS(scheme));

	// Initialization
	int cont = 1, last = 0;
	double h = rkp->h0;
	rko.retval = 0;

	double *f;

	rko.x = (double *)malloc(  (NNSTEP + 1)*sizeof(double));
	rko.y = (double *)malloc(n*(NNSTEP + 1)*sizeof(double));
	f     = (double *)malloc((n*rkm.nstep) *sizeof(double));

	rko.n        = 0; // Start at iteration zero
	rko.err      = 0.;
	rko.x[rko.n] = xspan[0];
	memcpy(rko.y + n*rko.n, y0, n*sizeof(double));

	// Check tableau
	#ifdef RK_CHECK_TABLEAU
	printf("Checking tableau for %s...", scheme);
	if (!CheckTableau(&rkm))
	{
		cont = 0;
		rko.retval = -1;
	}
	#endif

	// Definitions
	int dir[1];
	double ylow[n], yhigh[n], dydx[n], val[1], g_ant = 0.;

	// Runge-Kutta loop
	while (cont)
	{
		// Exit criteria
		if (rko.x[rko.n] + h > xspan[1])
		{
			last = 1;
			// Arrange the step so it finishes at xspan[1]
			h = fabs(rko.x[rko.n] - xspan[1]);
		}

		// Initialise
		memcpy(ylow,  rko.y + n*rko.n, n*sizeof(double));
		memcpy(yhigh, rko.y + n*rko.n, n*sizeof(double));

		// Calculus loop
		for (int ii = 0; ii < rkm.nstep; ii++)
		{
			double xint = rko.x[rko.n] + h*rkm.C[ii];
			double yint[n];

			memcpy(yint, rko.y + n*rko.n, n*sizeof(double));

			for (int kk = 0; kk < n; kk++)
			{
				for (int jj = 0; jj < ii; jj++)
					yint[kk] += AIND(ii,jj)*f[n*jj + kk];
			}

			// Call function
			(*odefun)(xint, yint, n, dydx);

			for (int kk = 0; kk < n; kk++)
			{
				f[n*ii + kk] = h * dydx[kk];
				ylow[kk]    += rkm.Bhat[ii] * f[n*ii + kk];
				yhigh[kk]   += rkm.B[ii]    * f[n*ii + kk];
			}
		}

		// Compute the total error
		// Work with both relative and absolute errors
		double rel_err = 1.e-20, abs_err = 1.e-20; // Avoid division by zero

		for (int kk = 0; kk < n; ++kk) {
			rel_err = fmax(fabs(1. - ylow[kk]/yhigh[kk]),rel_err);
			abs_err = fmax(fabs(yhigh[kk] - ylow[kk])   ,abs_err);
		}
		double error = fmin(rel_err, abs_err);

		// Step size control
		// Source: Ketcheson, David, and Umair bin Waheed.
		//         "A comparison of high-order explicit Runge–Kutta, extrapolation, and deferred correction methods in serial and parallel."
		//         Communications in Applied Mathematics and Computational Science 9.2 (2014): 175-200.
		double hest = rkp->secfact*h*pow(rkp->eps/error, 0.7/rkm.alpha);
		if (rkp->eventfcn)
		{
			// Run event function
			int ccont = rkp->eventfcn(rko.x[rko.n], yhigh, n, val, dir);
			// Naive approximation to a root finding algorithm
			// using a crude mid-point rule
			if (rko.n != 0 && val[0] * g_ant < 0)
			{
				h /= 2.;
				cont = 1;
				continue;
			}
			// Update continue according to the output of the
			// event function
			cont  = (fabs(val[0]) < rkp->epsevf) ? ccont : cont;
			g_ant = (fabs(val[0]) < rkp->epsevf) ?     0 : val[0]; // Also restart g_ant
		}

		if (error < rkp->eps)
		{
			// rkm is a successful step
			rko.retval = (rko.retval <= 0) ? 1 : rko.retval;
			rko.err    = fmax(error, rko.err);
			rko.n++;
			if (last) cont = 0;

			// Reallocate
			if (rko.n % NNSTEP == 0)
			{
				int mult = (int)(rko.n / NNSTEP) + 1;
				rko.x = (double *)realloc(rko.x, mult*(NNSTEP + 1)*sizeof(double)); 
				rko.y = (double *)realloc(rko.y, mult*n*(NNSTEP + 1)*sizeof(double));
			}

			// Set output values
			rko.x[rko.n] = rko.x[rko.n - 1] + h;
			memcpy(rko.y + n*rko.n, yhigh, n*sizeof(double));

			// Output function
			if (rkp->outputfcn)
				cont = (cont == 1) ? rkp->outputfcn(rko.x[rko.n - 1], yhigh, n) : 0;

			// Set the new step
			// Source: Ketcheson, David, and Umair bin Waheed.
			//         "A comparison of high-order explicit Runge–Kutta, extrapolation, and deferred correction methods in serial and parallel."
			//         Communications in Applied Mathematics and Computational Science 9.2 (2014): 175-200.
			h = fmin(rkp->secfact_max * h, fmax(rkp->secfact_min * h, hest));
			//h = fmax(h,hest);
		}
		else
		{
			// rkm is a failed step
			if (h < hmin)
			{
				// Check if our step has decreased too much
				cont = 0;
				rko.retval = -1;
			}
			else
			{
				h = hest;
			}
		}
	}

	// Allocate memory to output structure and copy the data
	rko.x  = (double *)realloc(rko.x,(rko.n + 1)   * sizeof(double));
	rko.y  = (double *)realloc(rko.y,n*(rko.n + 1) * sizeof(double));
    rko.dy = NULL;

	// Deallocate memory from the output structure
	free(f);
    RKMethodFree(&rkm);

	return rko;
}

/*
	ODERKN
	A generic Runge-Kutta-Nystrom variable step integrator.
	Inputs:
		> scheme: Runge-Kutta scheme to use
		> odefun: must be a function so that dydx = odefun(x,y,n)
		  and n is the number of derivatives to compute.
		> xspan: integration start and end point.
		> y0: initial conditions (must be size n).
		> dy0: initial conditions for the derivative (must be size n).
		> n: number of initial conditions.
		> rkp: parameters to the Runge-Kutta integrator.
	Outputs:
		> RK_OUT structure containing the exit code and the solution.
*/

RK_OUT odeRKN(const char *scheme, void (*odefun)(double, double *, int, double *),
			   double xspan[2], double *y0, double *dy0, const int n, const RK_PARAM *rkp)
{

	RK_OUT rko;

	const double hmin = fabs(xspan[1] - xspan[0]) * rkp->minstep;
	// Select Runge-Kutta method
	RKMethod rkm = RKMethodSelection(str2RKS(scheme));

	// Initialization
	int cont = 1, last = 0;
	double h = rkp->h0;
	rko.retval = 0;

	double *f, *f2;

	rko.x  = (double *)malloc((NNSTEP + 1)   * sizeof(double));
	rko.y  = (double *)malloc(n*(NNSTEP + 1) * sizeof(double));
	rko.dy = (double *)malloc(n*(NNSTEP + 1) * sizeof(double));
	f      = (double *)malloc((n*rkm.nstep)  * sizeof(double));
	f2     = (double *)malloc((n*rkm.nstep)  * sizeof(double));
	
	rko.n        = 0; // Start at iteration zero
	rko.err      = 0.;
	rko.x[rko.n] = xspan[0];
	memcpy(rko.y  + n*rko.n,  y0, n*sizeof(double));
	memcpy(rko.dy + n*rko.n, dy0, n*sizeof(double));

	// Definitions
	int dir[1];
	double ylow[n], yhigh[n], dylow[n], dyhigh[n], dy2dx[n], val[1], g_ant = 0.;
	
	// Runge-Kutta loop
	while (cont)
	{
		// Exit criteria
		if (rko.x[rko.n] + h > xspan[1])
		{
			last = 1;
			// Arrange the step so it finishes at xspan[1]
			h = fabs(rko.x[rko.n] - xspan[1]);
		}

		// Initialise
		memcpy(ylow,   rko.y  + n*rko.n, n*sizeof(double));
		memcpy(yhigh,  rko.y  + n*rko.n, n*sizeof(double));
		memcpy(dylow,  rko.dy + n*rko.n, n*sizeof(double));
		memcpy(dyhigh, rko.dy + n*rko.n, n*sizeof(double));

		for (int kk = 0; kk < n; kk++)
		{
			ylow[kk]  += h*dylow[kk];
			yhigh[kk] += h*dyhigh[kk];
		}

		// Calculus loop
		for (int ii = 0; ii < rkm.nstep; ii++)
		{
			double xint = rko.x[rko.n] + h*rkm.C[ii];
			double yint[n], dyint[n];

			memcpy(yint,  rko.y  + n*rko.n, n*sizeof(double));
			memcpy(dyint, rko.dy + n*rko.n, n*sizeof(double));
			
			for (int kk = 0; kk < n; kk++)
			{
				yint[kk] += h*rkm.C[ii]*dyint[kk];
				for (int jj = 0; jj < ii; jj++)
					yint[kk] += AIND(ii,jj)*f2[n*jj + kk];
			}

			// Call function
			(*odefun)(xint, yint, n, dy2dx);
			
			for (int kk = 0; kk < n; kk++)
			{
				f[n*ii + kk] = h*dy2dx[kk];
				f2[n*ii + kk] = h*h*dy2dx[kk];
				dylow[kk]  += rkm.Bphat[ii] * f[n * ii + kk];
				dyhigh[kk] += rkm.Bp[ii]    * f[n * ii + kk];
				ylow[kk]   += rkm.Bhat[ii]  * f2[n * ii + kk];
				yhigh[kk]  += rkm.B[ii]     * f2[n * ii + kk];
			}
		}
		// Compute the total error
		// Work with both relative and absolute errors
		double error = 1.e-20, rel_err[2] = {1.e-20, 1.e-20}, abs_err[2] = {1.e-20, 1.e-20};
		for (int kk = 0; kk < n; kk++)
		{
			rel_err[0] = fmax(fabs(1. - dylow[kk]/dyhigh[kk]), rel_err[0]);
			rel_err[1] = fmax(fabs(1. - ylow[kk]/yhigh[kk]),   rel_err[1]);
			abs_err[0] = fmax(fabs(dyhigh[kk] - dylow[kk]),    abs_err[0]);
			abs_err[1] = fmax(fabs(yhigh[kk] - ylow[kk]),      abs_err[1]);
		}
		error = fmin(fmax(rel_err[0], rel_err[1]), fmax(abs_err[0], abs_err[1]));
		// Step size control
		// Source: Ketcheson, David, and Umair bin Waheed.
		//         "A comparison of high-order explicit Runge–Kutta, extrapolation, and deferred correction methods in serial and parallel."
		//         Communications in Applied Mathematics and Computational Science 9.2 (2014): 175-200.
		double hest = rkp->secfact*h*pow(rkp->eps / error, 0.7 / rkm.alpha);
		if (rkp->eventfcn)
		{
			// Run event function
			int ccont = rkp->eventfcn(rko.x[rko.n], yhigh, n, val, dir);
			// Naive approximation to a root finding algorithm
			// using a crude mid-point rule
			if (rko.n != 0 && val[0] * g_ant < 0)
			{
				h /= 2.;
				cont = 1;
				continue;
			}
			// Update continue according to the output of the
			// event function
			cont = (fabs(val[0])  < rkp->epsevf) ? ccont : cont;
			g_ant = (fabs(val[0]) < rkp->epsevf) ? 0     : val[0]; // Also restart g_ant
		}
		if (error < rkp->eps)
		{
			// rkm is a successful step
			rko.retval = (rko.retval <= 0) ?     1 : rko.retval;
			rko.err    = fmax(error,rko.err);
			rko.n++;
			if (last) cont = 0;

			// Reallocate
			if (rko.n % NNSTEP == 0)
			{
				int mult = (int)(rko.n / NNSTEP) + 1;
				rko.x =  (double *)realloc(rko.x,  mult*(NNSTEP + 1)*sizeof(double));
				rko.y =  (double *)realloc(rko.y,  mult*n*(NNSTEP + 1)*sizeof(double));
				rko.dy = (double *)realloc(rko.dy, mult*n*(NNSTEP + 1)*sizeof(double));
			}
			// Set output values
			rko.x[rko.n] = rko.x[rko.n - 1] + h;
			memcpy(rko.y  + n*rko.n,  yhigh, n*sizeof(double));
			memcpy(rko.dy + n*rko.n, dyhigh, n*sizeof(double));
			// Output function
			if (rkp->outputfcn)
				cont = (cont == 1) ? rkp->outputfcn(rko.x[rko.n - 1], yhigh, n) : 0;
			// Set the new step
			// Source: Ketcheson, David, and Umair bin Waheed.
			//         "A comparison of high-order explicit Runge–Kutta, extrapolation, and deferred correction methods in serial and parallel."
			//         Communications in Applied Mathematics and Computational Science 9.2 (2014): 175-200.
			h = fmin(rkp->secfact_max * h, fmax(rkp->secfact_min * h, hest));
			//h = fmax(h,hest);
		}
		else
		{
			// rkm is a failed step
			if (h < hmin)
			{
				// Check if our step has decreased too much
				cont = 0;
				rko.retval = -1;
			}
			else
			{
				h = hest;
			}
		}
	}

	// Allocate memory to output structure and copy the data
	rko.x  = (double *)realloc(rko.x ,  (rko.n + 1)*sizeof(double));
	rko.y  = (double *)realloc(rko.y ,n*(rko.n + 1)*sizeof(double));
	rko.dy = (double *)realloc(rko.dy,n*(rko.n + 1)*sizeof(double));

	// Free memory
	free(f); free(f2);
	RKMethodFree(&rkm);

	return rko;
}

/*
	RUNGE-KUTTA OUTPUT
	Output for the Runge-Kutta integrator.
*/

void freerkout(const RK_OUT *rko)
{
    if (rko->x  != NULL) free(rko->x);
	if (rko->y  != NULL) free(rko->y);
	if (rko->dy != NULL) free(rko->dy);
}

void writerkout(const char *outfile, const RK_OUT *rko, const int nvars)
{
	FILE *myfile;
	myfile = fopen(outfile, "w");

	for (int ii = 0; ii < rko->n; ++ii)
	{
		fprintf(myfile, "%f\t", rko->x[ii]); // time variable (first column)
		for (int jj = 0; jj < nvars; ++jj)
			fprintf(myfile, "%f\t", rko->y[nvars * ii + jj]);
		fprintf(myfile, "\n");
	}
	fclose(myfile);
}

const char *RKS2str(const RK_SCHEME rks)
{
	switch (rks) {
		case NONE:				return "none";
		case EULERHEUN12:       return "eulerheun12";
		case BOGACKISHAMPINE23: return "bogackishampine23";
		case DORMANDPRINCE34A:  return "dormandprince34a";
		case FEHLBERG45:        return "fehlberg45";
		case CASHKARP45:        return "cashkarp45";
		case DORMANDPRINCE45:   return "dormandprince45";
		case DORMANDPRINCE45A:  return "dormandprince45a";
		case CALVO56:           return "calvo56";
		case DORMANDPRINCE78:   return "dormandprince78";
		case CURTIS810:         return "curtis810";
		case HIROSHI912:        return "hiroshi912";
		case RKN34:             return "rnk34";
		case RKN46:             return "rnk46";
		case RKN68:             return "rnk68";
		case RKN1012:           return "rnk1012";
	}
	return "notfound";
}

RK_SCHEME str2RKS(const char *str)
{
	if (EQSTR(str,"eulerheun12"))       return EULERHEUN12;
	if (EQSTR(str,"bogackishampine23")) return BOGACKISHAMPINE23;
	if (EQSTR(str,"dormandprince34a"))  return DORMANDPRINCE34A;
	if (EQSTR(str,"fehlberg45"))        return FEHLBERG45;
	if (EQSTR(str,"cashkarp45"))        return CASHKARP45;
	if (EQSTR(str,"dormandprince45"))   return DORMANDPRINCE45;
	if (EQSTR(str,"dormandprince45a"))  return DORMANDPRINCE45A;
	if (EQSTR(str,"calvo56"))           return CALVO56;
	if (EQSTR(str,"dormandprince78"))   return DORMANDPRINCE78;
	if (EQSTR(str,"curtis810"))         return CURTIS810;
	if (EQSTR(str,"hiroshi912"))        return HIROSHI912;
	if (EQSTR(str,"rkn34"))             return RKN34;
	if (EQSTR(str,"rkn46"))             return RKN46;
	if (EQSTR(str,"rkn68"))             return RKN68;
	if (EQSTR(str,"rkn1012"))           return RKN1012;
    return NONE;
}

RK_PARAM rkparams(const double xspan[2])
{
    RK_PARAM rkp;
    // Initial and max step based on xspan
    rkp.h0          = (xspan[1] - xspan[0]) / 10.;
    // Tolerance
    rkp.eps         = 1e-8;
    rkp.epsevf      = 1e-4;
    rkp.minstep     = 1.e-12;
    // Security factors
    rkp.secfact     = 0.9;
    rkp.secfact_max = 5.;
    rkp.secfact_min = 0.2;
    // No event or output function
    rkp.eventfcn    = NULL;
    rkp.outputfcn   = NULL;
    return rkp;
}

/*
	CHECK TABLEAU
	Checks for errors in the Butcher's tableau
*/
int check_tableau(const char *scheme)
{
	RKMethod rkm = RKMethodSelection(str2RKS(scheme));
	int retval   = CheckTableau(&rkm);
	RKMethodFree(&rkm);
	return retval;
}

/*
	RUNGE KUTTA METHOD
	Class constructor and destructor
*/
RKMethod RKMethodSelection(const RK_SCHEME rks) 
{
	RKMethod rkm;
	rkm.alloc = 0;

	// Selection of the Runge-Kutta method according
	// to the input scheme name
	switch (rks) {
		case EULERHEUN12:       
			EulerHeun12(&rkm);
			break;
		case BOGACKISHAMPINE23: 
			BogackiShampine23(&rkm);
			break;
		case DORMANDPRINCE34A:  
			DormandPrince34A(&rkm);
			break;
		case FEHLBERG45:        
			Fehlberg45(&rkm);
			break;
		case CASHKARP45:        
			CashKarp45(&rkm);
			break;
		case DORMANDPRINCE45:   
			DormandPrince45(&rkm);
			break;
		case DORMANDPRINCE45A:  
			DormandPrince45A(&rkm);
			break;
		case CALVO56:           
			Calvo56(&rkm);
			break;
		case DORMANDPRINCE78:   
			DormandPrince78(&rkm);
			break;
		case CURTIS810:         
			Curtis810(&rkm);
			break;
		case HIROSHI912:        
			Hiroshi912(&rkm);
			break;
		case RKN34:             
			RungeKuttaNystrom34(&rkm);
			break;
		case RKN46:             
			RungeKuttaNystrom46(&rkm);
			break;
		case RKN68:             
			RungeKuttaNystrom68(&rkm);
			break;
		case RKN1012:           
			RungeKuttaNystrom1012(&rkm);
			break;
		default:
			printf("Warning! Scheme <%s> not found\n",RKS2str(rks));
	}
	return rkm;
}

void RKMethodFree(RKMethod *rkm)
{
	if (rkm->alloc)
	{
		free(rkm->C);
		free(rkm->A);
		free(rkm->B);
		free(rkm->Bhat);
		if (rkm->Bp != NULL)    free(rkm->Bp);
		if (rkm->Bphat != NULL) free(rkm->Bphat);
	}
}

/*
	EULERHEUN12

	Runge-Kutta Euler-Heun 1-2 integrator.
	Source: Wikipedia
*/
void EulerHeun12(RKMethod *rkm) {

	rkm->nstep = 2;
	rkm->alpha = 1.;

	// Allocate Butcher's tableau
	rkm->C     = (double *)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double *)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double *)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double *)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = 1.;

	// 2nd order solution
	rkm->B[0] = .5;
	rkm->B[1] = .5;

	// 1st order solution
	rkm->Bhat[0] = 1.;
	rkm->Bhat[1] = 0.;
}

/*
	BOGACKISHAMPINE23

	Runge-Kutta Bogacki-Shampine 2-3 integrator.
	Source: Wikipedia
*/
void BogackiShampine23(RKMethod *rkm) {

	rkm->nstep = 4;
	rkm->alpha = 2.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = .25;
	rkm->C[2] = 3./4.;
	rkm->C[3] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = .25;
	AINDP(2,0) = 0.;
	AINDP(2,1) = 3./4.;
	AINDP(3,0) = 2./9.;
	AINDP(3,1) = 1./3.;
	AINDP(3,2) = 4./9.;

	// 3rd order solution
	rkm->B[0] = 2./9.;
	rkm->B[1] = 1./3.;
	rkm->B[2] = 4./9.;
	rkm->B[3] = 0.;

	// 2nd order solution
	rkm->Bhat[0] = 7./24.;
	rkm->Bhat[1] = 1./4.;
	rkm->Bhat[2] = 1./3.;
	rkm->Bhat[3] = 1./8.;
}

/*
	DORMANDPRINCE34A

	Runge-Kutta Dormand-Prince 3-4 for Astrodynamics.
	Source: “New Runge-Kutta algorithms for numerical simulation in dynamical astronomy”
			J. R. Dormand and P. J. Prince, 
			Celest. Mech., vol. 18, no. 3, pp. 223–232, Oct. 1978. 
*/
void DormandPrince34A(RKMethod *rkm) {

	const double lambda = 0.1;

	rkm->nstep = 5;
	rkm->alpha = 3.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = 0.5;
	rkm->C[2] = 0.5;
	rkm->C[3] = 1.;
	rkm->C[4] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = 0.5;
	AINDP(2,0) = 0.;
	AINDP(2,1) = 0.5;
	AINDP(3,0) = 0.;
	AINDP(3,1) = 0.;
	AINDP(3,2) = 1.;
	AINDP(4,0) = 1./6.;
	AINDP(4,1) = 1./3.;
	AINDP(4,2) = 1./3.;
	AINDP(4,3) = 1./6.;

	// 5th order solution
	rkm->B[0] = 1./6.;
	rkm->B[1] = 1./3.;
	rkm->B[2] = 1./3.;
	rkm->B[3] = 1./6.;
	rkm->B[4] = 0.;

	// 4th order solution
	rkm->Bhat[0] = 1./6.;
	rkm->Bhat[1] = 1./3.;
	rkm->Bhat[2] = 1./3.;
	rkm->Bhat[3] = 1./6.-lambda;
	rkm->Bhat[4] = lambda;
}

/*
	FEHLBERG45

	Runge-Kutta Fehlberg 4-5.
	Source: Wikipedia
*/
void Fehlberg45(RKMethod *rkm) {

	rkm->nstep = 6;
	rkm->alpha = 4.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = .25;
	rkm->C[2] = 3./8.;
	rkm->C[3] = 12./13.;
	rkm->C[4] = 1.;
	rkm->C[5] = 1./2.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = .25;
	AINDP(2,0) = 3./32.;
	AINDP(2,1) = 9./32.;
	AINDP(3,0) = 1932./2197.;
	AINDP(3,1) = -7200./2197.;
	AINDP(3,2) = 7296./ 2197.;
	AINDP(4,0) = 439./216.;
	AINDP(4,1) = -8.;
	AINDP(4,2) = 3680./513.;
	AINDP(4,3) = -845./4104.;
	AINDP(5,0) = -8./27.;
	AINDP(5,1) = 2.;
	AINDP(5,2) = -3544./2565.;
	AINDP(5,3) = 1859./4104.;
	AINDP(5,4) = -11./40.;

	// 5th order solution
	rkm->B[0] = 16./135.;
	rkm->B[1] = 0.;
	rkm->B[2] = 6656./12825.;
	rkm->B[3] = 28561./56430;
	rkm->B[4] = -9./50.;
	rkm->B[5] = 2./55.;

	// 4th order solution
	rkm->Bhat[0] = 25./216.;
	rkm->Bhat[1] = 0.;
	rkm->Bhat[2] = 1408./2565.; 
	rkm->Bhat[3] = 2197./4104.;
	rkm->Bhat[4] = -1./5.;
	rkm->Bhat[5] = 0.;
}

/*
	CASHKARP45

	Runge-Kutta Cash-Karp 4-5.
	Source: Wikipedia
*/
void CashKarp45(RKMethod *rkm) {

	rkm->nstep = 6;
	rkm->alpha = 4.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = 1./5.;
	rkm->C[2] = 3./10.;
	rkm->C[3] = 3./5.;
	rkm->C[4] = 1.;
	rkm->C[5] = 7./8.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = 1./5.;
	AINDP(2,0) = 3./40.;
	AINDP(2,1) = 9./40.;
	AINDP(3,0) = 3./10.;
	AINDP(3,1) = -9./10.;
	AINDP(3,2) = 6./5.;
	AINDP(4,0) = -11./54.;
	AINDP(4,1) = 5./2.;
	AINDP(4,2) = -70./27.;
	AINDP(4,3) = 35./27.;
	AINDP(5,0) = 1631./55296.;
	AINDP(5,1) = 175./512.; 
	AINDP(5,2) = 575./13824.;
	AINDP(5,3) = 44275./110592.;
	AINDP(5,4) = 253./4096.;

	// 5th order solution
	rkm->B[0] = 37./378.;
	rkm->B[1] = 0.;
	rkm->B[2] = 250./621.;
	rkm->B[3] = 125./594.;
	rkm->B[4] = 0.;
	rkm->B[5] = 512./1771.;

	// 4th order solution
	rkm->Bhat[0] = 2825./27648.;
	rkm->Bhat[1] = 0.;
	rkm->Bhat[2] = 18575./48384.;
	rkm->Bhat[3] = 13525./55296.;
	rkm->Bhat[4] = 277./14336.;
	rkm->Bhat[5] = 1./4.;
}

/*
	DORMANDPRINCE45

	Runge-Kutta Dormand-Prince 4-5.
	Source: Wikipedia
*/
void DormandPrince45(RKMethod *rkm) {

	rkm->nstep = 7;
	rkm->alpha = 4.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = .2;
	rkm->C[2] = .3;
	rkm->C[3] = 4./5.;
	rkm->C[4] = 8./9.;
	rkm->C[5] = 1.;
	rkm->C[6] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = .2;
	AINDP(2,0) = .3/4.;
	AINDP(2,1) = .9/4.;
	AINDP(3,0) = 44./45.;
	AINDP(3,1) = -56./15.;
	AINDP(3,2) = 32./ 9.;
	AINDP(4,0) = 19372./6561.;
	AINDP(4,1) = -25360./2187.;
	AINDP(4,2) = 64448./6561.;
	AINDP(4,3) = -212./729.;
	AINDP(5,0) = 9017./3168.;
	AINDP(5,1) = -355./33.;
	AINDP(5,2) = 46732./5247.;
	AINDP(5,3) = 49./176.;
	AINDP(5,4) = -5103./18656.;
	AINDP(6,0) = 35./384.;
	AINDP(6,1) = 0.;
	AINDP(6,2) = 500./1113.;
	AINDP(6,3) = 125./192.;
	AINDP(6,4) = -2187./6784.;
	AINDP(6,5) = 11./84.;

	// 5th order solution
	rkm->B[0] = 35./384.;
	rkm->B[1] = 0.;
	rkm->B[2] = 500./1113.;
	rkm->B[3] = 125./192.;
	rkm->B[4] = -2187./6784.;
	rkm->B[5] = 11./84.;
	rkm->B[6] = 0.;

	// 4th order solution
	rkm->Bhat[0] = 5179./57600.;
	rkm->Bhat[1] = 0.;
	rkm->Bhat[2] = 7571./16695.;
	rkm->Bhat[3] = 393./640.;
	rkm->Bhat[4] = -92097./339200.;
	rkm->Bhat[5] = 187./2100.;
	rkm->Bhat[6] = 1./40.;
}

/*
	DORMANDPRINCE45A

	Runge-Kutta Dormand-Prince 4-5 for Astrodynamics.
	Source: “New Runge-Kutta algorithms for numerical simulation in dynamical astronomy”
			J. R. Dormand and P. J. Prince, 
			Celest. Mech., vol. 18, no. 3, pp. 223–232, Oct. 1978.
*/
void DormandPrince45A(RKMethod *rkm) {

	const double lambda = 1./60.;

	rkm->nstep = 7;
	rkm->alpha = 4.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = .125;
	rkm->C[2] = .25;
	rkm->C[3] = 4./9.;
	rkm->C[4] = 4./5.;
	rkm->C[5] = 1.;
	rkm->C[6] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = .125;
	AINDP(2,0) = 0.;
	AINDP(2,1) = 1./4.;
	AINDP(3,0) = 196./729.;
	AINDP(3,1) = -320./729.;
	AINDP(3,2) = 448./729.;
	AINDP(4,0) = 836./2875.;
	AINDP(4,1) = 64./575.;
	AINDP(4,2) = -13376./20125.;
	AINDP(4,3) = 21384./20125.;
	AINDP(5,0) = -73./48.;
	AINDP(5,1) = 0.;
	AINDP(5,2) = 1312./231.;
	AINDP(5,3) = -2025./448.;
	AINDP(5,4) = 2875./2112.;
	AINDP(6,0) = 17./192.;
	AINDP(6,1) = 0.;
	AINDP(6,2) = 64./231.;
	AINDP(6,3) = 2187./8960.;
	AINDP(6,4) = 2875./8448.;
	AINDP(6,5) = 1./20.;

	// 5th order solution
	rkm->B[0] = 17./192.;
	rkm->B[1] = 0.;
	rkm->B[2] = 64./231.;
	rkm->B[3] = 2187./8960.;
	rkm->B[4] = 2875./8448.;
	rkm->B[5] = 1./20.;
	rkm->B[6] = 0.;

	// 4th order solution
	rkm->Bhat[0] = 17./192.;
	rkm->Bhat[1] = 0.;
	rkm->Bhat[2] = 64./231.;
	rkm->Bhat[3] = 2187./8960.;
	rkm->Bhat[4] = 2875./8448.;
	rkm->Bhat[5] = 1./20.-lambda;
	rkm->Bhat[6] = lambda;
}

/*
	CALVO56

	Runge-Kutta Calvo 5-6.
	Source: "A new embedded pair of Runge-Kutta formulas of orders 5 and 6"
	M. Calvo, J. I. Montijano, and L. Randez, Comput. Math. with Appl., vol. 20, no. 1, pp. 15–24, 1990.
*/
void Calvo56(RKMethod *rkm) {

	rkm->nstep = 9;
	rkm->alpha = 5.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = 2./15.;
	rkm->C[2] = 1./5.;
	rkm->C[3] = 3./10.;
	rkm->C[4] = 14./25.;
	rkm->C[5] = 19./25.;
	rkm->C[6] = 35226607./35688279.;
	rkm->C[7] = 1.;
	rkm->C[8] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = 2./15.;
	AINDP(2,0) = 1./20.;
	AINDP(2,1) = 3./20.;
	AINDP(3,0) = 3./40.;
	AINDP(3,1) = 0.;
	AINDP(3,2) = 9./40.;
	AINDP(4,0) = 86727015./196851553.;
	AINDP(4,1) = -60129073./52624712.;
	AINDP(4,2) = 957436434./1378352377.;
	AINDP(4,3) = 83886832./147842441.;
	AINDP(5,0) = -86860849./45628967.;
	AINDP(5,1) = 111022885./25716487.;
	AINDP(5,2) = 108046682./101167669.;
	AINDP(5,3) = -141756746./36005461.;
	AINDP(5,4) = 73139862./60170633.;
	AINDP(6,0) = 77759591./16096467.;
	AINDP(6,1) = -49252809./6452555.;
	AINDP(6,2) = -381680111./51572984.;
	AINDP(6,3) = 879269579./66788831.;
	AINDP(6,4) = -90453121./33722162.;
	AINDP(6,5) = 111179552./157155827.;
	AINDP(7,0) = 237564263./39280295.;
	AINDP(7,1) = -100523239./10677940.;
	AINDP(7,2) = -265574846./27330247.;
	AINDP(7,3) = 317978411./18988713.;
	AINDP(7,4) = -124494385./35453627.;
	AINDP(7,5) = 86822444./100138635.;
	AINDP(7,6) = -12873523./724232625.;
	AINDP(8,0) = 17572349./289262523.;
	AINDP(8,1) = 0.;
	AINDP(8,2) = 57513011./201864250.;
	AINDP(8,3) = 15587306./354501571.;
	AINDP(8,4) = 71783021./234982865.;
	AINDP(8,5) = 29672000./180480167.;
	AINDP(8,6) = 65567621./127060952.;
	AINDP(8,7) = -79074570./210557597.;

	// 8th order solution
	rkm->B[0] = 17572349./289262523.;
	rkm->B[1] = 0.;
	rkm->B[2] = 57513011./201864250.;
	rkm->B[3] = 15587306./354501571.;
	rkm->B[4] = 71783021./234982865.;
	rkm->B[5] = 29672000./180480167.;
	rkm->B[6] = 65567621./127060952.;
	rkm->B[7] = -79074570./210557597.;
	rkm->B[8] = 0.;

	// 5th order solution
	rkm->Bhat[0] = 15231665./510830334.;
	rkm->Bhat[1] = 0.;
	rkm->Bhat[2] = 59452991./116050448.;
	rkm->Bhat[3] = -28398517./122437738.;
	rkm->Bhat[4] = 56673824./137010559.;
	rkm->Bhat[5] = 68003849./426673583.;
	rkm->Bhat[6] = 7097631./37564021.;
	rkm->Bhat[7] = -71226429./583093742.;
	rkm->Bhat[8] = 1./20.;
}

/*
	DORMANDPRINCE78

	Runge-Kutta Dormand-Prince 7-8.
	Source: "High Order Embedded Runge-Kutta Formulae"
	P.J. Prince and J.R. Dormand, J. Comp. Appl. Math.,7, pp. 67-75, 1981
*/
void DormandPrince78(RKMethod *rkm) {

	rkm->nstep = 13;
	rkm->alpha = 7.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0]  = 0.;
	rkm->C[1]  = .5555555555555555555555555555555555555555555555555555555555555555555555555555555555556e-1;
	rkm->C[2]  = .8333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1;
	rkm->C[3]  = .125;
	rkm->C[4]  = .3125;
	rkm->C[5]  = .375;
	rkm->C[6]  = .1475;
	rkm->C[7]  = .465;
	rkm->C[8]  = .5648654513822595753983585014261682587385670087264133215107858541861772043626893634869;
	rkm->C[9]  = .65;
	rkm->C[10] = .9246562776405044467450135743183695426492034467027398177693057974193090932038137681734;
	rkm->C[11] = 1.;
	rkm->C[12] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0)   = .5555555555555555555555555555555555555555555555555555555555555555555555555555555555556e-1;
	AINDP(2,0)   = .2083333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1;
	AINDP(2,1)   = .625e-1;
	AINDP(3,0)   = .3125e-1;
	AINDP(3,1)   = 0.;
	AINDP(3,2)   = .9375e-1;
	AINDP(4,0)   = .3125;
	AINDP(4,1)   = 0.;
	AINDP(4,2)   = -1.171875;
	AINDP(4,3)   = 1.171875;
	AINDP(5,0)   = .375e-1;
	AINDP(5,1)   = 0.;
	AINDP(5,2)   = 0.;
	AINDP(5,3)   = .1875;
	AINDP(5,4)   = .15;
	AINDP(6,0)   = .4791013711111111111111111111111111111111111111111111111111111111111111111111111111111e-1;
	AINDP(6,1)   = 0.;
	AINDP(6,2)   = 0.;
	AINDP(6,3)   = .1122487127777777777777777777777777777777777777777777777777777777777777777777777777778;
	AINDP(6,4)   = -.2550567377777777777777777777777777777777777777777777777777777777777777777777777777778e-1;
	AINDP(6,5)   = .1284682388888888888888888888888888888888888888888888888888888888888888888888888888889e-1;
	AINDP(7,0)   = .1691798978729228118143110713603823606551492879543441068765183916901985069878116840258e-1;
	AINDP(7,1)   = 0.;
	AINDP(7,2)   = 0.;
	AINDP(7,3)   = .3878482784860431695265457441593733533707275558526027137301504136188413357699752956243;
	AINDP(7,4)   = .3597736985150032789670088963477236800815873945874968299566918115320174598617642813621e-1;
	AINDP(7,5)   = .1969702142156660601567152560721498881281698021684239074332818272245101975612112344987;
	AINDP(7,6)   = -.1727138523405018387613929970023338455725710262752107148467532611655731300161441266618;
	AINDP(8,0)   = .6909575335919230064856454898454767856104137244685570409618590205329379141801779261271e-1;
	AINDP(8,1)   = 0.;
	AINDP(8,2)   = 0.;
	AINDP(8,3)   = -.6342479767288541518828078749717385465426336031120747443029278238678259156558727343870;
	AINDP(8,4)   = -.1611975752246040803668769239818171234422237864808493434053355718725564004275765826294;
	AINDP(8,5)   = .1386503094588252554198669501330158019276654889495806914244800957316809692178564181463;
	AINDP(8,6)   = .9409286140357562697242396841302583434711811483145028697758469500622977999086075976730;
	AINDP(8,7)   = .2116363264819439818553721171319021047635363886083981439225363020792869599016568720713;
	AINDP(9,0)   = .1835569968390453854898060235368803084497642516277329033567822939550366893470341427061;
	AINDP(9,1)   = 0.;
	AINDP(9,2)   = 0.;
	AINDP(9,3)   = -2.468768084315592452744315759974107457777649753547732495587510208118948846864075671103;
	AINDP(9,4)   = -.2912868878163004563880025728039519800543380294081727678722180378521228988619909525814;
	AINDP(9,5)   = -.2647302023311737568843979946594614325963055282676866851254669985184857900430792761417e-1;
	AINDP(9,6)   = 2.847838764192800449164518254216773770231582185011953158945518598604677368771006006947;
	AINDP(9,7)   = .2813873314698497925394036418267117820980705455360140173438806414700945017562775967023;
	AINDP(9,8)   = .1237448998633146576270302126636397203122013536069738523260934117931117648560568049432;
	AINDP(10,0)  = -1.215424817395888059160510525029662994880024988044112261492933435562347094075572196306;
	AINDP(10,1)  = 0.;
	AINDP(10,2)  = 0.;
	AINDP(10,3)  = 16.67260866594577243228041328856410774858907460078758981907659463146817850515953150318;
	AINDP(10,4)  = .9157418284168179605957186504507426331593736334298114556675054711691705693180672015549;
	AINDP(10,5)  = -6.056605804357470947554505543091634004081967083324413691952297768849489422976536037188;
	AINDP(10,6)  = -16.00357359415617811184170641007882303068079304063907122528682690677051264091851070681;
	AINDP(10,7)  = 14.84930308629766255754539189802663208272299893302745636963476380528935538206408294404;
	AINDP(10,8)  = -13.37157573528984931829304139618159579089195289469740214314819172008509426917249413996;
	AINDP(10,9)  = 5.134182648179637933173253611658602898712494286162881495270691720760048063805245199662;
	AINDP(11,0)  = .2588609164382642838157309322317577667296307766301063163257807024364127266506000943372;
	AINDP(11,1)  = 0.;
	AINDP(11,2)  = 0.;
	AINDP(11,3)  = -4.774485785489205112310117509706042746829391853746564335432052648850273207666153414325;
	AINDP(11,4)  = -.4350930137770325094407004118103177819323551661617974116300885410959889587870529265443;
	AINDP(11,5)  = -3.049483332072241509560512866312031613982854911220735245095857885009959537455491093975;
	AINDP(11,6)  = 5.577920039936099117423676634464941858623588944531498679305832747426169795721583100328;
	AINDP(11,7)  = 6.155831589861040097338689126688954481197754937462937466615262374558053507207577588692;
	AINDP(11,8)  = -5.062104586736938370077406433910391644990220712141673880668482097753119682588189892664;
	AINDP(11,9)  = 2.193926173180679061274914290465806019788262707389033759251111926114017568440058142981;
	AINDP(11,10) = .1346279986593349415357262378873236613955852772571946513284934221746877884770684011715;
	AINDP(12,0)  = .8224275996265074779631682047726665909572303617765850630130165537026064642659285212794;
	AINDP(12,1)  = 0.;
	AINDP(12,2)  = 0.;
	AINDP(12,3)  = -11.65867325727766428397655303545841477547369082638864247596731917545919361630644756876;
	AINDP(12,4)  = -.7576221166909361958811161540882449653663757591954118634754438929149012424414834787423;
	AINDP(12,5)  = .7139735881595815279782692827650546753142248878566301594716125352788068216039494192152;
	AINDP(12,6)  = 12.07577498689005673956617044860067967095705800972197897084832370633043082304810801069;
	AINDP(12,7)  = -2.127659113920402656390820858969398635427927973275119750336406473203376323696576615278;
	AINDP(12,8)  = 1.990166207048955418328071698344314152176173017359794982793519250552799420955298804232;
	AINDP(12,9)  = -.2342864715440402926602946918568015314512425817004394700255034276765120670064339827205;
	AINDP(12,10) = .1758985777079422650731051058901448183145508638446243836782009233893397195776568900819;
	AINDP(12,11) = 0.;

	// 8th order solution
	rkm->B[0]  = 4.17474911415302462220859284685e-2;
	rkm->B[1]  = 0.;
	rkm->B[2]  = 0.;
	rkm->B[3]  = 0.;
	rkm->B[4]  = 0.;
	rkm->B[5]  = -5.54523286112393089615218946547e-2;
	rkm->B[6]  = 2.39312807201180097046747354249e-1;
	rkm->B[7]  = 7.0351066940344302305804641089e-1;
	rkm->B[8]  = -7.59759613814460929884487677085e-1;
	rkm->B[9]  = 6.60563030922286341461378594838e-1;
	rkm->B[10] = 1.58187482510123335529614838601e-1;
	rkm->B[11] = -2.38109538752862804471863555306e-1;
	rkm->B[12] = 2.5e-1;

	// 7th order solution
	rkm->Bhat[0]  = 2.9553213676353496981964883112e-2;
	rkm->Bhat[1]  = 0.;
	rkm->Bhat[2]  = 0.;
	rkm->Bhat[3]  = 0.;
	rkm->Bhat[4]  = 0.;
	rkm->Bhat[5]  = -8.28606276487797039766805612689e-1;
	rkm->Bhat[6]  = 3.11240900051118327929913751627e-1;
	rkm->Bhat[7]  = 2.46734519059988698196468570407;
	rkm->Bhat[8]  = -2.54694165184190873912738007542;
	rkm->Bhat[9]  = 1.44354858367677524030187495069;
	rkm->Bhat[10] = 7.94155958811272872713019541622e-2;
	rkm->Bhat[11] = 4.44444444444444444444444444445e-2;
	rkm->Bhat[12] = 0.;
}

/*
	CURTIS810
	
	Runge-Kutta Curtis 8-10
	Source: "High-order explicit Runge–Kutta formulae, their uses, and limitations"
	A. R. Curtis, J. Inst. Math. Appl. 16 (1975), no. 1, 35–55. MR
*/
void Curtis810(RKMethod *rkm) {

	rkm->nstep = 21;
	rkm->alpha = 8.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0]  = 0.;
	rkm->C[1]  = .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401;
	rkm->C[2]  = .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401;
	rkm->C[3]  = .2178778440474225776426322793155049980471766876994920591257703112975964683410769510101;
	rkm->C[4]  = .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252;
	rkm->C[5]  = .6536335321422677329278968379465149941415300630984761773773109338927894050232308530303;
	rkm->C[6]  = .2746594919905254008808021630247618520892150865127407293922085868737635475402543533498;
	rkm->C[7]  = .7735775201106609448405825008093973718589542913426807556412662673054607938029043386501;
	rkm->C[8]  = .5801831400829957086304368756070480288942157185070105667309497004790955953521782539876;
	rkm->C[9]  = .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383;
	rkm->C[10] = .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092;
	rkm->C[11] = .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908;
	rkm->C[12] = .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383;
	rkm->C[13] = .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617;
	rkm->C[14] = .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092;
	rkm->C[15] = .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908;
	rkm->C[16] = .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617;
	rkm->C[17] = 1.;
	rkm->C[18] = .3510848126232741617357001972386587771203155818540433925049309664694280078895463510848;
	rkm->C[19] = .6157407407407407407407407407407407407407407407407407407407407407407407407407407407407;
	rkm->C[20] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0)   = .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401;
	AINDP(2,0)   = .7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1;
	AINDP(2,1)   = .7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1;
	AINDP(3,0)   = .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252e-1;
	AINDP(3,1)   = 0.;
	AINDP(3,2)   = .1634083830355669332319742094866287485353825157746190443443277334731973512558077132576;
	AINDP(4,0)   = .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252;
	AINDP(4,1)   = 0.;
	AINDP(4,2)   = -2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720;
	AINDP(4,3)   = 2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720;
	AINDP(5,0)   = .6536335321422677329278968379465149941415300630984761773773109338927894050232308530303e-1;
	AINDP(5,1)   = 0.;
	AINDP(5,2)   = 0.;
	AINDP(5,3)   = .3268167660711338664639484189732574970707650315492380886886554669463947025116154265151;
	AINDP(5,4)   = .2614534128569070931711587351786059976566120252393904709509243735571157620092923412121;
	AINDP(6,0)   = .8233707757482716585173454344310125296066814318521742241762319051772963627695955263034e-1;
	AINDP(6,1)   = 0.;
	AINDP(6,2)   = 0.;
	AINDP(6,3)   = .2119171963202803561687843468555305553175658807629274312902985594840086570224567152664;
	AINDP(6,4)   = -.3997343508054218311577932550061320162379840049816347807630118786107674477850206579628e-1;
	AINDP(6,5)   = .2037865317596006197606259822674324543477946306275935376058802473310199901934015124941e-1;
	AINDP(7,0)   = .8595305779007343831562027786771081909543936570474230618236291858949564375587825985001e-1;
	AINDP(7,1)   = 0.;
	AINDP(7,2)   = 0.;
	AINDP(7,3)   = 0.;
	AINDP(7,4)   = 0.;
	AINDP(7,5)   = .2911769478058850960337179621761553399856026049598393013981874594942289837064329700000;
	AINDP(7,6)   = .3964475145147024104912442607655312127779123206780991480607158892217361663405931088001;
	AINDP(8,0)   = .8612093485606967549983047372292119178898514571588438099912534616486575243508895957628e-1;
	AINDP(8,1)   = 0.;
	AINDP(8,2)   = 0.;
	AINDP(8,3)   = 0.;
	AINDP(8,4)   = 0.;
	AINDP(8,5)   = .1397464826824442089036313891001189801074425314582326737716288563521183595455090268480;
	AINDP(8,6)   = .3951098495815674599900526056001284215294125840404176924334653987770478924197803010468;
	AINDP(8,7)   = -.4079412703708563576307759281612056453162454270752418047326990081493640904820003348350e-1;
	AINDP(9,0)   = .7233144422337948077616348229119326315582930871089020733092900891206129381937795204778e-1;
	AINDP(9,1)   = 0.;
	AINDP(9,2)   = 0.;
	AINDP(9,3)   = 0.;
	AINDP(9,4)   = 0.;
	AINDP(9,5)   = .2200276284689998102140972735735070061373242800181187459951219347361114857342828430157;
	AINDP(9,6)   = .8789533425436734013369780264792573637952226487753296416823846876217040795688489371334e-1;
	AINDP(9,7)   = -.4445383996260350863990674880611108986832860648196030000580004690002268108984238641730e-1;
	AINDP(9,8)   = -.2183282289488754689095532966861839909872150913926337371522805434288481649401165594213;
	AINDP(10,0)  = .8947100936731114228785441966773836169071038390882857211057269158522704971585365845223e-1;
	AINDP(10,1)  = 0.;
	AINDP(10,2)  = 0.;
	AINDP(10,3)  = 0.;
	AINDP(10,4)  = 0.;
	AINDP(10,5)  = .3946008170285561860741397654755022300929434262701385530048127140223687993778661654316;
	AINDP(10,6)  = .3443011367963333487713764986067104675654371857504670290688086760696354596195596354011;
	AINDP(10,7)  = -.7946682664292661290694938113119430997053815140863772328764150866582492425892231395780e-1;
	AINDP(10,8)  = -.3915218947895966123834967996391962853380545808840091268064277812752553499114569444180;
	AINDP(10,9)  = 0.;
	AINDP(11,0)  = .3210006877963209212945282736072241886741425314298532400216927262619488479186214523312e-1;
	AINDP(11,1)  = 0.;
	AINDP(11,2)  = 0.;
	AINDP(11,3)  = 0.;
	AINDP(11,4)  = 0.;
	AINDP(11,5)  = 0.;
	AINDP(11,6)  = 0.;
	AINDP(11,7)  = -.1846375997512050141835163881753227910996323204749769226655464078048769505209525299752e-3;
	AINDP(11,8)  = .1560894025313219860759149162557283383430181475726228517203663063649626288079337909898;
	AINDP(11,9)  = .1934496857654560252749984220385188727138526287670744309970093278715606577140084022992;
	AINDP(11,10) = .2611612387636636496908928477536452288263163392010050661129958478089356710938164130987;
	AINDP(12,0)  = .4423749328524996327035388417792688154433173133294892285295756457561276315648477233732e-1;
	AINDP(12,1)  = 0.;
	AINDP(12,2)  = 0.;
	AINDP(12,3)  = 0.;
	AINDP(12,4)  = 0.;
	AINDP(12,5)  = 0.;
	AINDP(12,6)  = 0.;
	AINDP(12,7)  = .4640774434539039636406222168781981616534115643208114455689698789119941732444857047798e-2;
	AINDP(12,8)  = .4704660282615136532130927218172390570903230981414159347904277946537920001824903276586e-1;
	AINDP(12,9)  = .8620749948011488160369445167416002799205317397013619044391270706339561700281526529703e-1;
	AINDP(12,10) = -.2607983024682138093233254079066687623148682426317395111719299641390118652802949600035e-1;
	AINDP(12,11) = -.3858020174396621532493277639159499581333235076531298977820093139813399390137768850940e-1;
	AINDP(13,0)  = .2318046717429411567006043539613275607940758021709332569729352990777336390158311630529e-1;
	AINDP(13,1)  = 0.;
	AINDP(13,2)  = 0.;
	AINDP(13,3)  = 0.;
	AINDP(13,4)  = 0.;
	AINDP(13,5)  = 0.;
	AINDP(13,6)  = 0.;
	AINDP(13,7)  = .3197856784116367067302124322582100058864027838197120089129330601737324659881765852593;
	AINDP(13,8)  = .5933233331841898686063939886797828376866051205773280426848164018120869674204443797948;
	AINDP(13,9)  = -1.937519548878479314706815782408229952008442222624773168771865465659822020582450444783;
	AINDP(13,10) = .1803950557030502357344063195737827904476240180662764468232042537858892203518134072359;
	AINDP(13,11) = -.4554014298857220726863505256926549022316460712353658688873150702827663762861750674926;
	AINDP(13,12) = 2.158764106255762807077594619172645539322916635447781333204724468181634037726021280742;
	AINDP(14,0)  = .2624364325798105891527733985858552391723553030719144065844544880498188553839263944447e-1;
	AINDP(14,1)  = 0.;
	AINDP(14,2)  = 0.;
	AINDP(14,3)  = 0.;
	AINDP(14,4)  = 0.;
	AINDP(14,5)  = 0.;
	AINDP(14,6)  = 0.;
	AINDP(14,7)  = .4863139423867266106526843913609225996253073727381961544415263239431571586043622332760e-1;
	AINDP(14,8)  = .4274382538346478867636942429421724367591866585774144180215122660980822123988151132213e-1;
	AINDP(14,9)  = -.4862259869465547771298976981868643277396586803130813159599600102115609499827986711663;
	AINDP(14,10) = .1326047194917652331781527125743684254490968718259563958293167893998110899691451568372;
	AINDP(14,11) = -.9402962152946515651634831658142934852383791641671387741034606371378082209616938685225e-1;
	AINDP(14,12) = .6993864679941022534190304512277131176659196396138275832136258135631963192299339871223;
	AINDP(14,13) = -.1197020013028860976492784934312243036670658451195397948726104511062042521592125912599e-1;
	AINDP(15,0)  = .5568066641536216461090823068917803436066365804361903532125349474551476120813558125830e-1;
	AINDP(15,1)  = 0.;
	AINDP(15,2)  = 0.;
	AINDP(15,3)  = 0.;
	AINDP(15,4)  = 0.;
	AINDP(15,5)  = 0.;
	AINDP(15,6)  = 0.;
	AINDP(15,7)  = -.4324853319508358432896036654421685136736530810118924113940744870078036705505610668088;
	AINDP(15,8)  = -.9979726994172038714656907882931844552238093285811791155499130927685987422432191170216;
	AINDP(15,9)  = 2.707893755718926115778725270396739994070337972517006747100005607751792006959604868323;
	AINDP(15,10) = -1.024823023512132929313567156576969954855232272749038347671818195935585095295127839150;
	AINDP(15,11) = 1.334565206642246959252239602313589265188981560552694580059808406200559397799055652161;
	AINDP(15,12) = -2.587748998830690939658228913150922979184368065866213469477796089200252812362701917187;
	AINDP(15,13) = .8992773696348355846430438306111181223414632598285854300924423251352733205187087732678e-1;
	AINDP(15,14) = 1.497578446211167333777988534023066333042434967475357134513165331964695787890042760189;
	AINDP(16,0)  = -.8434891199686377639125188391985671318383858641413517143104162188088468627447515172982e-3;
	AINDP(16,1)  = 0.;
	AINDP(16,2)  = 0.;
	AINDP(16,3)  = 0.;
	AINDP(16,4)  = 0.;
	AINDP(16,5)  = 0.;
	AINDP(16,6)  = 0.;
	AINDP(16,7)  = .7602144218856081893754106886111596435015500427480120290148318740899211421773423234728;
	AINDP(16,8)  = 1.769083927820959377467464871522349066447068428702073590698445112684989184432409492025;
	AINDP(16,9)  = -4.499239797622297101452915424261016593995695456495268863455643396071539024609271033574;
	AINDP(16,10) = 1.490558190212043468817221563278239942209691100326719140478588601720867838040211450448;
	AINDP(16,11) = -2.552203480132132516997563217309689292804518121743365818482497611667126218719069737195;
	AINDP(16,12) = 4.795167551528575994217413424533259845001657006088189480440731104737960266616292993321;
	AINDP(16,13) = -.9161854401769482236671414092387917470686251714192236693920061138984202381209109248553e-1;
	AINDP(16,14) = -1.525735678746850818217653470352135651821164556169070505816135230784807058389577753184;
	AINDP(16,15) = .7371445601564892133467497107205798584829803038168267854389817508169123996459113657504;
	AINDP(17,0)  = .1017366974111576638766809656369828971944080018220332809259398740674738807023371082700;
	AINDP(17,1)  = 0.;
	AINDP(17,2)  = 0.;
	AINDP(17,3)  = 0.;
	AINDP(17,4)  = 0.;
	AINDP(17,5)  = 0.;
	AINDP(17,6)  = 0.;
	AINDP(17,7)  = -1.696217553209432810711666838709742166182992092906177246174096517233561845662947862824;
	AINDP(17,8)  = -3.825235846211624254528740857512255693551264719132875740261231165548583482101116676418;
	AINDP(17,9)  = 9.754768979885866648856431516333641627109105703674164986615824197909762854575668793816;
	AINDP(17,10) = -2.520767789227152291196336314591227486393143379933686189126240710041836742414125694941;
	AINDP(17,11) = 5.472417145227780046950992000565734793413395536531652419585004300790370984185945495978;
	AINDP(17,12) = -9.781098113458736121002383874108051372067873053264954833376114258940736444388841687929;
	AINDP(17,13) = .3189152692455334369024560213486753019540464785641163242047782111839399471147176681561;
	AINDP(17,14) = 3.447227036527756718156475010324322155277035924051392880570525223655410460762027138915;
	AINDP(17,15) = -.6051983612219277832241707671295607127814820499715293613761402732652780120810041653591;
	AINDP(17,16) = .3334525350307787459202631378414806560287636505658634784117511174230383993073398823363;
	AINDP(18,0)  = -.1012987737478284424676828232882617689682012456457322189102956361570156443805900941944;
	AINDP(18,1)  = 0.;
	AINDP(18,2)  = 0.;
	AINDP(18,3)  = 0.;
	AINDP(18,4)  = 0.;
	AINDP(18,5)  = -.2409389328948775401304659380663043147167897928467308244359962659633933617326533285822e-1;
	AINDP(18,6)  = -.6679880790275182076676283582867036095782150170801495251932447614617249253864579543857;
	AINDP(18,7)  = 1.600262798493100648047998296908183265688507618079976446601985464092263571149154964705;
	AINDP(18,8)  = 3.706958893826695766827011000213884379914407774639901049574259778345288538246990591819;
	AINDP(18,9)  = -8.581755560147929325446798534254342948628755672447282004336563881429983605741487870996;
	AINDP(18,10) = .5607314974300953986559644699099897253584501767603091982484141468619493221310582281877e-1;
	AINDP(18,11) = -4.547761497422899514520768375507009011918601407646237921467449197008085790456674001879;
	AINDP(18,12) = 9.255775439941294621826928846245618922061242300726600002589630404152665447900428712156;
	AINDP(18,13) = -.3450876657451631707159097079770789925142348071643902737346329921538351794816584861003;
	AINDP(18,14) = 0.;
	AINDP(18,15) = 0.;
	AINDP(18,16) = 0.;
	AINDP(18,17) = 0.;
	AINDP(19,0)  = .3826909723812638609001259641818040193828105314579492422836388985468479567237561247336e-1;
	AINDP(19,1)  = 0.;
	AINDP(19,2)  = 0.;
	AINDP(19,3)  = 0.;
	AINDP(19,4)  = 0.;
	AINDP(19,5)  = .7786978965202527814624406274393101840018332461648638653990700950184871893714491273096;
	AINDP(19,6)  = .4859454140913448249612202172501868752761599132465501266008866131088163955018926230543;
	AINDP(19,7)  = 1.814925350154666364151014269029611427420766367555858499108920245656959783343309816408;
	AINDP(19,8)  = 4.551165245704657956889158854062833952834232753889932986749613143631480116805870313264;
	AINDP(19,9)  = -7.173770670344544101351160462586215092596352548535380880420409450623251883641801862305;
	AINDP(19,10) = -.3943009017000923237232456850787591816773705728833192412204243696911216045268772747196;
	AINDP(19,11) = -6.036544185898100312430357626685382432626027303329497026597513524312479466987506315664;
	AINDP(19,12) = 7.338904299721887701527380004651998686389416058019429466200740313593568240326087171554;
	AINDP(19,13) = -.4143158595971836110248598960027762194900538872022960061452263646470675916118824501965;
	AINDP(19,14) = 0.;
	AINDP(19,15) = 0.;
	AINDP(19,16) = 0.;
	AINDP(19,17) = 0.;
	AINDP(19,18) = -.3732349451502749258108621577582478607301443393311959731632798508493352335121760204375;
	AINDP(20,0)  = .2162339046022045866878628785550588026780578552494608097931198882276791962912244674840e-1;
	AINDP(20,1)  = 0.;
	AINDP(20,2)  = 0.;
	AINDP(20,3)  = 0.;
	AINDP(20,4)  = 0.;
	AINDP(20,5)  = .4611834700744369218866370212060318930941322187829670117414118166503940620998117275429;
	AINDP(20,6)  = .1940797759547798743610542713744618433967649025379792966207862125676964319674160574624;
	AINDP(20,7)  = .7041001229739959807963554405302474570280838416767002383409508232534658577705201658489;
	AINDP(20,8)  = 2.877431096792763528910415905652149398266490601780194388811216042455337979365709745445;
	AINDP(20,9)  = 0.;
	AINDP(20,10) = -.4332742088749107411735902392606181444105337491234912425673655059805456011404518143074;
	AINDP(20,11) = -2.234178753588834452567105459024473991729105867012210449973082203886376638514123583334;
	AINDP(20,12) = .2235678086885984010238782832657956960650576194069632574873732156942360146780276407657;
	AINDP(20,13) = .1293532338308457711442786069651741293532338308457711442786069651741293532338308457711;
	AINDP(20,14) = 0.;
	AINDP(20,15) = 0.;
	AINDP(20,16) = 0.;
	AINDP(20,17) = 0.;
	AINDP(20,18) = .1418136968194278394808045812385429206355105705182818920178205766092934777719870449624;
	AINDP(20,19) = -1.085699633131323582531514699802817081967439754938101617737029931360398856861850276906;

	// 10th order solution
	rkm->B[0]  = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1;
	rkm->B[1]  = 0.;
	rkm->B[2]  = 0.;
	rkm->B[3]  = 0.;
	rkm->B[4]  = 0.;
	rkm->B[5]  = 0.;
	rkm->B[6]  = 0.;
	rkm->B[7]  = 0.;
	rkm->B[8]  = 0.;
	rkm->B[9]  = 0.;
	rkm->B[10] = 0.;
	rkm->B[11] = .1387145942588715882541801312803271702142521598590204181697361204933422401935856968980;
	rkm->B[12] = .1892374781489234901583064041060123262381623469486258303271944256799821862794952728707;
	rkm->B[13] = .9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1;
	rkm->B[14] = .2774291885177431765083602625606543404285043197180408363394722409866844803871713937960;
	rkm->B[15] = .1387145942588715882541801312803271702142521598590204181697361204933422401935856968980;
	rkm->B[16] = .9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1;
	rkm->B[17] = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1;
	rkm->B[18] = 0.;
	rkm->B[19] = 0.;
	rkm->B[20] = 0.;

	// 8th order solution
	rkm->Bhat[0]  = .3339829895931337572271945815422988633728883413227543303554098429202731077409488318421e-1;
	rkm->Bhat[1]  = 0.;
	rkm->Bhat[2]  = 0.;
	rkm->Bhat[3]  = 0.;
	rkm->Bhat[4]  = 0.;
	rkm->Bhat[5]  = 0.;
	rkm->Bhat[6]  = 0.;
	rkm->Bhat[7]  = 0.;
	rkm->Bhat[8]  = .5024509803921568627450980392156862745098039215686274509803921568627450980392156862745e-1;
	rkm->Bhat[9]  = -.1423859191318858946753152353981644782061337055184060977838998119673893661279423564924;
	rkm->Bhat[10] = .2126013199429258434998789109063801828540550730541648287733608913970804891935883227446;
	rkm->Bhat[11] = .3254854965632843133622967470840062095221514741629108993207015882688341071771214986692;
	rkm->Bhat[12] = .3312629399585921325051759834368530020703933747412008281573498964803312629399585921325;
	rkm->Bhat[13] = .1887845809230650005639203350759631573314744764356665687950807917985316096487827551997;
	rkm->Bhat[14] = 0.;
	rkm->Bhat[15] = 0.;
	rkm->Bhat[16] = 0.;
	rkm->Bhat[17] = 0.;
	rkm->Bhat[18] = .6159811094287144604404508847679200761569154839698701533406779753675999506388462962070e-1;
	rkm->Bhat[19] = -.9440109660594088037957791636147830082275023098021053120999763935859315478744850552959e-1;
	rkm->Bhat[20] = .3341117040855897708234682470384970584684876341854831047975628586614323631403861184369e-1;
}

/*
	HIROSHI912
	Runge-Kutta Hiroshi 9-12.
	Source: "On the 25 stage 12th order explicit Runge–Kutta method"
	H. Ono, Trans. Japan Soc. Ind. Appl. Math. 16 (2006), no. 3, 177, In Japanese.
*/
void Hiroshi912(RKMethod *rkm)
{
	rkm->nstep = 29;
	rkm->alpha = 9.;

	// Allocate Butcher's tableau
	rkm->C     = (double *)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double *)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double *)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double *)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = NULL;
	rkm->Bphat = NULL;

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0]  = 0.;
	rkm->C[1]  = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852;
	rkm->C[2]  = .4429824561403508771929824561403508771929824561403508771929824561403508771929824561404;
	rkm->C[3]  = .6644736842105263157894736842105263157894736842105263157894736842105263157894736842105;
	rkm->C[4]  = .1069403994175161223216143124609943831911795298522987310172664863740378614520490950697;
	rkm->C[5]  = .1644736842105263157894736842105263157894736842105263157894736842105263157894736842105;
	rkm->C[6]  = .5843251088534107402031930333817126269956458635703918722786647314949201741654571843251;
	rkm->C[7]  = .6382358235823582358235823582358235823582358235823582358235823582358235823582358235824e-1;
	rkm->C[8]  = .2;
	rkm->C[9]  = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
	rkm->C[10] = .9446116054065563496847375720237340591364440847545279548637407620203329280139199114088;
	rkm->C[11] = .5179584680428461035835615460185207449756996828425551866437683588600807071775316523299e-1;
	rkm->C[12] = .8488805186071653506398389301626743020641481756400195420459339398355773991365476236893e-1;
	rkm->C[13] = .2655756032646428930981140590456168352972012641640776214486652703185222349414361456016;
	rkm->C[14] = .5;
	rkm->C[15] = .7344243967353571069018859409543831647027987358359223785513347296814777650585638543984;
	rkm->C[16] = .9151119481392834649360161069837325697935851824359980457954066060164422600863452376311;
	rkm->C[17] = .9446116054065563496847375720237340591364440847545279548637407620203329280139199114088;
	rkm->C[18] = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
	rkm->C[19] = .2;
	rkm->C[20] = .5843251088534107402031930333817126269956458635703918722786647314949201741654571843251;
	rkm->C[21] = .1644736842105263157894736842105263157894736842105263157894736842105263157894736842105;
	rkm->C[22] = .4429824561403508771929824561403508771929824561403508771929824561403508771929824561404;
	rkm->C[23] = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852;
	rkm->C[24] = 1.;
	rkm->C[25] = .4970267001007476028032930885363848318550815967070143979104342007017543082576391973337;
	rkm->C[26] = .8043478260869565217391304347826086956521739130434782608695652173913043478260869565217;
	rkm->C[27] = .8717948717948717948717948717948717948717948717948717948717948717948717948717948717949;
	rkm->C[28] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0)   = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852;
	AINDP(2,0)   = .2175227402212137286104398734798923400326123258875071216675507357419304139407870179368;
	AINDP(2,1)   = .2254597159191371485825425826604585371603701302528437555254317203984204632521954382036;
	AINDP(3,0)   = .1661184210526315789473684210526315789473684210526315789473684210526315789473684210526;
	AINDP(3,1)   = 0.;
	AINDP(3,2)   = .4983552631578947368421052631578947368421052631578947368421052631578947368421052631579;
	AINDP(4,0)   = .8681163193918508865080629713483674348617570159685880673878476405659784672427434203840e-1;
	AINDP(4,1)   = 0.;
	AINDP(4,2)   = .3456981948164897000739691452644149026866903085576649936607708344667115529230058847943e-1;
	AINDP(4,3)   = -.1444105200331793633658889920028385056366520260032657508759536112923114056452583544813e-1;
	AINDP(5,0)   = .3850951504524952575244726520324404860898823844950484490138343541601000605579674908780e-1;
	AINDP(5,1)   = 0.;
	AINDP(5,2)   = 0.;
	AINDP(5,3)   = .9889604363651382462812798900870548282288802143553903392828423536490053424426345577709e-4;
	AINDP(5,4)   = .1258652731216402762123982910182735616976625577395859318541619645591514091994326716669;
	AINDP(6,0)   = .5247404461891304721708365639844356354386624470737997647397329116437616959059854701778;
	AINDP(6,1)   = 0.;
	AINDP(6,2)   = 0.;
	AINDP(6,3)   = .7610651429965941990931946190296373660533922130208900643092560661920652146268815424885e-1;
	AINDP(6,4)   = -2.135538596825204678401302463372339475781988628449683891486081018353295131152366999963;
	AINDP(6,5)   = 2.119016745189825526524339470866652730733632823644186992594087231585247087949150559862;
	AINDP(7,0)   = .3572122856624484555412656871410309454657254402954350529356994466234212462699695879099e-1;
	AINDP(7,1)   = 0.;
	AINDP(7,2)   = 0.;
	AINDP(7,3)   = 0.;
	AINDP(7,4)   = .4596205641509305161863483311675360886293026632692961086374314175410528945129931429031e-1;
	AINDP(7,5)   = -.1800016957713219812376781513215261066354472582005125759344515744644655947135707070493e-1;
	AINDP(7,6)   = .1404669540301245333646491248782654898654978218139650184903068535815036288843799818693e-3;
	AINDP(8,0)   = .1888176809184108413680419706237505704281864015973418413153299364663524364341955734913e-1;
	AINDP(8,1)   = 0.;
	AINDP(8,2)   = 0.;
	AINDP(8,3)   = 0.;
	AINDP(8,4)   = 0.;
	AINDP(8,5)   = .8381198329740939383758467268684128106638290469173528975872944409244863365387981657651e-1;
	AINDP(8,6)   = .9031585241436450659339626877786556727134244569192456020771509670971435444573825121328e-5;
	AINDP(8,7)   = .9729721702550808557495179062390587533407132090396133365371679075124515126725605224924e-1;
	AINDP(9,0)   = -.4080067703469846387001599007585094416645285179656696882006645696781274941696191338570e-1;
	AINDP(9,1)   = 0.;
	AINDP(9,2)   = 0.;
	AINDP(9,3)   = 0.;
	AINDP(9,4)   = 0.;
	AINDP(9,5)   = -.6005539308646711394173617417282635784603029781694346062080009050806164919590318269156;
	AINDP(9,6)   = .1585222658367901823092810289995660365558277114283363946284113044047141269675809741481e-2;
	AINDP(9,7)   = .3026658851734428583799288358366834830670493200050903757742452562250879940584583879948;
	AINDP(9,8)   = .6704368334008921764176894190107687125274815661799611686408713261126274393811928758984;
	AINDP(10,0)  = 6.344326927733666873696263058980778771328558370914152084380122688263116113189631874466;
	AINDP(10,1)  = 0.;
	AINDP(10,2)  = 0.;
	AINDP(10,3)  = 0.;
	AINDP(10,4)  = 0.;
	AINDP(10,5)  = 0.;
	AINDP(10,6)  = 1.975263319684766813850955500728188106557880392559127720135666889111189136253813975729;
	AINDP(10,7)  = -13.82822337504897849061527131653164721900473695395858210726164571154383357718661916508;
	AINDP(10,8)  = 14.82423926991066347292103944504014293055447712090877566105421769783705551827529470930;
	AINDP(10,9)  = -8.370994536873562320168249116193728530299734845668945403444620801647194262518201483006;
	AINDP(11,0)  = -.9910783781470375735018424706990342484478259963427121190616343432259493286903501070750e-1;
	AINDP(11,1)  = 0.;
	AINDP(11,2)  = 0.;
	AINDP(11,3)  = 0.;
	AINDP(11,4)  = 0.;
	AINDP(11,5)  = -1.046319958132641377892211373672215159710875883372049858922517534871447317958500644882;
	AINDP(11,6)  = -.4662578801256029967826301701912208605727776822301916928051722653634090636623902615221e-3;
	AINDP(11,7)  = .4243518408197806737520879584209906999975132950474577957967837636354467847566731274758;
	AINDP(11,8)  = .8350767448298191909277817680995907515365477233544915021061336714983604939879298121048;
	AINDP(11,9)  = -.6204379562043795620437956204379562043795620437956204379562043795620437956204379562044e-1;
	AINDP(11,10) = .3051106025934401220442410373760488176964149504195270785659801678108314263920671243326e-3;
	AINDP(12,0)  = .1731635805893703684448829089734942070445562382116745420114094326722705203174775140332e-1;
	AINDP(12,1)  = 0.;
	AINDP(12,2)  = 0.;
	AINDP(12,3)  = 0.;
	AINDP(12,4)  = 0.;
	AINDP(12,5)  = 0.;
	AINDP(12,6)  = 0.;
	AINDP(12,7)  = 0.;
	AINDP(12,8)  = .8690027159880925811795093256330430220184719404902985010762237927304116502630885425916e-3;
	AINDP(12,9)  = -.9198044746158460357290146202828770811642343120551843217420771155968903404182588708338e-4;
	AINDP(12,10) = .1833594777349928706130754092439075395240645268074916590065144720493485758190652981058e-6;
	AINDP(12,11) = .6679448817377525524901838117990401028051762116902291244289142812068791591710992924480e-1;
	AINDP(13,0)  = .1497702285787817250603134142375551527881374815295532946701659971640303704469333256722e-1;
	AINDP(13,1)  = 0.;
	AINDP(13,2)  = 0.;
	AINDP(13,3)  = 0.;
	AINDP(13,4)  = 0.;
	AINDP(13,5)  = 0.;
	AINDP(13,6)  = 0.;
	AINDP(13,7)  = 0.;
	AINDP(13,8)  = .1299053198125687130806056561360945267563487713700848629960379786258472215118648828622;
	AINDP(13,9)  = .4444252090193421415880477613428889941728587016143899156135432416265141679149665438354e-2;
	AINDP(13,10) = -.2185380330434085597487818850115824768660155325792369439283645364044817815712239904765e-5;
	AINDP(13,11) = .6235770138025566019351747531432668514647398310437277172705875380249580436468314885283e-1;
	AINDP(13,12) = .5389349250407735998767659637686133399860483467584655047185578940287507515886082812093e-1;
	AINDP(14,0)  = .1956886038861791180930007029621719896471785620210633245797118818081004119873004962774;
	AINDP(14,1)  = 0.;
	AINDP(14,2)  = 0.;
	AINDP(14,3)  = 0.;
	AINDP(14,4)  = 0.;
	AINDP(14,5)  = 0.;
	AINDP(14,6)  = 0.;
	AINDP(14,7)  = 0.;
	AINDP(14,8)  = 1.132878041352190960246053936692850204846641933665191609652843227265803479941175715645;
	AINDP(14,9)  = .9754686357770750130309903534180324613492470735016258589778436907188769670738143543989;
	AINDP(14,10) = .4607633555393391875625052489818876130584589587566695260680043303699764560744838206512e-3;
	AINDP(14,11) = -.4833419929869440350491944569648386990066349467457051943375195864364996227787596872297;
	AINDP(14,12) = .2837255261146423361717978265203569959020643630988946562818139270625464117810614651744;
	AINDP(14,13) = -1.604879577498682731680210867877554840351555444499826924680761144749197624460666828087;
	AINDP(15,0)  = -.6628581660952109109950216641584960475848453023977418371554567406302394052049642993363;
	AINDP(15,1)  = 0.;
	AINDP(15,2)  = 0.;
	AINDP(15,3)  = 0.;
	AINDP(15,4)  = 0.;
	AINDP(15,5)  = 0.;
	AINDP(15,6)  = 0.;
	AINDP(15,7)  = 0.;
	AINDP(15,8)  = -5.301219753823165631138506283750981613687493380976143446419698251739201243786272433269;
	AINDP(15,9)  = -5.493744530005151771950832780438489668352801158207831962667055860430416893664761589330;
	AINDP(15,10) = .6448107716343851659685115053410322161204197332247306929859781244843876115320017653542e-2;
	AINDP(15,11) = 2.226911096857986220827171118161838780361655511952909244621862410327793805385281866476;
	AINDP(15,12) = -.8260094546883369994121948528168152647220748853112540011094140493170779782006783468229;
	AINDP(15,13) = 9.736973734199539440165209948098457305642565726332220268712929325968631361580234056276;
	AINDP(15,14) = 1.047923362573352907746375340805459350884588027111516805638308114257144242834404582751;
	AINDP(16,0)  = 9.451896878619703179615895147719645120861585484250220485534769741441232402797149004402;
	AINDP(16,1)  = 0.;
	AINDP(16,2)  = 0.;
	AINDP(16,3)  = 0.;
	AINDP(16,4)  = 0.;
	AINDP(16,5)  = 0.;
	AINDP(16,6)  = 0.;
	AINDP(16,7)  = 0.;
	AINDP(16,8)  = 74.07837676951819392708315469048815926453257791083353986411194358002277837032529121349;
	AINDP(16,9)  = 80.08971633421252928663933464653279322909492698803384552922422585968931770071180510457;
	AINDP(16,10) = -.1241702484260160323471778040612055855620194337087768444994146945534716737633295008233;
	AINDP(16,11) = -32.04108125365225923900659816438329678215226279318337501510105142207810131441066638158;
	AINDP(16,12) = 15.51919421000708125838831543558938623465692911684926717234376165413330242668122537342;
	AINDP(16,13) = -136.4444237346563024309541869705527862321862115141467999073204884062451779712319267938;
	AINDP(16,14) = -11.36109896858298168444633695622136642878717316645309025294277498926997344146696251530;
	AINDP(16,15) = 1.746701961099335199963616081872403749335232589961167014444435282876535760443759733215;
	AINDP(17,0)  = 1.059086740089530572084686602455620391446019980075022020682291957643879144548290466243;
	AINDP(17,1)  = 0.;
	AINDP(17,2)  = 0.;
	AINDP(17,3)  = 0.;
	AINDP(17,4)  = 0.;
	AINDP(17,5)  = 0.;
	AINDP(17,6)  = 1.975263319684766813850955500728188106557880392559127720135666889111189136253813975729;
	AINDP(17,7)  = -13.82822337504897849061527131653164721900473695395858210726164571154383357718661916508;
	AINDP(17,8)  = -26.72676722061492102175215840759265382960786498230207636418025793583833107088763635551;
	AINDP(17,9)  = -53.49798986553787152824955918158085193982743760760492627590775850823540448031315562431;
	AINDP(17,10) = .7179812487812968675250181644918743516316099586081748670146732996043783143974128607587e-1;
	AINDP(17,11) = 18.05559723664965139501142471462605880255358278669785143935434349294122071206578011157;
	AINDP(17,12) = -8.765163819232982113792324392383185125858432915348234208512916467037776994687248464210;
	AINDP(17,13) = 76.87522358555575651988028075659565616302805855797883225625477685501643635288672453023;
	AINDP(17,14) = 6.541007260506904941709851688460080123863085658071943624909334991972860537373864701434;
	AINDP(17,15) = -.8461837296739539683436767973109554382501241857205608248520628152174372445350689400749;
	AINDP(17,16) = .3096334815052354314802658810823658907325235844531318754050068324709258105543338930304e-1;
	AINDP(18,0)  = -.9776673260040047976710026469624134014469714141641816253185153958371390334954366206102e-1;
	AINDP(18,1)  = 0.;
	AINDP(18,2)  = 0.;
	AINDP(18,3)  = 0.;
	AINDP(18,4)  = 0.;
	AINDP(18,5)  = -.6005539308646711394173617417282635784603029781694346062080009050806164919590318269156;
	AINDP(18,6)  = .1585222658367901823092810289995660365558277114283363946284113044047141269675809741481e-2;
	AINDP(18,7)  = .3026658851734428583799288358366834830670493200050903757742452562250879940584583879948;
	AINDP(18,8)  = .4285248952296136656501504776128834912424841890282487604085652387338305228046948965999;
	AINDP(18,9)  = -.2494936049956423479331713062821820760626388763661853215936426940474748269791918051828e-1;
	AINDP(18,10) = .2011342149874705950638155217740310555133571768677965217932407885434178872356469888792e-1;
	AINDP(18,11) = .1265809640150575117740069284741749647793461966567090378701137188708013785394651655531;
	AINDP(18,12) = -.7625505511021937468044360160133448561674800889340550320120549846220783051800508863794e-2;
	AINDP(18,13) = .2257683215130023640661938534278959810874704491725768321513002364066193853427895981087;
	AINDP(18,14) = -.2745492489167708483879226183906024865124767290578864277433486220883760763766017543321e-1;
	AINDP(18,15) = .7783761260627937338705544952370637303437265442198231246348540202794422441986844495010e-2;
	AINDP(18,16) = .1680830476266175444138000910160343423421707378563384687234091879582812751321589095322e-4;
	AINDP(18,17) = -.2135549195295375067495229039525877007339581782873806009604806379884885927885993014667e-1;
	AINDP(19,0)  = .4918722777774214635454694806095471334564471116346374258920200467475606642975663633169e-1;
	AINDP(19,1)  = 0.;
	AINDP(19,2)  = 0.;
	AINDP(19,3)  = 0.;
	AINDP(19,4)  = 0.;
	AINDP(19,5)  = .8381198329740939383758467268684128106638290469173528975872944409244863365387981657651e-1;
	AINDP(19,6)  = .9031585241436450659339626877786556727134244569192456020771509670971435444573825121328e-5;
	AINDP(19,7)  = .9729721702550808557495179062390587533407132090396133365371679075124515126725605224924e-1;
	AINDP(19,8)  = -.7780700010332922911083090560520956762627768671838912755501179322945892692289272010345e-1;
	AINDP(19,9)  = .2322816011687352870985186646759446719321113455516507860765563208248738749278233611786;
	AINDP(19,10) = -.1069713522971685698540328329852038445947544770807964964519088069861375575524890552230e-1;
	AINDP(19,11) = -.1273357855008010164642613326804825893357151796749893640146110897003842493364412495521;
	AINDP(19,12) = .1195596306045770335798408097988199315656849089154404333183111529360783290389646594900;
	AINDP(19,13) = .1268882175226586102719033232628398791540785498489425981873111782477341389728096676737;
	AINDP(19,14) = .2044989775051124744376278118609406952965235173824130879345603271983640081799591002045e-1;
	AINDP(19,15) = -.3909856506984242164691157929655038995858933888334921275260056794557284704385445582377e-2;
	AINDP(19,16) = -.9117410119332842109125953032188840693633119317454401715027062207348542639139803659651e-5;
	AINDP(19,17) = .1126093151077456154585623418197514536662349793737574398626301058923382946852045056475e-1;
	AINDP(19,18) = -.3209868434922071245903287586373535845929558438862699119277778588606558307508436673461;
	AINDP(20,0)  = .5247404461891304721708365639844356354386624470737997647397329116437616959059854701778;
	AINDP(20,1)  = 0.;
	AINDP(20,2)  = 0.;
	AINDP(20,3)  = .7610651429965941990931946190296373660533922130208900643092560661920652146268815424885e-1;
	AINDP(20,4)  = -2.135538596825204678401302463372339475781988628449683891486081018353295131152366999963;
	AINDP(20,5)  = 2.119016745189825526524339470866652730733632823644186992594087231585247087949150559862;
	AINDP(20,6)  = 0.;
	AINDP(20,7)  = 0.;
	AINDP(20,8)  = -.6068669292751351984511141135476525247489731841233590933089036540462466843434512014723;
	AINDP(20,9)  = -.6975023816048118989159659127206041759560640169551283418432821453966450135301718130674;
	AINDP(20,10) = -.2539552135383387231958801508125376392283214913791972115937570999107334296436460657368e-1;
	AINDP(20,11) = 0.;
	AINDP(20,12) = 0.;
	AINDP(20,13) = 0.;
	AINDP(20,14) = 0.;
	AINDP(20,15) = 0.;
	AINDP(20,16) = 0.;
	AINDP(20,17) = .2539552135383387231958801508125376392283214913791972115937570999107334296436460657368e-1;
	AINDP(20,18) = .6975023816048118989159659127206041759560640169551283418432821453966450135301718130674;
	AINDP(20,19) = .6068669292751351984511141135476525247489731841233590933089036540462466843434512014723;
	AINDP(21,0)  = .3850951504524952575244726520324404860898823844950484490138343541601000605579674908780e-1;
	AINDP(21,1)  = 0.;
	AINDP(21,2)  = 0.;
	AINDP(21,3)  = .9889604363651382462812798900870548282288802143553903392828423536490053424426345577709e-4;
	AINDP(21,4)  = .1258652731216402762123982910182735616976625577395859318541619645591514091994326716669;
	AINDP(21,5)  = 0.;
	AINDP(21,6)  = .1148546657824708136802241649776446744025406100582655444423341878411637297467351936831;
	AINDP(21,7)  = 0.;
	AINDP(21,8)  = -.1299246462925958527077352114079744607750386324042347411360003531859939834460301561195;
	AINDP(21,9)  = -.3664591598580916638621456622089859363144583004356560608229075984337239720284817960193;
	AINDP(21,10) = 0.;
	AINDP(21,11) = 0.;
	AINDP(21,12) = 0.;
	AINDP(21,13) = 0.;
	AINDP(21,14) = 0.;
	AINDP(21,15) = 0.;
	AINDP(21,16) = 0.;
	AINDP(21,17) = 0.;
	AINDP(21,18) = .3664591598580916638621456622089859363144583004356560608229075984337239720284817960193;
	AINDP(21,19) = .1299246462925958527077352114079744607750386324042347411360003531859939834460301561195;
	AINDP(21,20) = -.1148546657824708136802241649776446744025406100582655444423341878411637297467351936831;
	AINDP(22,0)  = .2175227402212137286104398734798923400326123258875071216675507357419304139407870179368;
	AINDP(22,1)  = .2254597159191371485825425826604585371603701302528437555254317203984204632521954382036;
	AINDP(22,2)  = 0.;
	AINDP(22,3)  = 0.;
	AINDP(22,4)  = 0.;
	AINDP(22,5)  = -.7003676470588235294117647058823529411764705882352941176470588235294117647058823529412;
	AINDP(22,6)  = -.3841432262252079340469853296205751191707035527677198880717241749265840367609120672261;
	AINDP(22,7)  = 0.;
	AINDP(22,8)  = 0.;
	AINDP(22,9)  = 0.;
	AINDP(22,10) = 0.;
	AINDP(22,11) = 0.;
	AINDP(22,12) = 0.;
	AINDP(22,13) = 0.;
	AINDP(22,14) = 0.;
	AINDP(22,15) = 0.;
	AINDP(22,16) = 0.;
	AINDP(22,17) = 0.;
	AINDP(22,18) = 0.;
	AINDP(22,19) = 0.;
	AINDP(22,20) = .3841432262252079340469853296205751191707035527677198880717241749265840367609120672261;
	AINDP(22,21) = .7003676470588235294117647058823529411764705882352941176470588235294117647058823529412;
	AINDP(23,0)  = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852;
	AINDP(23,1)  = 0.;
	AINDP(23,2)  = -.4244806610219170956648818689597945762515945523075081224120777274585003570971510665439;
	AINDP(23,3)  = 0.;
	AINDP(23,4)  = 0.;
	AINDP(23,5)  = 0.;
	AINDP(23,6)  = 0.;
	AINDP(23,7)  = 0.;
	AINDP(23,8)  = 0.;
	AINDP(23,9)  = 0.;
	AINDP(23,10) = 0.;
	AINDP(23,11) = 0.;
	AINDP(23,12) = 0.;
	AINDP(23,13) = 0.;
	AINDP(23,14) = 0.;
	AINDP(23,15) = 0.;
	AINDP(23,16) = 0.;
	AINDP(23,17) = 0.;
	AINDP(23,18) = 0.;
	AINDP(23,19) = 0.;
	AINDP(23,20) = 0.;
	AINDP(23,21) = 0.;
	AINDP(23,22) = .4244806610219170956648818689597945762515945523075081224120777274585003570971510665439;
	AINDP(24,0)  = 14.54990971513478580914495526598245283940024810265087163382498896386800769812364397026;
	AINDP(24,1)  = -2.609444444444444444444444444444444444444444444444444444444444444444444444444444444444;
	AINDP(24,2)  = -2.016004609236637754870351028563643794559738431497207211298306162299623087053267335725;
	AINDP(24,3)  = 0.;
	AINDP(24,4)  = 0.;
	AINDP(24,5)  = -1.666875;
	AINDP(24,6)  = -1.840010137609049715480551028603992780751854184812582730319893278211480604766702319850;
	AINDP(24,7)  = 0.;
	AINDP(24,8)  = 112.8850026879393650927243496060343039004107661597559344830054269802155904451112935844;
	AINDP(24,9)  = 123.3942086822776167534198661922997568027612066045886506351026838994566451157760161079;
	AINDP(24,10) = -.7912126656078716671206667350992937018785263673658058027818457299775466324638082312137;
	AINDP(24,11) = -50.05149873558555185701853277950418028083920384013973228797374663569746882026502756132;
	AINDP(24,12) = 24.88778291494286023489265916235844707524170235927746532159915897878512972222657954044;
	AINDP(24,13) = -212.1164197057320847511194729602871445996396409264694299135694995667783941778631290781;
	AINDP(24,14) = -17.89082255740024165397920079197116480379826062999613139163493416384514698556677051724;
	AINDP(24,15) = 2.509716434086569985363077956674238286724267256774032757480361264001771307725305753519;
	AINDP(24,16) = .1162475886937088360535666476065014131326443901611275096836079807311057386735832665899;
	AINDP(24,17) = .5840328281597421751021178306186155103094570763962374837158281176580253234011237044692;
	AINDP(24,18) = 1.584417806712708397483834073813399023446091584920600008797782958694409008929749703075;
	AINDP(24,19) = 1.338635006378392645053446531474068534729248229446179562750186952887872256191439757181;
	AINDP(24,20) = 1.840010137609049715480551028603992780751854184812582730319893278211480604766702319850;
	AINDP(24,21) = 1.666875;
	AINDP(24,22) = 2.016004609236637754870351028563643794559738431497207211298306162299623087053267335725;
	AINDP(24,23) = 2.609444444444444444444444444444444444444444444444444444444444444444444444444444444444;
	AINDP(25,0) = .4213659219087082450668941175898637044576077571412870755532166081975824639155816025423;
	AINDP(25,1) = 0.;
	AINDP(25,2)  = 0.;
	AINDP(25,3)  = 0.;
	AINDP(25,4)  = 0.;
	AINDP(25,5)  = 2.360375290413766425107807321597798068298832670876354152068674897991709951684675057597;
	AINDP(25,6)  = .7887926811836902144270824477231365437442962856383491851349713361086542561144935964760e-1;
	AINDP(25,7)  = -1.881850641776530466652474803895333308262409829541779704822719782535259533895036456964;
	AINDP(25,8)  = -1.304700734906095391371228323883517431348033016096248869612640797091741378769699510853;
	AINDP(25,9)  = .1146971532060496506611311517441641422873299688900135018999537656144573134542692746612;
	AINDP(25,10) = -.5223613182942077907170609676910338480915071906000186351215047373740452623507436746018e-2;
	AINDP(25,11) = .7134840563194221964556259902880063405282394887795535106616674222878805188799073074497;
	AINDP(25,12) = 0.;
	AINDP(25,13) = 0.;
	AINDP(25,14) = 0.;
	AINDP(25,15) = 0.;
	AINDP(25,16) = 0.;
	AINDP(25,17) = 0.;
	AINDP(25,18) = 0.;
	AINDP(25,19) = 0.;
	AINDP(25,20) = 0.;
	AINDP(25,21) = 0.;
	AINDP(25,22) = 0.;
	AINDP(25,23) = 0.;
	AINDP(25,24) = 0.;
	AINDP(26,0)  = -1.016867684065179179311540011641152067739527559831057381007428484343486238202737995538;
	AINDP(26,1)  = 0.;
	AINDP(26,2)  = 0.;
	AINDP(26,3)  = 0.;
	AINDP(26,4)  = 0.;
	AINDP(26,5)  = -7.712044352285817603610736737545203003182107799250377304475316992902646516184627861475;
	AINDP(26,6)  = -.4034008409374858753410643280039779311650023076266296210030337872976327092081717454860;
	AINDP(26,7)  = 6.739165476490825275476530741799137001781411805688541615342712456896969009075069743117;
	AINDP(26,8)  = 6.014994643407224294180918860565568523540180411761624603152605084567461778121724411633;
	AINDP(26,9)  = -1.138427387973993086846707441740657423236451007119331997503516596250482644801941463866;
	AINDP(26,10) = .5009271973181599563449431362188397685579770251743883251930235916805543539186803311245e-1;
	AINDP(26,11) = -3.113250932564715585587456369457050044282973515833171733791753110672694235302567710668;
	AINDP(26,12) = 0.;
	AINDP(26,13) = 0.;
	AINDP(26,14) = 0.;
	AINDP(26,15) = 0.;
	AINDP(26,16) = 0.;
	AINDP(26,17) = 0.;
	AINDP(26,18) = 0.;
	AINDP(26,19) = 0.;
	AINDP(26,20) = 0.;
	AINDP(26,21) = 0.;
	AINDP(26,22) = 0.;
	AINDP(26,23) = 0.;
	AINDP(26,24) = 0.;
	AINDP(26,25) = 1.384086184284282287144691407184059663080846182736441247635994288225760468937471545692;
	AINDP(27,0)  = 1.131093475949031458408970675798323789793651098141584053672723684523337489654810813567;
	AINDP(27,1)  = 0.;
	AINDP(27,2)  = 0.;
	AINDP(27,3)  = 0.;
	AINDP(27,4)  = 0.;
	AINDP(27,5)  = -11.30475611955440577592346561419842170756276153921836803995513790246783417109406727002;
	AINDP(27,6)  = .8673508908529372037894544277195364499375277491231106383198729858086287127092489733738e-1;
	AINDP(27,7)  = 4.971317844154333915807514558966554931901059536397301577948574837963960265600948244094;
	AINDP(27,8)  = 14.86493772010299652718002500847699984298963300991290479767277341466986550270499895786;
	AINDP(27,9)  = -5.526130551905351405702373768620234518347747212552226033118137855980983127699224086148;
	AINDP(27,10) = .1017790491986200061558195486579246543940163889995857946318495764475820398967567810759;
	AINDP(27,11) = -5.412708567655345677389304794550103449135846140313886894023804835964789906470796078495;
	AINDP(27,12) = 0.;
	AINDP(27,13) = 0.;
	AINDP(27,14) = 0.;
	AINDP(27,15) = 0.;
	AINDP(27,16) = 0.;
	AINDP(27,17) = 0.;
	AINDP(27,18) = 0.;
	AINDP(27,19) = 0.;
	AINDP(27,20) = 0.;
	AINDP(27,21) = 0.;
	AINDP(27,22) = 0.;
	AINDP(27,23) = 0.;
	AINDP(27,24) = 0.;
	AINDP(27,25) = 2.119905903216124397337756706998226742489167533804374626293399155463227287948398852330;
	AINDP(27,26) = -.1603789707964253713820928925063521366431305782887091520824325014403564569409562398070;
	AINDP(28,0)  = 46.12864603958015905056850990838704062569763496465412763519385202873225169326091596022;
	AINDP(28,1)  = 0.;
	AINDP(28,2)  = 0.;
	AINDP(28,3)  = 0.;
	AINDP(28,4)  = 0.;
	AINDP(28,5)  = 27.91300163119399908845158457840426795358131126287389180101909720743727165524743997096;
	AINDP(28,6)  = 16.11362689862451240990975288339484000039234961533373225370151713993919933068137837655;
	AINDP(28,7)  = -125.4696763444318726329250646477825268481685879278990572088898486327607587168793188547;
	AINDP(28,8)  = 76.57182020120529497684089567511627659347626021577320126137924234653281598557298427556;
	AINDP(28,9)  = -48.97805558723490361747755876313897556229903489002597907918723603724236906416246304792;
	AINDP(28,10) = -1.242830487244052672528847080627776989497284066852925470069357906835696725179561610962;
	AINDP(28,11) = 18.85807213383620068645464308722546866214730809725606975205204214283143331463602851729;
	AINDP(28,12) = 0.;
	AINDP(28,13) = 0.;
	AINDP(28,14) = 0.;
	AINDP(28,15) = 0.;
	AINDP(28,16) = 0.;
	AINDP(28,17) = 0.;
	AINDP(28,18) = 0.;
	AINDP(28,19) = 0.;
	AINDP(28,20) = 0.;
	AINDP(28,21) = 0.;
	AINDP(28,22) = 0.;
	AINDP(28,23) = 0.;
	AINDP(28,24) = 0.;
	AINDP(28,25) = -8.871982194511738170929283011936752038683182063824370821380217475185779921363890310462;
	AINDP(28,26) = -2.069534982695615656321598541301254038890059866116256078922522295738656144906973495868;
	AINDP(28,27) = 2.046912691678016537956965912259391642243284658827565955103431482290288593093460219351;

	// 12th order solution
	rkm->B[0]  = .2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1;
	rkm->B[1]  = -.11;
	rkm->B[2]  = -.17;
	rkm->B[3]  = 0.;
	rkm->B[4]  = 0.;
	rkm->B[5]  = -.19;
	rkm->B[6]  = -.21;
	rkm->B[7]  = 0.;
	rkm->B[8]  = -.23;
	rkm->B[9]  = -.27;
	rkm->B[10] = -.29;
	rkm->B[11] = 0.;
	rkm->B[12] = .1384130236807829740053502031450331467488136400899412345912671194817223119377730668077;
	rkm->B[13] = .2158726906049313117089355111406811389654720741957730511230185948039919737765126474781;
	rkm->B[14] = .2438095238095238095238095238095238095238095238095238095238095238095238095238095238095;
	rkm->B[15] = .2158726906049313117089355111406811389654720741957730511230185948039919737765126474781;
	rkm->B[16] = .1384130236807829740053502031450331467488136400899412345912671194817223119377730668077;
	rkm->B[17] = .29;
	rkm->B[18] = .27;
	rkm->B[19] = .23;
	rkm->B[20] = .21;
	rkm->B[21] = .19;
	rkm->B[22] = .17;
	rkm->B[23] = .11;
	rkm->B[24] = .2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1;
	rkm->B[25] = 0.;
	rkm->B[26] = 0.;
	rkm->B[27] = 0.;
	rkm->B[28] = 0.;

	// 9th order solution
	rkm->Bhat[0]  = .1357267366422036691624508570375039921213961405197383540990729271585608575744367221422e-1;
	rkm->Bhat[1]  = 0.;
	rkm->Bhat[2]  = 0.;
	rkm->Bhat[3]  = 0.;
	rkm->Bhat[4]  = 0.;
	rkm->Bhat[5]  = 0.;
	rkm->Bhat[6]  = 0.;
	rkm->Bhat[7]  = 0.;
	rkm->Bhat[8]  = .1957242608025905233613155193355413222756712243739956539857744660017249157522413539338;
	rkm->Bhat[9]  = .6188866347435608661616060238380345321740579185698262925072441030810461130118224626650e-1;
	rkm->Bhat[10] = .2356461254963383884566009189531765712282459575681722654363888116993755654359289982043;
	rkm->Bhat[11] = .9356981277656948171659125355126822887878370575133833270918675400840988929847712528189e-1;
	rkm->Bhat[12] = 0.;
	rkm->Bhat[13] = 0.;
	rkm->Bhat[14] = 0.;
	rkm->Bhat[15] = 0.;
	rkm->Bhat[16] = 0.;
	rkm->Bhat[17] = 0.;
	rkm->Bhat[18] = 0.;
	rkm->Bhat[19] = 0.;
	rkm->Bhat[20] = 0.;
	rkm->Bhat[21] = 0.;
	rkm->Bhat[22] = 0.;
	rkm->Bhat[23] = 0.;
	rkm->Bhat[24] = 0.;
	rkm->Bhat[25] = .2788382624223597882496809755901865993371492355322230983832374605891024297949723722505;
	rkm->Bhat[26] = .4265887719284871852244002213010235951116369449606049493109912099005098392090245967637;
	rkm->Bhat[27] = -.2878025166474501962999477241233861487151201542084751225511851437798538088192303308430;
	rkm->Bhat[28] = -.1802605391747162424104685269536402054591231988681564193502526144322952773004003407203e-1;
}

/*
	RUNGEKUTTANYSTROM34

	Runge-Kutta-Nystrom 3-4.
	Source: “High Phase-Lag-Order Runge--Kutta and Nyström Pairs,”
			S. N. Papakostas and C. Tsitouras,  
			SIAM J. Sci. Comput., vol. 21, no. 2, pp. 747–763, Jan. 1999.
*/
void RungeKuttaNystrom34(RKMethod *rkm) {

	rkm->nstep = 4;
	rkm->alpha = 3.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bphat = (double*)malloc(rkm->nstep*sizeof(double));

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.0e0;
	rkm->C[1] = 5280320246./8739982487.;
	rkm->C[2] = 2488284716./4549257065.;
	rkm->C[3] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = 379540831./2079644486.;
	AINDP(2,0) = 4090178078./30784316359.;
	AINDP(2,1) = 724764885./43347833381.;
	AINDP(3,0) = 1382866212./8058611935.;
	AINDP(3,1) = -958895179./4232828403.;
	AINDP(3,2) = 3246882200./5850906049.;

	// 4th order solution
	rkm->B[0] = 1382866212./8058611935.;
	rkm->B[1] = -958895179./4232828403.;
	rkm->B[2] = 3246882200./5850906049.;
	rkm->B[3] = 0.0e0;

	// 3rd order solution
	rkm->Bhat[0] = 2788682442./15739802773.;
	rkm->Bhat[1] = 1966645306./8825920383.;
	rkm->Bhat[2] = 3./20.;
	rkm->Bhat[3] = -1./20.;

	// 4th order solution (derivative)
	rkm->Bp[0]  = 1382866212./8058611935.;
	rkm->Bp[1]  = -5103277582./8917268805.;
	rkm->Bp[2]  = 12283001905./10027502947.;
	rkm->Bp[3]  = 2640160123./15021458620.;

	// 3rd order solution (derivative)
	rkm->Bphat[0]  = 1153046258./2574473725.;
	rkm->Bphat[1]  = 27711564209./4540932828.;
	rkm->Bphat[2]  = -31509757256./6039640179.;
	rkm->Bphat[3]  = -1./3.;
}

/*
	RUNGEKUTTANYSTROM46

	Runge-Kutta-Nystrom 4-6.
	Source: “High Phase-Lag-Order Runge--Kutta and Nyström Pairs,”
			S. N. Papakostas and C. Tsitouras,  
			SIAM J. Sci. Comput., vol. 21, no. 2, pp. 747–763, Jan. 1999.
*/
void RungeKuttaNystrom46(RKMethod *rkm) {

	rkm->nstep = 6;
	rkm->alpha = 4.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bphat = (double*)malloc(rkm->nstep*sizeof(double));

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.;
	rkm->C[1] = 76064096./555208869.;
	rkm->C[2] = 61651457./172436989.;
	rkm->C[3] = 473./677.;
	rkm->C[4] = 1521284172./2494038851.;
	rkm->C[5] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = 148104835./15781657211.;
	AINDP(2,0) = 83570507./15621004272.;
	AINDP(2,1) = 1008730685./17224387836.;
	AINDP(3,0) = 313507335./5002407628.;
	AINDP(3,1) = 232561219./6632504445.;
	AINDP(3,2) = 2487592367./16999280447.;
	AINDP(4,0) = 497059253./7416116119.;
	AINDP(4,1) = -31814195./4301239521.;
	AINDP(4,2) = 666859859./4681708182.;
	AINDP(4,3) = -164695106./10269973361.;
	AINDP(5,0) = 1104491309./17344385380.;
	AINDP(5,1) = 2297852298./21296988487.;
	AINDP(5,2) = 1270882233./4862169760.;
	AINDP(5,3) = 1745513301./8827769149.;
	AINDP(5,4) = -854905921./6541617807.;

	// 6th order solution
	rkm->B[0] = 1104491309./17344385380.;
	rkm->B[1] = 2297852298./21296988487.;
	rkm->B[2] = 1270882233./4862169760.;
	rkm->B[3] = 1745513301./8827769149.;
	rkm->B[4] = -854905921./6541617807.;
	rkm->B[5] = 0.;

	// 4th order solution
	rkm->Bhat[0] = 390850314./4665518297.;
	rkm->Bhat[1] = 879866760./14012015573.;
	rkm->Bhat[2] = 1237986347./4111942715.;
	rkm->Bhat[3] = 1838896521./8824986790.;
	rkm->Bhat[4] = -1945509358./12470194255.;
	rkm->Bhat[5] = 0.;

	// 6th order solution (derivative)
	rkm->Bp[0] = 1104491309./17344385380.;
	rkm->Bp[1] = 928753894./7428602053.;
	rkm->Bp[2] = 1088487657./2675475233.;
	rkm->Bp[3] = 2281030107./3476164510.;
	rkm->Bp[4] = -5717085047./17062458528.;
	rkm->Bp[5] = 1./12.;
	
	// 4th order solution (derivative)
	rkm->Bphat[0]  = 390850314./4665518297.;
	rkm->Bphat[1]  = 831255784./11424277409.;
	rkm->Bphat[2]  = 5386054494./11493559817.;
	rkm->Bphat[3]  = 5380034471./7780066871.;
	rkm->Bphat[4]  = -2./5.;
	rkm->Bphat[5]  = 1./12.;
}

/*
	RUNGEKUTTANYSTROM68

	Runge-Kutta-Nystrom 6-8.
	Source: “High-order zero-dissipative Runge-Kutta-Nystrom methods”
			C. Tsitouras 
			J. Comput. Appl. Math., vol. 95, no. 1–2, pp. 157–161, 1998.
*/
void RungeKuttaNystrom68(RKMethod *rkm) {

	rkm->nstep = 9;
	rkm->alpha = 6.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bphat = (double*)malloc(rkm->nstep*sizeof(double));

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0] = 0.0e0;
	rkm->C[1] = 53510637./572121662.;
	rkm->C[2] = 357248860./1909448607.;
	rkm->C[3] = 59418235./171833923.;
	rkm->C[4] = 504173341./1087994762.;
	rkm->C[5] = 387843305./481838798.;
	rkm->C[6] = 862622967./973409174.;
	rkm->C[7] = 1.;
	rkm->C[8] = 1.;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0) = 2232527./510224080.;
	AINDP(2,0) = 2232527./382668060.;
	AINDP(2,1) = 37081334./3177977723.;
	AINDP(3,0) = 30566341./1760496148.;
	AINDP(3,1) = 8218085./734916954.;
	AINDP(3,2) = 20946667./670502569.;
	AINDP(4,0) = 8338478./576439607.;
	AINDP(4,1) = 261273394./6072168265.;
	AINDP(4,2) = 29462594./997634791.;
	AINDP(4,3) = 77204245./3795243917.;
	AINDP(5,0) = 1100589730./2256756531.;
	AINDP(5,1) = -2504085737./1611855801.;
	AINDP(5,2) = 7903856905./4158249776.;
	AINDP(5,3) = -2746138419./2848012111.;
	AINDP(5,4) = 131209591./289472778.;
	AINDP(6,0) = -5574225260./1895534583.;
	AINDP(6,1) = 18077935161./1758392437.;
	AINDP(6,2) = -9350911807./861383162;
	AINDP(6,3) = 1948485627./355652546.;
	AINDP(6,4) = -736004998./452939571;
	AINDP(6,5) = 30627356./562404525.;
	AINDP(7,0) = 1325770706./1152764413.;
	AINDP(7,1) = -4839505982./1296410665.;
	AINDP(7,2) = 3190857006./755292433.;
	AINDP(7,3) = -4573493832./2249363015.;
	AINDP(7,4) = 1400717651./1651772814.;
	AINDP(7,5) = 12114323./331609812.;
	AINDP(7,6) = 3175573./456323904.;
	AINDP(8,0) = 13093931./251590199.;
	AINDP(8,1) = 0.;
	AINDP(8,2) = 453730323./1723073015.;
	AINDP(8,3) = -174729768./1541576453.;
	AINDP(8,4) = 154302939./629373676.;
	AINDP(8,5) = 89112923./1899373863.;
	AINDP(8,6) = 3795253./644583773.;
	AINDP(8,7) = 0.;

	// 8th order solution
	rkm->B[0] = 13093931./251590199.;
	rkm->B[1] = 0.;
	rkm->B[2] = 453730323./1723073015.;
	rkm->B[3] = -174729768./1541576453.;
	rkm->B[4] = 154302939./629373676.;
	rkm->B[5] = 89112923./1899373863.;
	rkm->B[6] = 3795253./644583773.;
	rkm->B[7] = 0.;
	rkm->B[8] = 0.;

	// 6th order solution
	rkm->Bhat[0] = 106709197./1949525980.;
	rkm->Bhat[1] = 0.;
	rkm->Bhat[2] = 189148257./765523804.;
	rkm->Bhat[3] = -85259905./1101513319.;
	rkm->Bhat[4] = 281276824./1286133115.;
	rkm->Bhat[5] = 25824341./461156990.;
	rkm->Bhat[6] = 968741./1095289242.;
	rkm->Bhat[7] = 0.;
	rkm->Bhat[8] = 0.;

	// 8th order solution (derivative)
	rkm->Bp[0] = 13093931./251590199.;
	rkm->Bp[1] = 0.;
	rkm->Bp[2] = 236119875./728916961.;
	rkm->Bp[3] = -145423380./839364157.;
	rkm->Bp[4] = 56276012./123171693.;
	rkm->Bp[5] = 168458860./700436693.;
	rkm->Bp[6] = 36197530./699693487.;
	rkm->Bp[7] = 49704941./1032349655.;
	rkm->Bp[8] = 0.;

	// 6th order solution (derivative)
	rkm->Bphat[0] = 106709197./1949525980.;
	rkm->Bphat[1] = 0.;
	rkm->Bphat[2] = 463818847./1525964332.;
	rkm->Bphat[3] = -379038881./3203661980.;
	rkm->Bphat[4] = 461048333./1131231912.;
	rkm->Bphat[5] = 203285488./708159667.;
	rkm->Bphat[6] = 10565889./1359619685.;
	rkm->Bphat[7] = -146565653./1579900056.;
	rkm->Bphat[8] = 3./20.;
}

/*
	RUNGEKUTTANYSTROM1012

	Runge-Kutta-Nystrom 10-12.
	Source: http://www.tampa.phys.ucl.ac.uk/rmat/test/rknint.f
*/
void RungeKuttaNystrom1012(RKMethod *rkm) {

	rkm->nstep = 17;
	rkm->alpha = 10.;

	// Allocate Butcher's tableau
	rkm->C     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->A     = (double*)malloc(rkm->nstep*rkm->nstep*sizeof(double));
	rkm->B     = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bhat  = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bp    = (double*)malloc(rkm->nstep*sizeof(double));
	rkm->Bphat = (double*)malloc(rkm->nstep*sizeof(double));

	rkm->alloc = 1;

	// C coefficient
	rkm->C[0]  = 0.0e0;
	rkm->C[1]  = 2.0e-2;
	rkm->C[2]  = 4.0e-2;
	rkm->C[3]  = 1.0e-1;
	rkm->C[4]  = 1.33333333333333333333333333333e-1;
	rkm->C[5]  = 1.6e-1;
	rkm->C[6]  = 5.0e-2;
	rkm->C[7]  = 2.0e-1;
	rkm->C[8]  = 2.5e-1;
	rkm->C[9]  = 3.33333333333333333333333333333e-1;
	rkm->C[10] = 5.0e-1;
	rkm->C[11] = 5.55555555555555555555555555556e-1;
	rkm->C[12] = 7.5e-1;
	rkm->C[13] = 8.57142857142857142857142857143e-1;
	rkm->C[14] = 9.45216222272014340129957427739e-1;
	rkm->C[15] = 1.0e0;
	rkm->C[16] = 1.0e0;

	// A matrix
	memset(rkm->A,0,rkm->nstep*rkm->nstep*sizeof(double));
	AINDP(1,0)   = 2.0e-4;
	AINDP(2,0)   = 2.66666666666666666666666666667e-4;
	AINDP(2,1)   = 5.33333333333333333333333333333e-4;
	AINDP(3,0)   = 2.91666666666666666666666666667e-3;
	AINDP(3,1)   = -4.16666666666666666666666666667e-3;
	AINDP(3,2)   = 6.25e-3;
	AINDP(4,0)   = 1.64609053497942386831275720165e-3;
	AINDP(4,1)   = 0.0e0;
	AINDP(4,2)   = 5.48696844993141289437585733882e-3;
	AINDP(4,3)   = 1.75582990397805212620027434842e-3;
	AINDP(5,0)   = 1.9456e-3;
	AINDP(5,1)   = 0.0e0;
	AINDP(5,2)   = 7.15174603174603174603174603175e-3;
	AINDP(5,3)   = 2.91271111111111111111111111111e-3;
	AINDP(5,4)   = 7.89942857142857142857142857143e-4;
	AINDP(6,0)   = 5.6640625e-4;
	AINDP(6,1)   = 0.0e0;
	AINDP(6,2)   = 8.80973048941798941798941798942e-4;
	AINDP(6,3)   = -4.36921296296296296296296296296e-4;
	AINDP(6,4)   = 3.39006696428571428571428571429e-4;
	AINDP(6,5)   = -9.94646990740740740740740740741e-5;
	AINDP(7,0)   = 3.08333333333333333333333333333e-3;
	AINDP(7,1)   = 0.0e0;
	AINDP(7,2)   = 0.0e0;
	AINDP(7,3)   = 1.77777777777777777777777777778e-3;
	AINDP(7,4)   = 2.7e-3;
	AINDP(7,5)   = 1.57828282828282828282828282828e-3;
	AINDP(7,6)   = 1.08606060606060606060606060606e-2;
	AINDP(8,0)   = 3.65183937480112971375119150338e-3;
	AINDP(8,1)   = 0.0e0;
	AINDP(8,2)   = 3.96517171407234306617557289807e-3;
	AINDP(8,3)   = 3.19725826293062822350093426091e-3;
	AINDP(8,4)   = 8.22146730685543536968701883401e-3;
	AINDP(8,5)   = -1.31309269595723798362013884863e-3;
	AINDP(8,6)   = 9.77158696806486781562609494147e-3;
	AINDP(8,7)   = 3.75576906923283379487932641079e-3;
	AINDP(9,0)   = 3.70724106871850081019565530521e-3;
	AINDP(9,1)   = 0.0e0;
	AINDP(9,2)   = 5.08204585455528598076108163479e-3;
	AINDP(9,3)   = 1.17470800217541204473569104943e-3;
	AINDP(9,4)   = -2.11476299151269914996229766362e-2;
	AINDP(9,5)   = 6.01046369810788081222573525136e-2;
	AINDP(9,6)   = 2.01057347685061881846748708777e-2;
	AINDP(9,7)   = -2.83507501229335808430366774368e-2;
	AINDP(9,8)   = 1.48795689185819327555905582479e-2;
	AINDP(10,0)  = 3.51253765607334415311308293052e-2;
	AINDP(10,1)  = 0.0e0;
	AINDP(10,2)  = -8.61574919513847910340576078545e-3;
	AINDP(10,3)  = -5.79144805100791652167632252471e-3;
	AINDP(10,4)  = 1.94555482378261584239438810411e0;
	AINDP(10,5)  = -3.43512386745651359636787167574e0;
	AINDP(10,6)  = -1.09307011074752217583892572001e-1;
	AINDP(10,7)  = 2.3496383118995166394320161088e0;
	AINDP(10,8)  = -7.56009408687022978027190729778e-1;
	AINDP(10,9)  = 1.09528972221569264246502018618e-1;
	AINDP(11,0)  = 2.05277925374824966509720571672e-2;
	AINDP(11,1)  = 0.0e0;
	AINDP(11,2)  = -7.28644676448017991778247943149e-3;
	AINDP(11,3)  = -2.11535560796184024069259562549e-3;
	AINDP(11,4)  = 9.27580796872352224256768033235e-1;
	AINDP(11,5)  = -1.65228248442573667907302673325e0;
	AINDP(11,6)  = -2.10795630056865698191914366913e-2;
	AINDP(11,7)  = 1.20653643262078715447708832536e0;
	AINDP(11,8)  = -4.13714477001066141324662463645e-1;
	AINDP(11,9)  = 9.07987398280965375956795739516e-2;
	AINDP(11,10) = 5.35555260053398504916870658215e-3;
	AINDP(12,0)  = -1.43240788755455150458921091632e-1;
	AINDP(12,1)  = 0.0e0;
	AINDP(12,2)  = 1.25287037730918172778464480231e-2;
	AINDP(12,3)  = 6.82601916396982712868112411737e-3;
	AINDP(12,4)  = -4.79955539557438726550216254291e0;
	AINDP(12,5)  = 5.69862504395194143379169794156e0;
	AINDP(12,6)  = 7.55343036952364522249444028716e-1;
	AINDP(12,7)  = -1.27554878582810837175400796542e-1;
	AINDP(12,8)  = -1.96059260511173843289133255423e0;
	AINDP(12,9)  = 9.18560905663526240976234285341e-1;
	AINDP(12,10) = -2.38800855052844310534827013402e-1;
	AINDP(12,11) = 1.59110813572342155138740170963e-1;
	AINDP(13,0)  = 8.04501920552048948697230778134e-1;
	AINDP(13,1)  = 0.0e0;
	AINDP(13,2)  = -1.66585270670112451778516268261e-2;
	AINDP(13,3)  = -2.1415834042629734811731437191e-2;
	AINDP(13,4)  = 1.68272359289624658702009353564e1;
	AINDP(13,5)  = -1.11728353571760979267882984241e1;
	AINDP(13,6)  = -3.37715929722632374148856475521e0;
	AINDP(13,7)  = -1.52433266553608456461817682939e1;
	AINDP(13,8)  = 1.71798357382154165620247684026e1;
	AINDP(13,9)  = -5.43771923982399464535413738556e0;
	AINDP(13,10) = 1.38786716183646557551256778839e0;
	AINDP(13,11) = -5.92582773265281165347677029181e-1;
	AINDP(13,12) = 2.96038731712973527961592794552e-2;
	AINDP(14,0)  = -9.13296766697358082096250482648e-1;
	AINDP(14,1)  = 0.0e0;
	AINDP(14,2)  = 2.41127257578051783924489946102e-3;
	AINDP(14,3)  = 1.76581226938617419820698839226e-2;
	AINDP(14,4)  = -1.48516497797203838246128557088e1;
	AINDP(14,5)  = 2.15897086700457560030782161561e0;
	AINDP(14,6)  = 3.99791558311787990115282754337e0;
	AINDP(14,7)  = 2.84341518002322318984542514988e1;
	AINDP(14,8)  = -2.52593643549415984378843352235e1;
	AINDP(14,9)  = 7.7338785423622373655340014114e0;
	AINDP(14,10) = -1.8913028948478674610382580129e0;
	AINDP(14,11) = 1.00148450702247178036685959248e0;
	AINDP(14,12) = 4.64119959910905190510518247052e-3;
	AINDP(14,13) = 1.12187550221489570339750499063e-2;
	AINDP(15,0)  = -2.75196297205593938206065227039e-1;
	AINDP(15,1)  = 0.0e0;
	AINDP(15,2)  = 3.66118887791549201342293285553e-2;
	AINDP(15,3)  = 9.7895196882315626246509967162e-3;
	AINDP(15,4)  = -1.2293062345886210304214726509e1;
	AINDP(15,5)  = 1.42072264539379026942929665966e1;
	AINDP(15,6)  = 1.58664769067895368322481964272e0;
	AINDP(15,7)  = 2.45777353275959454390324346975e0;
	AINDP(15,8)  = -8.93519369440327190552259086374e0;
	AINDP(15,9)  = 4.37367273161340694839327077512e0;
	AINDP(15,10) = -1.83471817654494916304344410264e0;
	AINDP(15,11) = 1.15920852890614912078083198373e0;
	AINDP(15,12) = -1.72902531653839221518003422953e-2;
	AINDP(15,13) = 1.93259779044607666727649875324e-2;
	AINDP(15,14) = 5.20444293755499311184926401526e-3;
	AINDP(16,0)  = 1.30763918474040575879994562983e0;
	AINDP(16,1)  = 0.0e0;
	AINDP(16,2)  = 1.73641091897458418670879991296e-2;
	AINDP(16,3)  = -1.8544456454265795024362115588e-2;
	AINDP(16,4)  = 1.48115220328677268968478356223e1;
	AINDP(16,5)  = 9.38317630848247090787922177126e0;
	AINDP(16,6)  = -5.2284261999445422541474024553e0;
	AINDP(16,7)  = -4.89512805258476508040093482743e1;
	AINDP(16,8)  = 3.82970960343379225625583875836e1;
	AINDP(16,9)  = -1.05873813369759797091619037505e1;
	AINDP(16,10) = 2.43323043762262763585119618787e0;
	AINDP(16,11) = -1.04534060425754442848652456513e0;
	AINDP(16,12) = 7.17732095086725945198184857508e-2;
	AINDP(16,13) = 2.16221097080827826905505320027e-3;
	AINDP(16,14) = 7.00959575960251423699282781988e-3;
	AINDP(16,15) = 0.0e0;

	// 12th order solution
	rkm->B[0]  = 1.21278685171854149768890395495e-2;
	rkm->B[1]  = 0.0e0;
	rkm->B[2]  = 0.0e0;
	rkm->B[3]  = 0.0e0;
	rkm->B[4]  = 0.0e0;
	rkm->B[5]  = 0.0e0;
	rkm->B[6]  = 8.62974625156887444363792274411e-2;
	rkm->B[7]  = 2.52546958118714719432343449316e-1;
	rkm->B[8]  = -1.97418679932682303358307954886e-1;
	rkm->B[9]  = 2.03186919078972590809261561009e-1;
	rkm->B[10] = -2.07758080777149166121933554691e-2;
	rkm->B[11] = 1.09678048745020136250111237823e-1;
	rkm->B[12] = 3.80651325264665057344878719105e-2;
	rkm->B[13] = 1.16340688043242296440927709215e-2;
	rkm->B[14] = 4.65802970402487868693615238455e-3;
	rkm->B[15] = 0.0e0;
	rkm->B[16] = 0.0e0;

	// 10th order solution
	rkm->Bhat[0]  = 1.70087019070069917527544646189e-2;
	rkm->Bhat[1]  = 0.0e0;
	rkm->Bhat[2]  = 0.0e0;
	rkm->Bhat[3]  = 0.0e0;
	rkm->Bhat[4]  = 0.0e0;
	rkm->Bhat[5]  = 0.0e0;
	rkm->Bhat[6]  = 7.22593359308314069488600038463e-2;
	rkm->Bhat[7]  = 3.72026177326753045388210502067e-1;
	rkm->Bhat[8]  = -4.01821145009303521439340233863e-1;
	rkm->Bhat[9]  = 3.35455068301351666696584034896e-1;
	rkm->Bhat[10] = -1.31306501075331808430281840783e-1;
	rkm->Bhat[11] = 1.89431906616048652722659836455e-1;
	rkm->Bhat[12] = 2.68408020400290479053691655806e-2;
	rkm->Bhat[13] = 1.63056656059179238935180933102e-2;
	rkm->Bhat[14] = 3.79998835669659456166597387323e-3;
	rkm->Bhat[15] = 0.0e0;
	rkm->Bhat[16] = 0.0e0;

	// 12th order solution (derivative)
	rkm->Bp[0]  = 1.21278685171854149768890395495e-2;
	rkm->Bp[1]  = 0.0e0;
	rkm->Bp[2]  = 0.0e0;
	rkm->Bp[3]  = 0.0e0;
	rkm->Bp[4]  = 0.0e0;
	rkm->Bp[5]  = 0.0e0;
	rkm->Bp[6]  = 9.08394342270407836172412920433e-2;
	rkm->Bp[7]  = 3.15683697648393399290429311645e-1;
	rkm->Bp[8]  = -2.63224906576909737811077273181e-1;
	rkm->Bp[9]  = 3.04780378618458886213892341513e-1;
	rkm->Bp[10] = -4.15516161554298332243867109382e-2;
	rkm->Bp[11] = 2.46775609676295306562750285101e-1;
	rkm->Bp[12] = 1.52260530105866022937951487642e-1;
	rkm->Bp[13] = 8.14384816302696075086493964505e-2;
	rkm->Bp[14] = 8.50257119389081128008018326881e-2;
	rkm->Bp[15] = -9.15518963007796287314100251351e-3;
	rkm->Bp[16] = 2.5e-2;

	// 10th order solution (derivative)
	rkm->Bphat[0]  = 1.70087019070069917527544646189e-2;
	rkm->Bphat[1]  = 0.0e0;
	rkm->Bphat[2]  = 0.0e0;
	rkm->Bphat[3]  = 0.0e0;
	rkm->Bphat[4]  = 0.0e0;
	rkm->Bphat[5]  = 0.0e0;
	rkm->Bphat[6]  = 7.60624588745593757356421093119e-2;
	rkm->Bphat[7]  = 4.65032721658441306735263127583e-1;
	rkm->Bphat[8]  = -5.35761526679071361919120311817e-1;
	rkm->Bphat[9]  = 5.03182602452027500044876052344e-1;
	rkm->Bphat[10] = -2.62613002150663616860563681567e-1;
	rkm->Bphat[11] = 4.26221789886109468625984632024e-1;
	rkm->Bphat[12] = 1.07363208160116191621476662322e-1;
	rkm->Bphat[13] = 1.14139659241425467254626653171e-1;
	rkm->Bphat[14] = 6.93633866500486770090602920091e-2;
	rkm->Bphat[15] = 2.0e-2;
	rkm->Bphat[16] = 0.0e0;
}

/*
CHECK_TABLEAU
Checks if the Butcher's tableau is consistent.
*/
int CheckTableau(RKMethod *rkm)
{
	int retval = 1;
	const double checksum = 200. * 1E-6;

	// Check steps
	double sum1 = 0., sum2 = 0.;
	for (int ii = 0; ii < rkm->nstep; ii++)
	{
		double sum = 0.;
		for (int jj = 0; jj < ii; jj++)
		{
			sum += AINDP(ii, jj);
		}
		if (fabs(rkm->C[ii] - sum) > checksum)
		{
			printf("Error in row %d of Butcher tableau!\n", ii);
			retval = 0;
		}
		sum1 += rkm->B[ii];
		sum2 += rkm->Bhat[ii];
	}

	if (fabs(1. - sum1) > checksum)
	{
		printf("Error in solution high of Butcher tableau!\n");
		retval = 0;
	}

	if (fabs(1. - sum2) > checksum)
	{
		printf("Error in solution low of Butcher tableau!\n");
		retval = 0;
	}

	return retval;
}