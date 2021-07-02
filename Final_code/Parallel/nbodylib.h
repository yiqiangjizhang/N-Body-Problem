
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
#define R4    227943824 	     // [km] Mars distance from Sun
#define R5    778340821 	     // [km] Jupiter distance from Sun
#define R6    1426666422 	     // [km] Saturn distance from Sun
#define R7    2870658186 	     // [km] Uranus distance from Sun
#define R8    4498396441 	     // [km] Neptune distance from Sun
#define R9    5906.4e6 	       // [km] Pluto distance from Sun

#define NBODY 3                // Number of bodies
#define N     6*NBODY          // Number of initial conditions

#define Vx(Mx,Rx) sqrt(G * Mx / Rx)
#define Tx(Mx,Rx) 2*PI/sqrt(G * Mx / Rx / Rx / Rx)

// Macros
#define VAR(b,p) var[6*(b) + (p)]                             // b: body, p: parameter
#define VARP(b,p) varp[6*(b) + (p)]                           // b: body, p: parameter
#define Y(b,p) y[3*(b) + (p)]          		                  // b: body, p: parameter
#define DIST(b1,b2,p) (Y(b2,p) - Y(b1,p))                     // b1: body 1, b2: body 2 p: (0,1,2) = (x,y,z)
#define IGLOB(i,start) ((i) + (start))                        // global planet number, i: local planet number, start: starting planet of the processor
#define BODIES(proc) planet2proc[proc]-planet2proc[proc - 1]  // number of bodies in processor (proc)

// Global vector of mass
extern double MASS[NBODY];

// Prototypes
void checkr(int r, char *txt);
int proc();
int nproc();
void worksplit(int *mystart, int *myend, int proc, int nproc, int start, int end);

// Import variable to this file
extern int* planet2proc;

void NBody(double t, double* var, int n, double* varp);
