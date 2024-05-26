/**
# Numerical code
Id 1 is for the Viscoelastic liquid pool, and Id 2 is Newtonian gas.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
#include "two-phase.h"

// #include "log-conform.h" // main version
#define ViscoElastic
#include "log-conform-elastic_v5.h" // elastic version -> With Deborah number -> \infty
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"

#define tsnap (0.01) // 0.001 only for some cases. 
// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-2)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define OmegaErr (1e-3)                             // error tolerances in vorticity inside the liquid

// Numbers!
#define Ldomain 8

// boundary conditions
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

int MAXlevel;
// Oh -> Solvent Ohnesorge number
// Oha -> air Ohnesorge number
// De -> Deborah number
// Ec -> Elasto-capillary number

double Oh, Oha, Ec, Bond, tmax;
char nameOut[80], dumpFile[80];
scalar Gpd[];

int  main(int argc, char const *argv[]) {
  dtmax = 1e-4; //  BEWARE of this for stability issues. 

  L0 = Ldomain;
  origin (-L0/2., 0.);
  init_grid (1 << 6);

  // Values taken from the terminal
  MAXlevel = atoi(argv[1]);
  Ec = atof(argv[2]);
  Oh = atof(argv[3]);
  Oha = atof(argv[4]);
  Bond = atof(argv[5]);
  tmax = atof(argv[6]);

  // Ensure that all the variables were transferred properly from the terminal or job script.
  if (argc < 7){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 7-argc);
    return 1;
  }
  fprintf(ferr, "Level %d, De Infty, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, Ec, Oh, Oha, Bond);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "dump");


  rho1 = 1., rho2 = 1e-3;
  mu1 = Oh, mu2 = Oha;

  f.sigma = 1.0;

  // polymers
  Gp = Gpd;

  run();
}

event properties (i++) {
  foreach () {
    Gpd[] = (f[] > 1-1e-6 ? Ec: 0.); // this is an artificial patch for now. The code has issues with VE terms in the interfacial cells!
  }
}

int refRegion(double x, double y, double z){
  return (y < 1.28 ? MAXlevel: y < 2.56 ? MAXlevel-1 : y < 5.12 ? MAXlevel-2 : MAXlevel-3);
}

event init (t = 0) {
  if (!restore (file = dumpFile)){

    char filename[60];
    sprintf(filename,"Bo%5.4f.dat",Bond);
    FILE * fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    while (adapt_wavelet_limited ((scalar *){f, d}, (double[]){1e-8, 1e-8}, refRegion).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fraction of the domain. */
    fractions (phi, f);
  }
  // return 1;
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  scalar KAPPA[], omega[];
  curvature(f, KAPPA);
  vorticity (u, omega);
  foreach(){
    omega[] *= f[];
  }
  boundary ((scalar *){KAPPA, omega});
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA, omega},
     (double[]){fErr, VelErr, VelErr, KErr, OmegaErr},
     refRegion);
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Ending Simulation
*/
event end (t = end) {
  fprintf(ferr, "Level %d, De Infty, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, Ec, Oh, Oha, Bond);
}

/**
## Log writing
*/
event logWriting (i++) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf(fp, "Level %d, De Infty, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, Ec, Oh, Oha, Bond);
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

  if (ke > 1e2 && i > 1e2){
    fprintf(ferr, "The kinetic energy blew up. Stopping simulation\n");
    return 1;
  }
  if (ke < 1e-6 && i > 1e2){
    fprintf(ferr, "The kinetic energy is too small. Perhaps you reached a stationary state. Check! Stopping simulation\n");
    return 1;
  }
}
