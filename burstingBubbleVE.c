/*
This is OLD code and will need to be updated!!!

# Introduction:
We investigate the classical problem of bubble bursting at a liquid-gas interface, but now in the presence of a viscoelastic liquid medium. Here are the schematics of the problem. This simualtion will start from Figure 1(c). 

# Numerical code
Id 1 is for the Viscoelastic liquid pool, and Id 2 is Newtonian gas.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
/**
To model Viscoelastic liquids, we use a modified version of [two-phase.h](http://basilisk.fr/src/two-phase.h). [two-phaseVP.h](two-phaseVP.h) contains these modifications.
*/
#include "two-phase.h"
#include "log-conform.h"
// /**
//  You can use: conserving.h as well. Even without it, I was still able to conserve the total energy (also momentum?) of the system if I adapt based on curvature and vorticity/deformation tensor norm (see the adapt even). I had to smear the density and viscosity anyhow because of the sharp ratios in liquid (Bingham) and the gas.
// */
// #include "navier-stokes/conserving.h"
#include "tension.h"
// #include "reduced.h"
#include "distance.h"
/**
We use a modified adapt-wavelet algorithm available [(here)](http://basilisk.fr/sandbox/pairetti/bag_mode/adapt_wavelet_limited.h). It is written by *César Pairetti* (Thanks :)).
*/
#include "adapt_wavelet_limited_v3.h"

#define tsnap (0.001)
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
// Ohs -> Solvent Ohnesorge number
// Ohp -> Polymer Ohnesorge number
// Oha -> air Ohnesorge number
// De -> Deborah number

double Ohs, Ohp, Oha, Bond, De, tmax;
char nameOut[80], dumpFile[80];
scalar mupd[], lam[];

int  main(int argc, char const *argv[]) {
  L0 = Ldomain;
  origin (-L0/2., 0.);
  init_grid (1 << 6);
  // Values taken from the terminal
  MAXlevel = atoi(argv[1]);
  De = atof(argv[2]);
  Ohs = atof(argv[3]);
  Ohp = atof(argv[4]);
  Oha = atof(argv[5]);
  Bond = atof(argv[6]);
  tmax = atof(argv[7]);

  // Ensure that all the variables were transferred properly from the terminal or job script.
  if (argc < 8){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",8-argc);
    return 1;
  }
  fprintf(ferr, "Level %d, De %4.3f, Ohs %2.1e, Ohp %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, De, Ohs, Ohp, Oha, Bond);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "dump");


  rho1 = 1., rho2 = 1e-3;
  mu1 = Ohs, mup = mupd, mu2 = Oha;
  lambda = lam;
  f.sigma = 1.0;
  // G.x = -Bond; // Gravity switched off
  run();
}

event properties (i++) {
  foreach () {
    lam[] = De*clamp(f[], 0.,1.);
    mupd[] = Ohp*clamp(f[], 0.,1.);
  }
}

/**
This event is specific to César's adapt_wavelet_limited.
*/
int refRegion(double x, double y, double z){
  return (y < 1.28 ? MAXlevel : y < 2.56 ? MAXlevel-1 : y < 5.12 ? MAXlevel-2 : MAXlevel-3);
}

/**
## Initial Condition

The initial shape of the bubble at the liquid-gas interface can be calculated by solving the Young-Laplace equations ([Lhuissier & Villermaux, 2012](https://doi.org/10.1017/jfm.2011.418​)) and it depends on the $\mathcal{B}o$ number. 
Resources: 

* [Alex Berny's Sandbox](http://www.basilisk.fr/sandbox/aberny/bubble/bubble.c)
* [My Matlab code:](https://github.com/VatsalSy/Bursting-Bubble-In-a-Viscoelastic-Medium/blob/main/InitialCondition.m) Also see the results for different $\mathcal{B}o$ number [here](https://youtu.be/Z_vdsOW5fsg).

<p align="center">
  <img src="VillermauxComparision.png" width="25%">
  <caption><p align="center">Comparision of the initial shape calculated using the Young-Laplace equations. Ofcourse, for this study, we only need: $\mathcal{B}o = 10^{-3}$. In the figure, $a$ is the capillary length, $a = \sqrt{\gamma/(\rho_lg)}$</caption>
</p>

Since we do not focus on the influence of $\mathcal{B}o$ number in the present study, I am not elaborating on it here. For all the simulations, I use the interfacial shape calculated for $\mathcal{B}o = 10^{-3}$. The resultant data file is available [here](https://raw.githubusercontent.com/VatsalSy/Bursting-Bubble-In-a-Viscoelastic-Medium/main/Bo0.0010.dat).

**Note:** The curvature diverges at the cavity-free surface intersection. We fillet this corner to circumvent this singularity, introducing a rim with a finite curvature that connects the bubble to the free surface. We ensured that the curvature of the rim is high enough such that the subsequent dynamics are independent of its finite value.
*/
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
  return 1;
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  /**
  We adapt based on curvature, $\kappa$ and vorticity $\omega$. 
  Adaptation based on $\kappa$ ensures a constant grid resolution across the interface. See [this](http://basilisk.fr/sandbox/Antoonvh/rc.c) for further reading. 

  We also adapt based on vorticity in the liquid domain. I have noticed that this refinement helps resolve the fake-yield surface accurately (see the black regions in the videos below). 
  */
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
  fprintf(ferr, "Level %d, De %4.3f, Ohs %2.1e, Ohp %2.1e, Oha %2.1e, Bo %4.3f\n", MAXlevel, De, Ohs, Ohp, Oha, Bond);
}

/**
## Log writing
*/
event logWriting (i+=100) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
  // Ensure that the cut-off Kinetic energy is smaller than or equal to 1e-6 times the maximum kinetic energy of the system.
  if (ke > 1e3 || ke < 1e-6){
    if (i > 1e2){
      return 1;
    }
  }
}

/**
## Running the code
~~~bash
#!/bin/bash
qcc -fopenmp -Wall -O2 burstingBubble.c -o burstingBubble -lm
export OMP_NUM_THREADS=8
./burstingBubble 10 0.25 1e-3 1e-2 5.0
~~~
*/