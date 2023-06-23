/* 
   Metropolis Monte Carlo simulation of a 2D Ising system

   Cameron F. Abrams

   Written for the course CHE 800-002, Adv. Chem. Engr. Thermodynamics
   Winter 1011

   compile using "gcc -o ising ising.c -lm"

   runs as "./ising -L <sidelength(20)> -T <temperature(1.0)> \
                    -nc <numcycles(1e6)> -fs <samplefreq(100)> \
		    -s <seed(?)>"

   For example, to run a 20x20 system at T = 0.5 for 1e7 cycles
   sampling every 100 cycles, the command looks like
           
          ./ising -L 20 -T 0.5 -nc 1e7 -fs 100

   Append the flag "-novis" if you do NOT want to watch the magnet!
   
   The default values are shown in parentheses above.

   The Hamiltonian is 

   H = -J sum_<ij> s_i * s_j,

   where "sum_<ij>" is the sum over all unique
   nearest neighbor pairs, s_i = +/- 1, and J 
   is an arbitrary "coupling" parameter having
   units of energy.  We work in a unit system
   where energy is measured in units of J and
   temperature is nondimensionalized as (k_B*T)/J.

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2011
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This function computes and returns the change in
   system energy when spin (i,j) is flipped.  The
   modulo arithmetic (the % operator) ensures 
   periodic boundaries. The syntax "i?(i-1):(L-1)"
   performs the following check:  If i is non-zero,
   return i-1, otherwise return L-1.  This also 
   ensures periodic boundaries.  */
int E ( int ** F, int L, int i, int j ) {
  return -2*(F[i][j])*(F[i?(i-1):(L-1)][j]+F[(i+1)%L][j]+
		      F[i][j?(j-1):(L-1)]+F[i][(j+1)%L]);
}

/* Sample the system; compute the average magnetization
   and the average energy per spin */
double samp ( int ** F, int L, double * s, double * e ) {
  int i,j;

  *s=0.0;
  *e=0.0;
  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      *s+=(double)F[i][j];
      *e-=(double)(F[i][j])*(F[i][(j+1)%L]+F[(i+1)%L][j]);
    }
  }
  *s/=(L*L);
  *e/=(L*L);
}

/* Randomly assigns all spins */
void init ( int ** F, int L ) {
  int i,j;

  /* Visit each position (i,j) in the lattice */
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      /* 2*x-1, where x is randomly 0,1.
         lrand48() returns a non-negative integer, and
         %2 takes its modulo-2; if even, this zero, if odd, 1. */
      F[i][j]=2*(int)(lrand48()%2)-1;
    }
  }
}

void draw ( int ** F, int L, int c, double s, double e, char u, char d, int visual ) {
  int i,j;

  fprintf(stdout,"%c[3;0H",27);
  fprintf(stdout,"cycle % 8d  <s> % 8.5lf   <E>/N % 8.5lf\n",
	  c,s,e);
  if (!visual) return;
  for (i=0;i<L;i++) {
    for (j=0;j<L;j++) {
      fprintf(stdout,"% 2c",(F[i][j]==-1)?d:u);
    }
    fprintf(stdout,"\n");
  }
}

int main (int argc, char * argv[]) {

  /* System parameters */
  int ** F;       /* The 2D array of spins; i.e., the "magnet" */
  int L = 20;     /* The sidelength of the magnet */
  int N;          /* The total number of spins = L*L */
  double T = 1.0; /* Dimensionless temperature = (T*k)/J */

  /* Run parameters */
  int nCycles = 1000000; /* number of MC cycles to run; one cycle is N 
			    consecutive attempted spin flips */
  int fSamp = 1000;      /* Frequency with which samples are taken */

  /* Computational variables */
  int nSamp;      /* Number of samples taken */
  int de;         /* energy change due to flipping a spin */
  double b;       /* Boltzman factor */
  double x;       /* random number */
  int i,j,a,c;    /* loop counters */

  /* Observables */
  double s=0.0, ssum=0.0;    /* average magnetization */
  double e=0.0, esum=0.0;    /* average energy per spin */

  unsigned long int Seed = 23410981;

  char upchar = '#';
  char downchar = ' ';

  int visual = 1;

  /* Seed the pseudorandom number generator */
  srand48(Seed);

  /* Parse command-line arguments */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-L")) L=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-nc")) nCycles = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-s")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-u")) upchar = argv[++i][0];
    else if (!strcmp(argv[i],"-d")) downchar = argv[++i][0];
    else if (!strcmp(argv[i],"-novis")) visual = 0;
  }
  
  /* Output some initial information */
  fprintf(stdout,"%c[2J",27);
  fprintf(stdout,"%c[0;0H",27);
  fprintf(stdout,"# Ising simulation, NVT Metropolis Monte Carlo -- cfa 2011\n");
  fprintf(stdout,"# L = %i, T = %.3lf, nCycles %i, fSamp %i, Seed %u\n",
	  L,T,nCycles,fSamp,Seed);

  /* Compute the number of spins */
  N=L*L;

  /* Allocate memory for the magnet */
  F=(int**)malloc(L*sizeof(int*));
  for (i=0;i<L;i++) F[i]=(int*)malloc(L*sizeof(int));

  /* Generate an initial state */
  init(F,L);

  /* For computational efficiency, convert T to reciprocal T */
  T=1.0/T;

  s = 0.0;
  e = 0.0;
  nSamp = 0;
  for (c=0;c<nCycles;c++) {
    /* Make N flip attempts */
    for (a=0;a<N;a++) {
      /* randomly select a spin */
      i=(int)(drand48()*L);
      j=(int)(drand48()*L);
      /* get the "new" energy as the incremental change due
         to flipping spin (i,j) */
      de = E(F,L,i,j);
      /* compute the Boltzmann factor; recall T is now
         reciprocal temperature */
      b = exp(de*T);
      /* pick a random number between 0 and 1 */
      x = drand48();
      /* accept or reject this flip */
      if (x<b) { /* accept */
	/* flip it */
	F[i][j]*=-1;
      }
    }
    /* Sample and accumulate averages */
    if (!(c%fSamp)) {
      samp(F,L,&s,&e);
      ssum+=s;
      esum+=e;
      nSamp++;
      draw(F,L,c,ssum/nSamp,esum/nSamp,upchar,downchar,visual);
      fflush(stdout);
    }
  }
  fprintf(stdout,"# The average magnetization is %.5lf\n",ssum/nSamp);
  fprintf(stdout,"# The average energy per spin is %.5lf\n",esum/nSamp);
  fprintf(stdout,"# Program ends.\n");
}
