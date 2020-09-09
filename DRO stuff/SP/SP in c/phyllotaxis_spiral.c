/*
 * Copyright (c) 2016 by General Electric Company. All Rights Reserved.
 */

/**
 * \file phyllotaxis_spiral.c *
 * This is the source file for Cartesian spiral
 *
 * @author R. Marc Lebel
 * @since 25.0
 */

/*
 * Comments:
 * Date (dd mmm yyyy) Author (or Initials)
 * Comment
 *
 * 13 May 2016 RML
 *  Initial version.
 */


#ifndef PI
#define PI 3.141592653589793115997963468544
#endif

#ifndef TWOPI
#define TWOPI 6.283185307179586231995926937088
#endif

#ifndef MGAN
#define MGAN 0.618033988749894902525738871191
#endif

/* Fully sampled centre */
#ifndef NFULLC
#define NFULLC 40
#endif

/* Eddy current calibration, averages and encodes */
#ifndef NEDDYR
#define NEDDYR 3
#endif

#ifndef NEDDYE
#define NEDDYE 8
#endif

#ifndef NSTST
#define NSTST 512
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "phyllotaxis_spiral.h"


/* Simple but adequate and reproducible random number generator */
long long myrand(long long seed)
{
    return( (seed*16807LL) % 2147483647LL );
}


int gen_phyllotaxis_spiral(const char *filename, int yres, int zres,
                              int N, int pts, double target_rot, double density)
{

    /* Declare variables */
    FILE *fid;
    const int fibonacci[] = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765};
    int fib_base = 0, Nact, Nactr, arms, armsneeded, ptsr, ovs, shft;
    long long seed=1;
    double rnum;
    double slope_rot, actual_rot = 0.0, cur_rot;
    double diff_rot = 1000000.0;
    int indx, indx2, indx3;
    int ploc;
    int sloc;
    double cur_phs,cur_rad;
    int *ky, *kz;
    double kyt, kzt;
    /*double kys, kzs, sphs, samp;*/
    int nfully = NFULLC;
    int nfullz = NFULLC;

    /* Estimate oversampling factor */
    ovs = ceil(3.0 * (double)((yres+zres)/2) / (double)(pts));
    pts = ovs * pts;

    /* Reduce number of points per arm (will manually add ky=0 kz=0 to each arm */
    ptsr = pts-1;

    /* Compute total rotation using each Fib. number (a bit brute force) */
    for (indx=0; indx<19; indx++) {
        Nact = ptsr*fibonacci[indx] * (N*ovs/(ptsr*fibonacci[indx]) + 1);
        arms = Nact/ptsr;

        slope_rot = fmod(MGAN*(double)(arms),1.0);
        if (slope_rot > 0.5) {
            slope_rot = 1.0 - slope_rot;
        }

        cur_rot = slope_rot * ptsr;
        if ( fabs(cur_rot - target_rot) < diff_rot) {
            diff_rot = fabs(cur_rot - target_rot);
            fib_base = fibonacci[indx];
            actual_rot = cur_rot;
        }
    }

    /* Compute number of total points with chosen Fibonacci number */
    Nact = pts*fib_base * (N*ovs/(pts*fib_base) + 1);
    Nactr = ptsr*fib_base * (N*ovs/(pts*fib_base) +1);
    arms = Nactr/ptsr;
    armsneeded = ceil(2 * N/(pts/ovs+1));

    /* Disply some debug info */
    printf("\nMaking phyllotaxic spiral:\n");
    printf("\tFibonacci number: %d\n",fib_base);
    printf("\tInput points: %d\n",N);
    printf("\tArm oversampling: %d\n",ovs);
    printf("\tTotal points computed: %d\n",Nact);
    printf("\tTotal points computed (reduced): %d\n",Nactr);
    printf("\tPoints per spiral: %d\n",pts);
    printf("\tSpiral arms: %d\n",arms);
    printf("\tRotations per arm: %f\n",actual_rot);
    fflush(stdout);
    
    /* Allocate memory */
    ky = malloc(sizeof(int) * (armsneeded*(pts/ovs+1) + nfully*nfullz));
    kz = malloc(sizeof(int) * (armsneeded*(pts/ovs+1) + nfully*nfullz));


    /* Start phase encode counter */
    sloc = 0;

    /* Define dummy samples */
    for (indx=0; indx<NSTST; indx++) {
        ky[sloc] = yres/2 + 1;
        kz[sloc] = zres/2 + 1;
        sloc++;
    }


    /* Define eddy current estimation steps */
    /* z-encodes */
    for (indx=0; indx<NEDDYR; indx++) {
        for (indx2=0; indx2<NEDDYE; indx2++) {
            for (indx3=0; indx3<NEDDYE; indx3++) {
                /* Donor */
                kz[sloc] = 1;
                ky[sloc] = yres/2 + 1;
                sloc++;

                /* Encode */
                kz[sloc] = zres/2 + 1 - NEDDYE/2 + indx2;
                ky[sloc] = yres/2 + 1 - NEDDYE/2 + indx3;
                sloc++;

                /* Reset */
                kz[sloc] = zres/2 + 1;
                ky[sloc] = yres/2 + 1;
                sloc++;

                /* Donor (baseline) */
                kz[sloc] = zres/2 + 1;
                ky[sloc] = yres/2 + 1;
                sloc++;

                /* Encode */
                kz[sloc] = zres/2 + 1 - NEDDYE/2 + indx2;
                ky[sloc] = yres/2 + 1 - NEDDYE/2 + indx3;
                sloc++;

                /* Reset */
                kz[sloc] = zres/2 + 1;
                ky[sloc] = yres/2 + 1;
                sloc++;
            }
        }
    }
    /* y-encodes */
    for (indx=0; indx<NEDDYR; indx++) {
        for (indx2=0; indx2<NEDDYE; indx2++) {
            for (indx3=0; indx3<NEDDYE; indx3++) {
                /* Donor */
                kz[sloc] = zres/2 + 1;
                ky[sloc] = 1;
                sloc++;

                /* Encode */
                kz[sloc] = zres/2 + 1 - NEDDYE/2 + indx2;
                ky[sloc] = yres/2 + 1 - NEDDYE/2 + indx3;
                sloc++;

                /* Reset */
                kz[sloc] = zres/2 + 1;
                ky[sloc] = yres/2 + 1;
                sloc++;

                /* Donor (baseline) */
                kz[sloc] = zres/2 + 1;
                ky[sloc] = yres/2 + 1;
                sloc++;

                /* Encode */
                kz[sloc] = zres/2 + 1 - NEDDYE/2 + indx2;
                ky[sloc] = yres/2 + 1 - NEDDYE/2 + indx3;
                sloc++;
            
                /* Reset */
                kz[sloc] = zres/2 + 1;
                ky[sloc] = yres/2 + 1;
                sloc++;
            }
        }
    }

    

    /* Define fully sampled centre */
    for (indx = 0; indx<nfully; indx++) {
        kyt = indx  - nfully/2 + 1 + yres/2;
        for (indx2 = 0; indx2<nfullz; indx2++) {
            kzt = indx2 - nfullz/2 + 1 + zres/2;
            ky[sloc] = kyt;
            kz[sloc] = kzt;
            sloc++;
        }
    }

    /* SPIRAL SECTION */
    /* Loop through samples (in time order) */
    for (indx = 0; indx<armsneeded; indx++) { /*Arm loop */
        seed = myrand(seed);
        rnum = (seed%1000) / 1000.0;
        shft = (int)(rnum*ovs);
        for (indx2 = (pts-1); indx2>=0; indx2--) { /* Point loop */
            
            if (indx2 == 0) { /* First point in arm */
                ploc = 0;
            } else {
                ploc = (indx2-1)*arms + indx;/*/ptsr;*/
            }

            /* Compute phase and radius at this location */
            cur_phs = TWOPI * MGAN * ploc;
            cur_rad = sqrt(2.0f) * pow((double)(ploc)/(double)(Nactr),density);

            /* Compute ky and kz locations */
            kyt = cur_rad * cos(cur_phs);
            kzt = cur_rad * sin(cur_phs);

            /* Convert to lookup table format (from 1 to y/zres) */
            kyt = floor((double)(yres)*(kyt/2.0 + 0.5)+1.0);
            kzt = floor((double)(zres)*(kzt/2.0 + 0.5)+1.0);

            /* Store in output */
            if ( indx2 ==0 || 
               ( (((indx2+shft) % ovs) == 0) &&
                kyt>=1.0 && kyt<=yres &&
                kzt>=1.0 && kzt<=zres ) )
            {
                ky[sloc] = (int)(kyt);
                kz[sloc] = (int)(kzt);
                
                /* Increment array index */
                sloc++;
            }
        }
    }


    /* write text file */
    fid = fopen(filename,"w");
    for (indx=0; indx<N; indx++)
    {
        fprintf(fid, "%d %d\n", ky[indx], kz[indx]);
    }
    fclose(fid);


    /* Free memory */
    free(ky);
    free(kz);

    return 0;
}


int main() {
    
    /* Hardcode for demo purposes */
    char filename[]="pe_table.txt";
    int yres=240;
    int zres=120;
    int N=240*120*4; /* total points to generate */
    int pts=45; /* points per arm */
    double target_rot=0.333; /* some rotation, not too much */
    double density=0.65; /* 1.0=uniform arm; 0.5=uniform k-space */
    int status;
    
    status = gen_phyllotaxis_spiral(filename, yres, zres, N, pts, target_rot, density);
}
