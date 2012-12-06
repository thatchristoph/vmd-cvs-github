/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2011 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: QuickSurf.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.63 $	$Date: 2012/03/20 15:41:03 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Fast gaussian surface representation
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "QuickSurf.h"
#include "CUDAQuickSurf.h"
#include "Measure.h"
#include "Inform.h"
#include "utilities.h"
#include "WKFUtils.h"
#include "VolumetricData.h"

#include "VMDDisplayList.h"
#include "Displayable.h"
#include "DispCmds.h"

#define MIN(X,Y) (((X)<(Y))? (X) : (Y))
#define MAX(X,Y) (((X)>(Y))? (X) : (Y))


/*
 * David J. Hardy
 * 12 Dec 2008
 *
 * aexpfnx() - Approximate expf() for negative x.
 *
 * Assumes that x <= 0.
 *
 * Assumes IEEE format for single precision float, specifically:
 * 1 sign bit, 8 exponent bits biased by 127, and 23 mantissa bits.
 *
 * Interpolates exp() on interval (-1/log2(e), 0], then shifts it by
 * multiplication of a fast calculation for 2^(-N).  The interpolation
 * uses a linear blending of 3rd degree Taylor polynomials at the end
 * points, so the approximation is once differentiable.
 *
 * The error is small (max relative error per interval is calculated
 * to be 0.131%, with a max absolute error of -0.000716).
 *
 * The cutoff is chosen so as to speed up the computation by early
 * exit from function, with the value chosen to give less than the
 * the max absolute error.  Use of a cutoff is unnecessary, except
 * for needing to shift smallest floating point numbers to zero,
 * i.e. you could remove cutoff and replace by:
 *
 * #define MINXNZ  -88.0296919311130  // -127 * log(2)
 *
 *   if (x < MINXNZ) return 0.f;
 *
 * Use of a cutoff causes a discontinuity which can be eliminated
 * through the use of a switching function.
 *
 * We can obtain arbitrarily smooth approximation by taking k+1 nodes on
 * the interval and weighting their respective Taylor polynomials by the
 * kth order Lagrange interpolant through those nodes.  The wiggle in the
 * polynomial interpolation due to equidistant nodes (Runge's phenomenon)
 * can be reduced by using Chebyshev nodes.
 */

#define MLOG2EF    -1.44269504088896f

/*
 * Interpolating coefficients for linear blending of the
 * 3rd degree Taylor expansion of 2^x about 0 and -1.
 */
#define SCEXP0     1.0000000000000000f
#define SCEXP1     0.6987082824680118f
#define SCEXP2     0.2633174272827404f
#define SCEXP3     0.0923611991471395f
#define SCEXP4     0.0277520543324108f

/* for single precision float */
#define EXPOBIAS   127
#define EXPOSHIFT   23

/* cutoff is optional, but can help avoid unnecessary work */
#define ACUTOFF    -10

typedef union flint_t {
  float f;
  int n;
} flint;

static float aexpfnx(float x) {
  /* assume x <= 0 */
  float mb;
  int mbflr;
  float d;
  float sy;
  flint scalfac;

  if (x < ACUTOFF) return 0.f;

  mb = x * MLOG2EF;    /* change base to 2, mb >= 0 */
  mbflr = (int) mb;    /* get int part, floor() */
  d = mbflr - mb;      /* remaining exponent, -1 < d <= 0 */
  sy = SCEXP0 + d*(SCEXP1 + d*(SCEXP2 + d*(SCEXP3 + d*SCEXP4)));
                       /* approx with linear blend of Taylor polys */
  scalfac.n = (EXPOBIAS - mbflr) << EXPOSHIFT;  /* 2^(-mbflr) */
  return (sy * scalfac.f);  /* scaled approx */
}


static void vmd_gaussdensity(int natoms, const float *xyzr,
                             const float *colors,
                             float *densitymap, float *voltexmap, 
                             const int *numvoxels, 
                             float radscale, float gridspacing, 
                             float isovalue, float gausslim) {
  int i, x, y, z;
  int maxvoxel[3];
  maxvoxel[0] = numvoxels[0]-1; 
  maxvoxel[1] = numvoxels[1]-1; 
  maxvoxel[2] = numvoxels[2]-1; 
  const float invgridspacing = 1.0f / gridspacing;

  // compute colors only if necessary, since they are costly
  if (voltexmap != NULL) {
    float invisovalue = 1.0f / isovalue;
    // compute both density map and floating point color texture map
    for (i=0; i<natoms; i++) {
#if !defined(ARCH_BLUEWATERS)
      if ((i & 0x3fff) == 0) {
        printf("."); 
        fflush(stdout);
      }
#endif

      int ind = i*4;
      float scaledrad = xyzr[ind + 3] * radscale;
      float arinv = 1.0f/(2.0f*scaledrad*scaledrad);
      float radlim = gausslim * scaledrad;
      float radlim2 = radlim * radlim;

      float tmp;
      radlim *= invgridspacing;
      tmp = xyzr[ind  ] * invgridspacing;
      int xmin = MAX((int) (tmp - radlim), 0);
      int xmax = MIN((int) (tmp + radlim), maxvoxel[0]);
      tmp = xyzr[ind+1] * invgridspacing;
      int ymin = MAX((int) (tmp - radlim), 0);
      int ymax = MIN((int) (tmp + radlim), maxvoxel[1]);
      tmp = xyzr[ind+2] * invgridspacing;
      int zmin = MAX((int) (tmp - radlim), 0);
      int zmax = MIN((int) (tmp + radlim), maxvoxel[2]);

      float dz = zmin*gridspacing - xyzr[ind+2];
      for (z=zmin; z<=zmax; z++,dz+=gridspacing) {
        float dy = ymin*gridspacing - xyzr[ind+1];
        for (y=ymin; y<=ymax; y++,dy+=gridspacing) {
          float dy2dz2 = dy*dy + dz*dz;

          // early-exit when outside the cutoff radius in the Y-Z plane
          if (dy2dz2 >= radlim2) 
            continue;

          int addr = z * numvoxels[0] * numvoxels[1] + y * numvoxels[0];
          float dx = xmin*gridspacing - xyzr[ind];
          for (x=xmin; x<=xmax; x++,dx+=gridspacing) {
            float r2 = dx*dx + dy2dz2;
            float expval = -r2 * arinv;
#if 0
            // use the math library exponential routine
            float density = exp(expval);
#else
            // use our (much faster) fast exponential approximation
            float density = aexpfnx(expval);
#endif

            // accumulate density value to density map
            densitymap[addr + x] += density;

            // Accumulate density-weighted color to texture map.
            // Pre-multiply colors by the inverse isovalue we will extract   
            // the surface on, to cause the final color to be normalized.
            density *= invisovalue;
            int caddr = (addr + x) * 3;

            // color by atom colors
            voltexmap[caddr    ] += density * colors[ind    ];
            voltexmap[caddr + 1] += density * colors[ind + 1];
            voltexmap[caddr + 2] += density * colors[ind + 2];
          }
        }
      }
    }
  } else {
    // compute density map only
    for (i=0; i<natoms; i++) {
#if !defined(ARCH_BLUEWATERS)
      if ((i & 0x3fff) == 0) {
        printf("."); 
        fflush(stdout);
      }
#endif

      int ind = i*4;
      float scaledrad = xyzr[ind + 3] * radscale;
      float arinv = 1.0f/(2.0f*scaledrad*scaledrad);
      float radlim = gausslim * scaledrad;
      float radlim2 = radlim * radlim;

      float tmp;
      radlim *= invgridspacing;
      tmp = xyzr[ind  ] * invgridspacing;
      int xmin = MAX((int) (tmp - radlim), 0);
      int xmax = MIN((int) (tmp + radlim), maxvoxel[0]);
      tmp = xyzr[ind+1] * invgridspacing;
      int ymin = MAX((int) (tmp - radlim), 0);
      int ymax = MIN((int) (tmp + radlim), maxvoxel[1]);
      tmp = xyzr[ind+2] * invgridspacing;
      int zmin = MAX((int) (tmp - radlim), 0);
      int zmax = MIN((int) (tmp + radlim), maxvoxel[2]);

      float dz = zmin*gridspacing - xyzr[ind+2];
      for (z=zmin; z<=zmax; z++,dz+=gridspacing) {
        float dy = ymin*gridspacing - xyzr[ind+1];
        for (y=ymin; y<=ymax; y++,dy+=gridspacing) {
          float dy2dz2 = dy*dy + dz*dz;

          // early-exit when outside the cutoff radius in the Y-Z plane
          if (dy2dz2 >= radlim2) 
            continue;

          int addr = z * numvoxels[0] * numvoxels[1] + y * numvoxels[0];
          float dx = xmin*gridspacing - xyzr[ind];
          for (x=xmin; x<=xmax; x++,dx+=gridspacing) {
            float r2 = dx*dx + dy2dz2;
            float expval = -r2 * arinv;
#if 0
            // use the math library exponential routine
            float density = exp(expval);
#else
            // use our (much faster) fast exponential approximation
            float density = aexpfnx(expval);
#endif
            densitymap[addr + x] += density;
          }
        }
      }
    }
  }
}


static void vmd_gaussdensity_opt(int natoms, const float *xyzr,
                                 const float *colors,
                                 float *densitymap, float *voltexmap, 
                                 const int *numvoxels, 
                                 float radscale, float gridspacing, 
                                 float isovalue, float gausslim) {
  int i, x, y, z;
  int maxvoxel[3];
  maxvoxel[0] = numvoxels[0]-1; 
  maxvoxel[1] = numvoxels[1]-1; 
  maxvoxel[2] = numvoxels[2]-1; 
  const float invgridspacing = 1.0f / gridspacing;

  // compute colors only if necessary, since they are costly
  if (voltexmap != NULL) {
    float invisovalue = 1.0f / isovalue;
    // compute both density map and floating point color texture map
    for (i=0; i<natoms; i++) {
#if !defined(ARCH_BLUEWATERS)
      if ((i & 0x3fff) == 0) {
        printf("."); 
        fflush(stdout);
      }
#endif

      int ind = i*4;
      float scaledrad = xyzr[ind + 3] * radscale;
      // negate, precompute reciprocal, and change to base 2 from the outset
      float arinv = -(1.0f/(2.0f*scaledrad*scaledrad)) * MLOG2EF;
      float radlim = gausslim * scaledrad;
      float radlim2 = radlim * radlim;

      float tmp;
      radlim *= invgridspacing;
      tmp = xyzr[ind  ] * invgridspacing;
      int xmin = MAX((int) (tmp - radlim), 0);
      int xmax = MIN((int) (tmp + radlim), maxvoxel[0]);
      tmp = xyzr[ind+1] * invgridspacing;
      int ymin = MAX((int) (tmp - radlim), 0);
      int ymax = MIN((int) (tmp + radlim), maxvoxel[1]);
      tmp = xyzr[ind+2] * invgridspacing;
      int zmin = MAX((int) (tmp - radlim), 0);
      int zmax = MIN((int) (tmp + radlim), maxvoxel[2]);

      float dz = zmin*gridspacing - xyzr[ind+2];
      for (z=zmin; z<=zmax; z++,dz+=gridspacing) {
        float dy = ymin*gridspacing - xyzr[ind+1];
        for (y=ymin; y<=ymax; y++,dy+=gridspacing) {
          float dy2dz2 = dy*dy + dz*dz;

          // early-exit when outside the cutoff radius in the Y-Z plane
          if (dy2dz2 >= radlim2) 
            continue;

          int addr = z * numvoxels[0] * numvoxels[1] + y * numvoxels[0];
          float dx = xmin*gridspacing - xyzr[ind];
          for (x=xmin; x<=xmax; x++,dx+=gridspacing) {
            float r2 = dx*dx + dy2dz2;

            // use our (much faster) fully inlined exponential approximation
            float mb = r2 * arinv;         /* already negated and in base 2 */
            int mbflr = (int) mb;          /* get int part, floor() */
            float d = mbflr - mb;          /* remaining exponent, -1 < d <= 0 */

            /* approx with linear blend of Taylor polys */
            float sy = SCEXP0 + d*(SCEXP1 + d*(SCEXP2 + d*(SCEXP3 + d*SCEXP4)));

            /* 2^(-mbflr) */
            flint scalfac;
            scalfac.n = (EXPOBIAS - mbflr) << EXPOSHIFT;  

            // XXX assume we are never beyond the cutoff value in this loop
            float density = (sy * scalfac.f);

            // accumulate density value to density map
            densitymap[addr + x] += density;

            // Accumulate density-weighted color to texture map.
            // Pre-multiply colors by the inverse isovalue we will extract   
            // the surface on, to cause the final color to be normalized.
            density *= invisovalue;
            int caddr = (addr + x) * 3;

            // color by atom colors
            voltexmap[caddr    ] += density * colors[ind    ];
            voltexmap[caddr + 1] += density * colors[ind + 1];
            voltexmap[caddr + 2] += density * colors[ind + 2];
          }
        }
      }
    }
  } else {
    // compute density map only
    for (i=0; i<natoms; i++) {
#if !defined(ARCH_BLUEWATERS)
      if ((i & 0x3fff) == 0) {
        printf("."); 
        fflush(stdout);
      }
#endif

      int ind = i*4;
      float scaledrad = xyzr[ind+3] * radscale;

      // negate, precompute reciprocal, and change to base 2 from the outset
      float arinv = -(1.0f/(2.0f*scaledrad*scaledrad)) * MLOG2EF;
      float radlim = gausslim * scaledrad;
      float radlim2 = radlim * radlim;

      float tmp;
      radlim *= invgridspacing;
      tmp = xyzr[ind  ] * invgridspacing;
      int xmin = MAX((int) (tmp - radlim), 0);
      int xmax = MIN((int) (tmp + radlim), maxvoxel[0]);
      tmp = xyzr[ind+1] * invgridspacing;
      int ymin = MAX((int) (tmp - radlim), 0);
      int ymax = MIN((int) (tmp + radlim), maxvoxel[1]);
      tmp = xyzr[ind+2] * invgridspacing;
      int zmin = MAX((int) (tmp - radlim), 0);
      int zmax = MIN((int) (tmp + radlim), maxvoxel[2]);

      float dz = zmin*gridspacing - xyzr[ind+2];
      for (z=zmin; z<=zmax; z++,dz+=gridspacing) {
        float dy = ymin*gridspacing - xyzr[ind+1];
        for (y=ymin; y<=ymax; y++,dy+=gridspacing) {
          float dy2dz2 = dy*dy + dz*dz;

          // early-exit when outside the cutoff radius in the Y-Z plane
          if (dy2dz2 >= radlim2) 
            continue;

          int addr = z * numvoxels[0] * numvoxels[1] + y * numvoxels[0];
          float dx = xmin*gridspacing - xyzr[ind];
          for (x=xmin; x<=xmax; x++,dx+=gridspacing) {
            float r2 = dx*dx + dy2dz2;

            // use our (much faster) fully inlined exponential approximation
            float mb = r2 * arinv;         /* already negated and in base 2 */
            int mbflr = (int) mb;          /* get int part, floor() */
            float d = mbflr - mb;          /* remaining exponent, -1 < d <= 0 */

            /* approx with linear blend of Taylor polys */
            float sy = SCEXP0 + d*(SCEXP1 + d*(SCEXP2 + d*(SCEXP3 + d*SCEXP4)));

            /* 2^(-mbflr) */
            flint scalfac;
            scalfac.n = (EXPOBIAS - mbflr) << EXPOSHIFT;  

            // XXX assume we are never beyond the cutoff value in this loop
            float density = (sy * scalfac.f);

            densitymap[addr + x] += density;
          }
        }
      }
    }
  }
}


typedef struct {
  int natoms;
  float radscale;
  float gridspacing;
  float isovalue;
  float gausslim;
  const int *numvoxels;
  const float *xyzr; 
  const float *colors;
  float **thrdensitymaps;
  float **thrvoltexmaps;
} densitythrparms;


static void * densitythread(void *voidparms) {
  wkf_tasktile_t tile;
  densitythrparms *parms = NULL;
  int threadid;

  wkf_threadlaunch_getid(voidparms, &threadid, NULL);
  wkf_threadlaunch_getdata(voidparms, (void **) &parms);

  while (wkf_threadlaunch_next_tile(voidparms, 16384, &tile) != WKF_SCHED_DONE) {
    int natoms = tile.end-tile.start;
    vmd_gaussdensity_opt(natoms, 
                         &parms->xyzr[4*tile.start], 
                         (parms->thrvoltexmaps[0]!=NULL) ? &parms->colors[4*tile.start] : NULL,
                         parms->thrdensitymaps[threadid], 
                         parms->thrvoltexmaps[threadid], 
                         parms->numvoxels, 
                         parms->radscale, 
                         parms->gridspacing, 
                         parms->isovalue, 
                         parms->gausslim);
  }

  return NULL;
}


static void * reductionthread(void *voidparms) {
  wkf_tasktile_t tile;
  densitythrparms *parms = NULL;
  int threadid, numthreads;

  wkf_threadlaunch_getid(voidparms, &threadid, &numthreads);
  wkf_threadlaunch_getdata(voidparms, (void **) &parms);

  while (wkf_threadlaunch_next_tile(voidparms, 16384, &tile) != WKF_SCHED_DONE) {
    // do a reduction over each of the individual density grids
    int i, x;
    for (x=tile.start; x<tile.end; x++) {
      float tmp = 0.0f;
      for (i=1; i<numthreads; i++) {
        tmp += parms->thrdensitymaps[i][x];
      }
      parms->thrdensitymaps[0][x] += tmp;
    }

    // do a reduction over each of the individual texture grids
    if (parms->thrvoltexmaps[0] != NULL) {
      for (x=tile.start*3; x<tile.end*3; x++) {
        float tmp = 0.0f;
        for (i=1; i<numthreads; i++) {
          tmp += parms->thrvoltexmaps[i][x];
        }
        parms->thrvoltexmaps[0][x] += tmp;
      }
    }
  }

  return NULL;
}


static int vmd_gaussdensity_threaded(int natoms, const float *xyzr,
                                     const float *colors,
                                     float *densitymap, float *voltexmap, 
                                     const int *numvoxels, 
                                     float radscale, float gridspacing, 
                                     float isovalue, float gausslim) {
  densitythrparms parms;
  memset(&parms, 0, sizeof(parms));

  parms.natoms = natoms;
  parms.radscale = radscale;
  parms.gridspacing = gridspacing;
  parms.isovalue = isovalue;
  parms.gausslim = gausslim;
  parms.numvoxels = numvoxels;
  parms.xyzr = xyzr;
  parms.colors = colors;
 
  int physprocs = wkf_thread_numprocessors();
  int maxprocs = physprocs;

  // We can productively use only a few cores per socket due to the
  // limited memory bandwidth per socket. Also, hyperthreading
  // actually hurts performance.  These two considerations combined
  // with the linear increase in memory use prevent us from using large
  // numbers of cores with this simple approach, so if we've got more 
  // than 8 CPU cores, we'll iteratively cutting the core count in 
  // half until we're under 8 cores.
  while (maxprocs > 8) 
    maxprocs /= 2;

  // Limit the number of CPU cores used so we don't run the 
  // machine out of memory during surface computation.
  // For now we'll set a practical maximum memory use limit to 
  // 2GB total for all cores. 
  long volsz = numvoxels[0] * numvoxels[1] * numvoxels[2];
  long volmemsz = sizeof(float) * volsz;
  long volmemszkb = volmemsz / 1024;
  long volmemtexszkb = (volmemszkb + (voltexmap != NULL) ? 3*volmemszkb : 0);
#if defined(ARCH_BLUEWATERS)
  while ((volmemtexszkb * maxprocs) > (4 * 1024 * 1024))
    maxprocs /= 2;
#else
  while ((volmemtexszkb * maxprocs) > (2 * 1024 * 1024))
    maxprocs /= 2;
#endif

  if (maxprocs < 1) 
    maxprocs = 1;

  // Loop over number of physical processors and try to create 
  // per-thread volumetric maps for each of them.
  parms.thrdensitymaps = (float **) calloc(1,maxprocs * sizeof(float *));
  parms.thrvoltexmaps = (float **) calloc(1, maxprocs * sizeof(float *));

  // first thread is already ready to go
  parms.thrdensitymaps[0] = densitymap;
  parms.thrvoltexmaps[0] = voltexmap;

  int i;
  int numprocs = maxprocs; // ever the optimist
  for (i=1; i<maxprocs; i++) {
    parms.thrdensitymaps[i] = (float *) calloc(1, volmemsz);
    if (parms.thrdensitymaps[i] == NULL) {
      numprocs = i;
      break;
    }
    if (voltexmap != NULL) {
      parms.thrvoltexmaps[i] = (float *) calloc(1, 3 * volmemsz);
      if (parms.thrvoltexmaps[i] == NULL) {
        free(parms.thrdensitymaps[i]);
        parms.thrdensitymaps[i] = NULL;
        numprocs = i;
        break;
      }
    }
  }

  // launch independent thread calculations
  wkf_tasktile_t tile;
  tile.start = 0;
  tile.end = natoms;
  wkf_threadlaunch(numprocs, &parms, densitythread, &tile);

  // do a parallel reduction of the resulting density maps
  tile.start = 0;
  tile.end = volsz;
  wkf_threadlaunch(numprocs, &parms, reductionthread, &tile);

  // free work area
  for (i=1; i<maxprocs; i++) {
    if (parms.thrdensitymaps[i] != NULL)
      free(parms.thrdensitymaps[i]);

    if (parms.thrvoltexmaps[i] != NULL)
      free(parms.thrvoltexmaps[i]);
  }
  free(parms.thrdensitymaps);
  free(parms.thrvoltexmaps);

  return 0;
}

QuickSurf::QuickSurf() {
  volmap = NULL;
  voltexmap = NULL;
  s.clear();
  isovalue = 0.5f;

  numvoxels[0] = 128;
  numvoxels[1] = 128;
  numvoxels[2] = 128;

  origin[0] = 0.0f;
  origin[1] = 0.0f;
  origin[2] = 0.0f;

  xaxis[0] = 1.0f;
  xaxis[1] = 0.0f;
  xaxis[2] = 0.0f;

  yaxis[0] = 0.0f;
  yaxis[1] = 1.0f;
  yaxis[2] = 0.0f;

  zaxis[0] = 0.0f;
  zaxis[1] = 0.0f;
  zaxis[2] = 1.0f;
   
  cudaqs = NULL;
#if defined(VMDCUDA)
  if (!getenv("VMDNOCUDA")) {
    cudaqs = new CUDAQuickSurf();
  }
#endif

  timer = wkf_timer_create();
}

int QuickSurf::calc_surf(AtomSel *atomSel, DrawMolecule *mol,
                         const float *atompos, const float *atomradii,
                         int quality, float radscale, float gridspacing,
                         float isoval, const int *colidx, const float *cmap,
                         VMDDisplayList *cmdList) {
  wkf_timer_start(timer);
  int colorperatom = (colidx != NULL && cmap != NULL);
  int usebeads=0;

  // clean up any existing CPU arrays before going any further...
  if (voltexmap != NULL)
    free(voltexmap);
  voltexmap = NULL;

  ResizeArray<float> beadpos(64 + (3 * atomSel->selected) / 20);
  ResizeArray<float> beadradii(64 + (3 * atomSel->selected) / 20);
  ResizeArray<float> beadcolors(64 + (3 * atomSel->selected) / 20);

  if (getenv("VMDQUICKSURFBEADS")) {
    usebeads=1;
#if !defined(ARCH_BLUEWATERS)
    printf("QuickSurf using residue beads representation...\n");
#endif
  }

  int numbeads = 0;
  if (usebeads) {
    int i, resid, numres;

    // draw a bead for each residue
    numres = mol->residueList.num();
    for (resid=0; resid<numres; resid++) {
      float com[3] = {0.0, 0.0, 0.0};
      const ResizeArray<int> &atoms = mol->residueList[resid]->atoms;
      int numatoms = atoms.num();
      int oncount = 0;
   
      // find COM for residue
      for (i=0; i<numatoms; i++) {
        int idx = atoms[i];
        if (atomSel->on[idx]) {
          oncount++;
          vec_add(com, com, atompos + 3*idx);
        }
      }

      if (oncount < 1)
        continue; // exit if there weren't any atoms

      vec_scale(com, 1.0f / (float) oncount, com);

      // find radius of bounding sphere and save last atom index for color
      int atomcolorindex=0; // initialize, to please compilers
      float boundradsq = 0.0f;
      for (i=0; i<numatoms; i++) {
        int idx = atoms[i];
        if (atomSel->on[idx]) {
          float tmpdist[3];
          atomcolorindex = idx;
          vec_sub(tmpdist, com, atompos + 3*idx);
          float distsq = dot_prod(tmpdist, tmpdist);
          if (distsq > boundradsq) {
            boundradsq = distsq;
          }
        }
      }
      beadpos.append(com[0]);
      beadpos.append(com[1]);
      beadpos.append(com[2]);

      beadradii.append(sqrtf(boundradsq) + 1.0f);

      if (colorperatom) {
        const float *cp = &cmap[colidx[atomcolorindex] * 3];
        beadcolors.append(cp[0]);
        beadcolors.append(cp[1]);
        beadcolors.append(cp[2]);
      }

      // XXX still need to add pick points...
    }

    numbeads = beadpos.num() / 3;
  }

  // initialize class variables
  isovalue=isoval;

  // If no volumetric texture will be computed we will use the cmap
  // parameter to pass in the solid color to be applied to all vertices
  vec_copy(solidcolor, cmap);

  // compute min/max atom radius, build list of selected atom radii,
  // and compute bounding box for the selected atoms
  float minx, miny, minz, maxx, maxy, maxz;
  float minrad, maxrad;
  int i;
  if (usebeads) {
    minx = maxx = beadpos[0];
    miny = maxy = beadpos[1];
    minz = maxz = beadpos[2];
    minrad = maxrad = beadradii[0];
    for (i=0; i<numbeads; i++) {
      int ind = i * 3;
      float tmpx = beadpos[ind  ];
      float tmpy = beadpos[ind+1];
      float tmpz = beadpos[ind+2];

      minx = (tmpx < minx) ? tmpx : minx;
      maxx = (tmpx > maxx) ? tmpx : maxx;

      miny = (tmpy < miny) ? tmpy : miny;
      maxy = (tmpy > maxy) ? tmpy : maxy;

      minz = (tmpz < minz) ? tmpz : minz;
      maxz = (tmpz > maxz) ? tmpz : maxz;
 
      // we always have to compute the rmin/rmax for beads
      // since these radii are defined on-the-fly
      float r = beadradii[i];
      minrad = (r < minrad) ? r : minrad;
      maxrad = (r > maxrad) ? r : maxrad;
    }
  } else {
    minx = maxx = atompos[atomSel->firstsel*3  ];
    miny = maxy = atompos[atomSel->firstsel*3+1];
    minz = maxz = atompos[atomSel->firstsel*3+2];

    // Query min/max atom radii for the entire molecule
    mol->get_radii_minmax(minrad, maxrad);

    // We only compute rmin/rmax for the actual group of selected atoms if 
    // (rmax/rmin > 2.5) for the whole molecule, otherwise it's a small 
    // enough range that we don't care since it won't hurt our performance. 
    if (minrad <= 0.001 || maxrad/minrad > 2.5) {
      minrad = maxrad = atomradii[atomSel->firstsel];
      for (i=atomSel->firstsel; i<=atomSel->lastsel; i++) {
        if (atomSel->on[i]) {
          int ind = i * 3;
          float tmpx = atompos[ind  ];
          float tmpy = atompos[ind+1];
          float tmpz = atompos[ind+2];

          minx = (tmpx < minx) ? tmpx : minx;
          maxx = (tmpx > maxx) ? tmpx : maxx;

          miny = (tmpy < miny) ? tmpy : miny;
          maxy = (tmpy > maxy) ? tmpy : maxy;

          minz = (tmpz < minz) ? tmpz : minz;
          maxz = (tmpz > maxz) ? tmpz : maxz;
  
          float r = atomradii[i];
          minrad = (r < minrad) ? r : minrad;
          maxrad = (r > maxrad) ? r : maxrad;
        }
      }
    } else {
      for (i=atomSel->firstsel; i<=atomSel->lastsel; i++) {
        if (atomSel->on[i]) {
          int ind = i * 3;
          float tmpx = atompos[ind  ];
          float tmpy = atompos[ind+1];
          float tmpz = atompos[ind+2];

          minx = (tmpx < minx) ? tmpx : minx;
          maxx = (tmpx > maxx) ? tmpx : maxx;

          miny = (tmpy < miny) ? tmpy : miny;
          maxy = (tmpy > maxy) ? tmpy : maxy;

          minz = (tmpz < minz) ? tmpz : minz;
          maxz = (tmpz > maxz) ? tmpz : maxz;
        }
      }
    }
  }

  float mincoord[3], maxcoord[3];
  mincoord[0] = minx;
  mincoord[1] = miny;
  mincoord[2] = minz;
  maxcoord[0] = maxx;
  maxcoord[1] = maxy;
  maxcoord[2] = maxz;

  // crude estimate of the grid padding we require to prevent the
  // resulting isosurface from being clipped
  float gridpadding = radscale * maxrad * 1.70f;
  float padrad = gridpadding;
  padrad = 0.65f * sqrtf(4.0f/3.0f*((float) VMD_PI)*padrad*padrad*padrad);
  gridpadding = MAX(gridpadding, padrad);

  // Handle coarse-grained structures and whole-cell models
  // XXX The switch at 4.0A from an assumed all-atom scale structure to 
  //     CG or cell models is a simple heuristic at a somewhat arbitrary 
  //     threshold value.  
  //     For all-atom models the units shown in the GUI are in Angstroms
  //     and are absolute, but for CG or cell models the units in the GUI 
  //     are relative to the atom with the minimum radius.
  //     This code doesn't do anything to handle structures with a minrad 
  //     of zero, where perhaps only one particle has an unset radius.
  if (minrad > 4.0f) {
    gridspacing *= minrad;
  }

#if !defined(ARCH_BLUEWATERS)
  printf("QuickSurf: R*%.1f, I=%.1f, H=%.1f Pad: %.1f minR: %.1f maxR: %.1f)\n",
         radscale, isovalue, gridspacing, gridpadding, minrad, maxrad);
#endif

  mincoord[0] -= gridpadding;
  mincoord[1] -= gridpadding;
  mincoord[2] -= gridpadding;
  maxcoord[0] += gridpadding;
  maxcoord[1] += gridpadding;
  maxcoord[2] += gridpadding;

  // compute the real grid dimensions from the selected atoms
  xaxis[0] = maxcoord[0]-mincoord[0];
  yaxis[1] = maxcoord[1]-mincoord[1];
  zaxis[2] = maxcoord[2]-mincoord[2];
  numvoxels[0] = (int) ceil(xaxis[0] / gridspacing);
  numvoxels[1] = (int) ceil(yaxis[1] / gridspacing);
  numvoxels[2] = (int) ceil(zaxis[2] / gridspacing);

  // recalc the grid dimensions from rounded/padded voxel counts
  xaxis[0] = (numvoxels[0]-1) * gridspacing;
  yaxis[1] = (numvoxels[1]-1) * gridspacing;
  zaxis[2] = (numvoxels[2]-1) * gridspacing;
  maxcoord[0] = mincoord[0] + xaxis[0];
  maxcoord[1] = mincoord[1] + yaxis[1];
  maxcoord[2] = mincoord[2] + zaxis[2];

#if !defined(ARCH_BLUEWATERS)
  printf("  GridSZ: (%4d %4d %4d)  BBox: (%.1f %.1f %.1f)->(%.1f %.1f %.1f)\n",
         numvoxels[0], numvoxels[1], numvoxels[2],
         mincoord[0], mincoord[1], mincoord[2],
         maxcoord[0], maxcoord[1], maxcoord[2]);
#endif

  vec_copy(origin, mincoord);

  // build compacted lists of bead coordinates, radii, and colors
  float *xyzr = NULL;
  float *colors = NULL;
  if (usebeads) { 
    int ind =0;
    int ind4=0; 
    xyzr = (float *) malloc(numbeads * sizeof(float) * 4);
    if (colorperatom) {
      colors = (float *) malloc(numbeads * sizeof(float) * 4);

      // build compacted lists of bead coordinates, radii, and colors
      for (i=0; i<numbeads; i++) {
        const float *fp = &beadpos[0] + ind;
        xyzr[ind4    ] = fp[0]-origin[0];
        xyzr[ind4 + 1] = fp[1]-origin[1];
        xyzr[ind4 + 2] = fp[2]-origin[2];
        xyzr[ind4 + 3] = beadradii[i];
 
        const float *cp = &beadcolors[0] + ind;
        colors[ind4    ] = cp[0];
        colors[ind4 + 1] = cp[1];
        colors[ind4 + 2] = cp[2];
        colors[ind4 + 3] = 1.0f;
        ind4 += 4;
        ind += 3;
      }
    } else {
      // build compacted lists of bead coordinates and radii only
      for (i=0; i<numbeads; i++) {
        const float *fp = &beadpos[0] + ind;
        xyzr[ind4    ] = fp[0]-origin[0];
        xyzr[ind4 + 1] = fp[1]-origin[1];
        xyzr[ind4 + 2] = fp[2]-origin[2];
        xyzr[ind4 + 3] = beadradii[i];
        ind4 += 4;
        ind += 3;
      }
    }
  } else {
    int ind = atomSel->firstsel * 3;
    int ind4=0; 
    xyzr = (float *) malloc(atomSel->selected * sizeof(float) * 4);
    if (colorperatom) {
      colors = (float *) malloc(atomSel->selected * sizeof(float) * 4);

      // build compacted lists of atom coordinates, radii, and colors
      for (i=atomSel->firstsel; i <= atomSel->lastsel; i++) {
        if (atomSel->on[i]) {
          const float *fp = atompos + ind;
          xyzr[ind4    ] = fp[0]-origin[0];
          xyzr[ind4 + 1] = fp[1]-origin[1];
          xyzr[ind4 + 2] = fp[2]-origin[2];
          xyzr[ind4 + 3] = atomradii[i];
 
          const float *cp = &cmap[colidx[i] * 3];
          colors[ind4    ] = cp[0];
          colors[ind4 + 1] = cp[1];
          colors[ind4 + 2] = cp[2];
          colors[ind4 + 3] = 1.0f;
          ind4 += 4;
        }
        ind += 3;
      }
    } else {
      // build compacted lists of atom coordinates and radii only
      for (i=atomSel->firstsel; i <= atomSel->lastsel; i++) {
        if (atomSel->on[i]) {
          const float *fp = atompos + ind;
          xyzr[ind4    ] = fp[0]-origin[0];
          xyzr[ind4 + 1] = fp[1]-origin[1];
          xyzr[ind4 + 2] = fp[2]-origin[2];
          xyzr[ind4 + 3] = atomradii[i];
          ind4 += 4;
        }
        ind += 3;
      }
    }
  }

  // set gaussian window size based on user-specified quality parameter
  float gausslim = 2.0f;
  switch (quality) {
    case 3: gausslim = 4.0f; break; // max quality

    case 2: gausslim = 3.0f; break; // high quality

    case 1: gausslim = 2.5f; break; // medium quality

    case 0: 
    default: gausslim = 2.0f; // low quality
      break;
  }

  pretime = wkf_timer_timenow(timer);

#if defined(VMDCUDA)
  if (!getenv("VMDNOCUDA")) {
    // compute both density map and floating point color texture map
    int pcount = (usebeads) ? numbeads : atomSel->selected; 
    int rc = cudaqs->calc_surf(pcount, &xyzr[0],
                               (colorperatom) ? &colors[0] : &cmap[0],
                               colorperatom, origin, numvoxels, maxrad,
                               radscale, gridspacing, isovalue, gausslim,
                               cmdList);

    if (rc == 0) {
      free(xyzr);
      if (colors)
        free(colors);
  
      voltime = wkf_timer_timenow(timer);
      return 0;
    }
  }
#endif

#if !defined(ARCH_BLUEWATERS)
  printf("  Computing density map grid on CPUs ");
#endif

  long volsz = numvoxels[0] * numvoxels[1] * numvoxels[2];
  volmap = new float[volsz];
  if (colidx != NULL && cmap != NULL) {
    voltexmap = (float*) calloc(1, 3 * sizeof(float) * numvoxels[0] * numvoxels[1] * numvoxels[2]);
  }

  fflush(stdout);
  memset(volmap, 0, sizeof(float) * volsz);
  if ((volsz * atomSel->selected) > 20000000) {
    vmd_gaussdensity_threaded(atomSel->selected, &xyzr[0],
                              (voltexmap!=NULL) ? &colors[0] : NULL,
                              volmap, voltexmap, numvoxels, radscale, 
                              gridspacing, isovalue, gausslim);
  } else {
    vmd_gaussdensity_opt(atomSel->selected, &xyzr[0],
                         (voltexmap!=NULL) ? &colors[0] : NULL,
                         volmap, voltexmap, 
                         numvoxels, radscale, gridspacing, isovalue, gausslim);
  }

  free(xyzr);
  if (colors)
    free(colors);

  voltime = wkf_timer_timenow(timer);

  // draw the surface
  draw_trimesh(cmdList);

#if !defined(ARCH_BLUEWATERS)
  printf(" Done.\n");
#endif
  return 0;
}


// Extract the isosurface from the QuickSurf density map
int QuickSurf::get_trimesh(int &numverts, 
                           float *&v3fv, float *&n3fv, float *&c3fv, 
                           int &numfacets, int *&fiv) {
#if !defined(ARCH_BLUEWATERS)
  printf("Running marching cubes on CPU...\n");
#endif

  VolumetricData *surfvol; ///< Container used to generate isosurface on the CPU
  surfvol = new VolumetricData("molecular surface",
                               origin, xaxis, yaxis, zaxis,
                               numvoxels[0], numvoxels[1], numvoxels[2],
                               volmap);

  // XXX we should calculate the volume gradient only for those
  //     vertices we extract, since for this rep any changes to settings
  //     will require recomputation of the entire volume
  surfvol->compute_volume_gradient(); // calc gradients: smooth vertex normals
  gradtime = wkf_timer_timenow(timer);

  // trimesh polygonalized surface, max of 6 triangles per voxel
  const int stepsize = 1;
  s.clear();                              // initialize isosurface data
  s.compute(surfvol, isovalue, stepsize); // compute the isosurface

  mctime = wkf_timer_timenow(timer);

  s.vertexfusion(surfvol, 9, 9);          // identify and eliminate duplicated vertices
  s.normalize();                          // normalize interpolated gradient/surface normals

  if (s.numtriangles > 0) {
    if (voltexmap != NULL) {
      // assign per-vertex colors by a 3-D texture map
      s.set_color_voltex_rgb3fv(voltexmap);
    } else {
      // use a single color for the entire mesh
      s.set_color_rgb3fv(solidcolor);
    }
  }

  numverts = s.v.num() / 3;
  v3fv=&s.v[0];
  n3fv=&s.n[0];
  c3fv=&s.c[0];

  numfacets = s.numtriangles;
  fiv=&s.f[0];

  delete surfvol;

  mcverttime = wkf_timer_timenow(timer);
  reptime = mcverttime;

#if !defined(ARCH_BLUEWATERS)
  char strmsg[1024];
  sprintf(strmsg, "QuickSurf: %.3f [pre:%.3f vol:%.3f gr:%.3f mc:%.2f mcv:%.3f]",
          reptime, pretime, voltime-pretime, gradtime-voltime, 
          mctime-gradtime, mcverttime-mctime);

  msgInfo << strmsg << sendmsg;
#endif
 
  return 0;
}


int QuickSurf::draw_trimesh(VMDDisplayList *cmdList) {
  DispCmdTriMesh cmdTriMesh;

  int numverts=0;
  float *v=NULL, *n=NULL, *c=NULL;
  int numfacets=0;
  int *f=NULL;

  get_trimesh(numverts, v, n, c, numfacets, f);

  // Create a triangle mesh
  if (numfacets > 0) {
    cmdTriMesh.putdata(v, n, c, numverts, f, numfacets, 0, cmdList);
  }
 
  return 0;
}


QuickSurf::~QuickSurf() {
#if defined(VMDCUDA)
  if (cudaqs) {
    delete cudaqs;
  }
#endif

  if (voltexmap != NULL)
    free(voltexmap);
  voltexmap = NULL;

  wkf_timer_destroy(timer);
}


