#ifndef __CS267_XCOMMON_H__
#define __CS267_XCOMMON_H__

#include "common.h"

#ifdef TWO_COORD
  #define THIS [!S]
  #define NEXT [S]
#else
  #define THIS
  #define NEXT
#endif

typedef struct
{
#ifdef TWO_COORD
  double x[2];
  double y[2];
#else
  double x;
  double y;
#endif
  double vx;
  double vy;
} Xparticle_t;

void Xinit_particles( int n, Xparticle_t *p );
void Xapply_force( int S, Xparticle_t &particle, Xparticle_t &neighbor , double & ax, double & ay, double *dmin, double *davg, int *navg);
void Xmove( int S, Xparticle_t &p, double ax, double ay );


void Xsave( int S, FILE *f, int n, Xparticle_t *p );

#endif
