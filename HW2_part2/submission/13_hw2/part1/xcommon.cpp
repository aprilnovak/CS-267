#include <stdio.h>
#include "xcommon.h"
#include "common.cpp"

//
//  Initialize the particle positions and velocities
//
void Xinit_particles( int n, Xparticle_t *p )
{
#ifdef RAND_SEED
    srand48( RAND_SEED );
#else
    srand48( time( NULL ) );
#endif

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;

    for( int i = 0; i < n; i++ )
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];

        //
        //  distribute particles evenly to ensure proper spacing
        //
#ifdef TWO_COORD
        p[i].x[0] = p[i].x[1] = size*(1.+(k%sx))/(1+sx);
        p[i].y[0] = p[i].y[1] = size*(1.+(k/sx))/(1+sy);
#else
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);
#endif

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void Xapply_force( int S, Xparticle_t &particle, Xparticle_t &neighbor, double & ax, double & ay, double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x NEXT - particle.x NEXT;
    double dy = neighbor.y NEXT - particle.y NEXT;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    if (r2 != 0)
    {
        if (r2/(cutoff*cutoff) < *dmin * (*dmin))
          *dmin = sqrt(r2)/cutoff;
        (*davg) += sqrt(r2)/cutoff;
        (*navg) ++;
    }

    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );



    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    ax += coef * dx;
    ay += coef * dy;
}

//
//  integrate the ODE
//
void Xmove( int S, Xparticle_t &p, double ax, double ay )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += ax * dt;
    p.vy += ay * dt;
    p.x NEXT  = p.x THIS + p.vx * dt;
    p.y NEXT  = p.y THIS + p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x NEXT < 0 || p.x NEXT > size )
    {
        p.x NEXT  = p.x NEXT < 0 ? -p.x NEXT : 2*size-p.x NEXT;
        p.vx = -p.vx;
    }
    while( p.y NEXT < 0 || p.y NEXT > size )
    {
        p.y NEXT  = p.y NEXT < 0 ? -p.y NEXT : 2*size-p.y NEXT;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void Xsave( int S, FILE *f, int n, Xparticle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x NEXT, p[i].y NEXT );
}

