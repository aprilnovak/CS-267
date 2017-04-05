#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <cassert>
#include <stdint.h>

#include "xcommon.h"
#include <omp.h>

#define xcutoff (1.5*0.01)

static double cutoff;

#if 1
  #define P(_p) *_p
  #define PP(_p) &_p
  static std::vector< std::vector<Xparticle_t> > bins[2];
#else
  #define P(_p) _p
  #define PP(_p) _p
  static std::vector< std::vector<Xparticle_t*> > bins[2];
#endif

static std::vector<omp_lock_t> locks;
static std::vector< std::vector< Xparticle_t * > > pending;

extern double size;
static int nbins;

#define CLAMP(_v, _l, _h) ( (_v) < (_l) ? (_l) : (((_v) > (_h)) ? (_h) : (_v)) )
#define XY(_x, _y) ((_x) + (_y) * nbins)

#define get_thread_for_bin(_b) CLAMP((_b) / bins_per_thread, 0, numthreads-1)

static inline int get_bin_for_particle (int S, Xparticle_t * p)
{
  int x = (int)(p->x NEXT / cutoff);
  int y = (int)(p->y NEXT / cutoff);

  x = CLAMP(x, 0, nbins - 1);
  y = CLAMP(y, 0, nbins - 1);

  return XY(x,y);
}


int main( int argc, char **argv )
{
    int numthreads;
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    Xparticle_t *particles = (Xparticle_t*) malloc( n * sizeof(Xparticle_t) );
    set_size( n );
    Xinit_particles( n, particles );

    nbins = (int)(size / xcutoff);
    cutoff = size / nbins;

    bins[0].resize( nbins*nbins );
    bins[1].resize( nbins*nbins );

    //  simulate a number of time steps
    double simulation_time = read_timer( );
    #pragma omp parallel
    {
      numthreads = omp_get_num_threads();
    }

    locks.resize( numthreads );
    for (int i = 0; i < locks.size(); i++)
      omp_init_lock(&locks[i]);

    pending.resize( numthreads );

    int bins_per_thread = (nbins*nbins) / numthreads;

    for (int i = 0; i < n; i++)
    {
      int new_bin = get_bin_for_particle(1, &particles[i]);
      int target_thread = get_thread_for_bin(new_bin);
      pending[target_thread].push_back(&particles[i]);
    }

    #pragma omp parallel
    {
      int this_thread = omp_get_thread_num();


      for (int p = 0; p < pending[this_thread].size(); p++)
      {
        Xparticle_t * pp = pending[this_thread][p];
        int new_bin = get_bin_for_particle(1, pp);
        bins[1][new_bin].push_back(P(pp));
      }
      pending[this_thread].resize(0);
      std::vector< Xparticle_t * >().swap( pending[this_thread] );

      #pragma omp barrier

      for( int step = 0; step < NSTEPS; step++ )
      {
          int S = step & 1;
          navg = 0;
          davg = 0.0;
          dmin = 1.0;

          #pragma omp for schedule(static, bins_per_thread)
          for (int xy = 0; xy < nbins*nbins; xy++)
          {
            bins[S][xy].resize(0);
          }
          #pragma omp for schedule(static, bins_per_thread) reduction (+:navg) reduction(+:davg)
          for (int xy = 0; xy < nbins*nbins; xy++)
          {
            int x = xy % nbins;
            int y = xy / nbins;

            assert (XY(x,y) == xy);
            int yl = CLAMP(y - 1, 0, nbins - 1);
            int yh = CLAMP(y + 1, 0, nbins - 1);
            assert (XY(x,y) < bins[!S].size());
            int xl = CLAMP(x - 1, 0, nbins - 1);
            int xh = CLAMP(x + 1, 0, nbins - 1);

            for (int p = 0; p < bins[!S][XY(x,y)].size(); p++)
            {
              Xparticle_t * pp = PP(bins[!S][XY(x,y)][p]);

              assert ( get_bin_for_particle(S, pp) == xy );
              assert ( get_thread_for_bin(get_bin_for_particle(S, pp)) == this_thread );

              double ax = 0, ay = 0;
              for (int y2 = yl; y2 <= yh; y2++)
              {
                for (int x2 = xl; x2 <= xh; x2++)
                {
                  assert (XY(x2,y2) < bins[!S].size());

                  for (int p2 = 0; p2 < bins[!S][XY(x2,y2)].size(); p2++)
                  {
                    Xparticle_t * pp2 = PP(bins[!S][XY(x2,y2)][p2]);
                    Xapply_force(S, *pp, *pp2, ax, ay, &dmin, &davg, &navg);
                  }
                }
              }
              Xmove(!S, *pp, ax, ay);
              int new_bin = get_bin_for_particle(!S, pp);
              if (new_bin == xy || get_thread_for_bin(new_bin) == this_thread)
              {
                // We own this bin!
                bins[S][new_bin].push_back(P(pp));
              }
              else
              {
                int target_thread = get_thread_for_bin(new_bin);
                omp_set_lock(&locks[target_thread]);
                pending[target_thread].push_back(pp);
                omp_unset_lock(&locks[target_thread]);
              }
            }
          }

          #pragma omp barrier


          for (int i = 0; i < pending[this_thread].size(); i++)
          {
            Xparticle_t * pp = pending[this_thread][i];
            int new_bin = get_bin_for_particle(!S, pp);
            bins[S][new_bin].push_back(P(pp));
          }
          pending[this_thread].resize(0);

          #pragma omp barrier

#if STATS_ENABLE
          if( find_option( argc, argv, "-no" ) == -1 )
          {
            // Computing statistical data
            #pragma omp master
            if (navg) {
              absavg +=  davg/navg;
              nabsavg++;
            }
            #pragma omp critical
            if (dmin < absmin) absmin = dmin;

            //  save if necessary
            #pragma omp master
            if( fsave && (step%SAVEFREQ) == 0 )
            {
              int i = 0;
              for (int xy = 0; xy < nbins*nbins; xy++)
              {
                int c = bins[S][xy].size();
                memcpy(particles + i, PP(bins[S][xy][0]), c * sizeof(Xparticle_t));
                i += c;
              }
              Xsave( S, fsave, n, particles );

#ifdef MOVIE_MODE
              // Easier for making movies...
              char buf[100];
              sprintf(buf, "coords.%03i", step/SAVEFREQ);
              FILE * fp = fopen(buf, "w");
              for( i = 0; i < n; i++ )
                fprintf( fp, "%g %g\n", particles[i].x NEXT, particles[i].y NEXT );
              fclose(fp);
#endif
            }
          }
#endif
      }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, threads = %d, simulation time = %g seconds", n, numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
      //
      //  -The minimum distance absmin between 2 particles during the run of
      //   the simulation
      //  -A Correct simulation will have particles stay at greater than 0.4
      //   (of cutoff) with typical values between .7-.8
      //  -A simulation where particles don't interact correctly will be less
      //   than 0.4 (of cutoff) with typical values between .01-.05
      //
      //  -The average distance absavg is ~.95 when most particles are
      //   interacting correctly and ~.66 when no particles are interacting
      //
      printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
      if (absmin < 0.4) printf ("\nsome particle is not interacting");
      if (absavg < 0.8) printf ("\nmost particles are not interacting");
    }
    printf("\n");

    // Printing summary data
    if( fsum) fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    // Clearing space
    if( fsum ) fclose( fsum );
    free( particles );
    if( fsave ) fclose( fsave );

    return 0;
}
