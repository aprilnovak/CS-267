#ifndef __GNUC__
#error Please compile with the GNU compiler
#endif
#ifdef __INTEL_COMPILER
#error Please compile with the GNU compiler
#endif
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <time.h>
#include <string.h>
#include <math.h>
#include "common.h"

#include <unistd.h>


#define ROOT 0
MPI_Datatype PARTICLE;


// This is exactly the init_particles() code, but we've commented out
// the seeding of the random generator.
extern double size;
static void init_particles_no_seed( int n, particle_t *p )
{
    //srand48( time( NULL ) );

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
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}


// This is much like a vector, but it can be resized without initializing
// the uncovered elements (which turns out to be expensive).
template <class T>
class Buffer
{
private:
    T * _buf;
    int _size;
    int _capacity;

public:
    Buffer (int size = 0)
    : _buf(new T[size]), _size(size), _capacity(size)
    {
    }

    ~Buffer (void)
    {
        delete [] _buf;
    }

    void reserve (int reserve)
    {
        if (reserve <= _capacity) return;
        T * nb = new T[reserve];
        if (_buf)
        {
            memcpy(nb, _buf, sizeof(_buf[0]) * _size);
            delete [] _buf;
        }
        _buf = nb;
        _capacity = reserve;
    }

    int capacity (void) const
    {
        return _capacity;
    }

    int size (void) const
    {
        return _size;
    }

    void resize (int size)
    {
        reserve(size);
        _size = size;
    }

    T & operator [] (const int e) const
    {
        return _buf[e];
    }

    void push_back (T & item)
    {
        if (_size >= _capacity) reserve(max(_size + 1, _capacity * 2));
        _buf[_size] = item;
        _size++;
    }
};

typedef Buffer<particle_t> ParticleBuffer;


MPI_Comm world;

static void init_random (int rank)
{
    long data[3] = {0}; // Seed plus some test values
    if (rank == 0)
    {
        long real_seed;
        #ifdef RAND_SEED
            real_seed = RAND_SEED;
        #else
            srand48(time(NULL));
            do
            {
                real_seed = lrand48();
            } while (real_seed == 0);
        #endif
        srand48(real_seed);
        data[0] = real_seed;
        for (int i = 1; i < sizeof(data)/sizeof(data[0]); i++)
        {
            data[i] = lrand48();
        }
    }

    MPI_Bcast(data, sizeof(data)/sizeof(data[0]), MPI_LONG, ROOT, world);

    if (data[0] == 0)
    {
        fprintf(stderr, "Random synchronization failure! (Bad seed)\n");
        exit(1);
    }
    srand48(data[0]);

    for (int i = 1; i < sizeof(data)/sizeof(data[0]); i++)
    {
        long v = lrand48();
        if (data[i] != v)
        {
            fprintf(stderr, "Random synchronization failure! (%li != %li)\n", data[i], v);
            exit(1);
        }
    }

    srand48(data[0]);
}


#define xcutoff (1.25*0.01)

static double cutoff;

#define P(_p) *(_p)
#define PP(_p) &(_p)
typedef std::vector<particle_t> cell_t;
static std::vector< cell_t > cells[2];

extern double size;
static int ncells;

#define CLAMP(_v, _l, _h) ( (_v) < (_l) ? (_l) : (((_v) > (_h)) ? (_h) : (_v)) )
#define XY(_x, _y) ((_x) + (_y) * ncells)
#define GET_Y(_xy) ((_xy) / ncells)
#define GET_X(_xy) ((_xy) % ncells)


static inline int get_cell_for_particle (particle_t * p)
{
  int x = (int)(p->x / cutoff);
  int y = (int)(p->y / cutoff);

  x = CLAMP(x, 0, ncells - 1);
  y = CLAMP(y, 0, ncells - 1);

  return XY(x,y);
}



static void put_in_cells (ParticleBuffer & particles, std::vector< cell_t > & cells)
{
    for (int p = 0; p < particles.size(); p++)
    {
        particle_t * pp = PP(particles[p]);
        int new_cell = get_cell_for_particle(pp);
        cells[new_cell].push_back(P(pp));
    }
}

static int count_particles (std::vector<cell_t> & cells)
{
    int c = 0;
    for (int i = 0; i < cells.size(); i++)
        c += cells[i].size();
    return c;
}

static int count_particles (std::vector<cell_t> & cells, int lo, int hi)
{
    int c = 0;
    for (int i = 0; i < cells.size(); i++)
    {
        if (GET_Y(i) < lo || GET_Y(i) > hi) continue;
        c += cells[i].size();
    }
    return c;
}


//  benchmarking program
int main( int argc, char **argv )
{
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;

    //  process command line parameters
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

    //  set up MPI
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    set_size( n );
    ncells = (int)(size / xcutoff);
    cutoff = size / ncells;
    assert (cutoff >= 0.01);
    cells[0].resize(ncells*ncells);
    cells[1].resize(ncells*ncells);

    MPI_Group everyone;
    MPI_Comm_group(MPI_COMM_WORLD, &everyone);
    if (n_proc > ncells)
    {
        // More processes than we need.
        int range[] = {0, ncells-1, 1};
        MPI_Group_range_incl(everyone, 1, &range, &everyone);
        n_proc =  ncells;
    }
    MPI_Comm_create(MPI_COMM_WORLD, everyone, &world);
    if (world == MPI_COMM_NULL)
    {
        // We weren't included
        MPI_Finalize();
        return 0;
    }

    //  allocate generic resources
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


    //MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    init_random(rank);

    int stripe_height = ncells / n_proc;
    int leftover_height = ncells % n_proc;

    if (rank < leftover_height) stripe_height += 1;

    int ymin = rank * stripe_height;
    int ymax = CLAMP(ymin + stripe_height - 1, 0, ncells-1);

    if (rank >= leftover_height)
    {
        ymin += leftover_height;
        ymax += leftover_height;
    }

    if (rank == n_proc - 1)
    {
        assert (ymax == ncells - 1);
    }

    bool has_uhalo = rank > 0; // Have upward halo?
    bool has_dhalo = rank < (n_proc - 1); // Have downward halo?
    int yminh = ymin - (has_uhalo ? 1 : 0); // min including halo
    int ymaxh = ymax + (has_dhalo ? 1 : 0); // max including halo

    int first_cell = ymin * ncells; // First cell we own
    int last_cell = ymax * ncells + ncells - 1; // Last cell we own


    #ifdef DEBUGGING
    for (int i = 0; i < n_proc; i++)
    {
        MPI_Barrier(world);
        if (i == rank) printf("rank:%i  leftover_height:%i  stripe_height:%i  ncells:%i  n_proc:%i  ymin:%i  ymax:%i  has_uhalo:%i  has_dhalo:%i  first_cell:%i  last_cell:%i\n",
                               rank, leftover_height, stripe_height, ncells, n_proc, ymin, ymax, has_uhalo, has_dhalo, first_cell, last_cell);
    }
    #endif


    // Particles to be sent
    // We need these to be double-buffered because we start filling it
    // in one step before we know that it's sent from the previous step.
    ParticleBuffer tx_uhalos[2];
    ParticleBuffer tx_dhalos[2];

    // Buffers for received halo particles
    ParticleBuffer rx_uhalo;
    ParticleBuffer rx_dhalo;

    // We don't know the number of particles that we'll receive ahead
    // of time.  So we'll allocate enough memory for lots of particles.
    // This burns memory, but makes things very easy and means we don't
    // need to send a "length" message.
    rx_uhalo.reserve(n);
    rx_dhalo.reserve(n);


    // Initialize the particles
    ParticleBuffer initial_particles(n);
    init_particles_no_seed( n, &initial_particles[0] );
    put_in_cells(initial_particles, cells[0]);

    //printf("%i: %i - %i (%i - %i)\n", rank, ymin, ymax, yminh, ymaxh);

    MPI_Request rx_uhalo_req;
    MPI_Request rx_dhalo_req;
    MPI_Request tx_uhalo_req;
    MPI_Request tx_dhalo_req;


    //  simulate a number of time steps
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        #ifdef DEBUGGING
        int count = 0;
        #endif
        int toggle_cur = step & 1;
        int toggle_nxt = (toggle_cur ^ 1);
        std::vector<cell_t> & cur_cells = cells[toggle_cur];
        std::vector<cell_t> & nxt_cells = cells[toggle_nxt];

        ParticleBuffer & tx_uhalo = tx_uhalos[toggle_cur];
        ParticleBuffer & tx_dhalo = tx_dhalos[toggle_cur];

        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        tx_uhalo.resize(0);
        tx_dhalo.resize(0);

        #ifdef DEBUGGING
        int going_up = 0;
        int going_down = 0;
        int up_halo = 0;
        int down_halo = 0;
        #endif

        if (has_uhalo)
        {
            if (step > 0)
            {
                MPI_Status status;
                MPI_Wait(&rx_uhalo_req, &status);
                int rx_size;
                MPI_Get_count(&status, PARTICLE, &rx_size);
                rx_uhalo.resize(rx_size);
                put_in_cells(rx_uhalo, cur_cells);
            }
            //rx_uhalo.resize(rx_uhalo.capacity());
            int r = MPI_Irecv(&rx_uhalo[0], rx_uhalo.capacity(), PARTICLE, rank - 1, 0, world, &rx_uhalo_req);
            assert (r == MPI_SUCCESS);
        }

        int need_dhalo_at = CLAMP(ymax - 1, ymin, ymax);
        bool uhalo_sent = false;

        //  compute all forces
        for (int y = ymin; y <= ymax; y++)
        {
            int yl = CLAMP(y - 1, yminh, ymaxh);
            int yh = CLAMP(y + 1, yminh, ymaxh);

            if (has_dhalo && y == need_dhalo_at)
            {
                if (step > 0)
                {
                    // Okay, now we really need our downward halo!
                    // We need it at ymax-1 because we may now own new
                    // particles that came from below, and they'll be in
                    // ymax, meaning they might affect ones in ymax-1.
                    MPI_Status status;
                    MPI_Wait(&rx_dhalo_req, &status);
                    int rx_size;
                    MPI_Get_count(&status, PARTICLE, &rx_size);
                    rx_dhalo.resize(rx_size);
                    put_in_cells(rx_dhalo, cur_cells);
                }
                //rx_dhalo.resize(rx_dhalo.capacity());
                int r = MPI_Irecv(&rx_dhalo[0], rx_dhalo.capacity(), PARTICLE, rank + 1, 0, world, &rx_dhalo_req);
                assert (r == MPI_SUCCESS);
            }

            for (int x = 0; x < ncells; x++)
            {
                int xl = CLAMP(x - 1, 0, ncells - 1);
                int xh = CLAMP(x + 1, 0, ncells - 1);

                for (int p = 0; p < cur_cells[XY(x,y)].size(); p++)
                {
                    #ifdef DEBUGGING
                        count++;
                    #endif
                    particle_t pp = cur_cells[XY(x,y)][p];
                    pp.ax = pp.ay = 0;
                    for (int y2 = yl; y2 <= yh; y2++)
                    {
                        for (int x2 = xl; x2 <= xh; x2++)
                        {
                            for (int p2 = 0; p2 < cur_cells[XY(x2,y2)].size(); p2++)
                            {
                                particle_t * pp2 = PP(cur_cells[XY(x2,y2)][p2]);
                                apply_force( pp, *pp2, &dmin, &davg, &navg );
                            }
                        }
                    }

                    move(pp);

                    int new_cell = get_cell_for_particle(&pp);


                    if (GET_Y(new_cell) <= ymin)
                        tx_uhalo.push_back(P(&pp));
                    if (GET_Y(new_cell) >= ymax)
                        tx_dhalo.push_back(P(&pp));

                    #ifdef DEBUGGING
                    if (GET_Y(new_cell) < ymin)
                        going_up++;
                    if (GET_Y(new_cell) > ymax)
                        going_down++;
                    if (GET_Y(new_cell) == ymin)
                        up_halo++;
                    if (GET_Y(new_cell) == ymax)
                        down_halo++;
                    #endif
                    if ( (GET_Y(new_cell) >= yminh) && (GET_Y(new_cell) <= ymaxh) )
                        nxt_cells[new_cell].push_back(P(&pp));
                }
            }

            if (y == ymin + 2 && has_uhalo)
            {
                uhalo_sent = true;
                MPI_Status status;
                //printf("%i tx uhalo %p\n", rank, &tx_uhalo_req);
                if (step > 0) MPI_Wait(&tx_uhalo_req, &status);
                int r = MPI_Isend(&tx_uhalo[0], tx_uhalo.size(), PARTICLE, rank - 1, 0, world, &tx_uhalo_req);
                assert (r == MPI_SUCCESS);
            }
        }


        //  save current step if necessary (slightly different semantics than in other codes)
        if( find_option( argc, argv, "-no" ) == -1 )
        {
            if( fsave && (step%SAVEFREQ) == 0 )
            {
                MPI_Barrier(world);

                // Collect all our local particles...
                std::vector<particle_t> particles;
                for (int y = ymin; y <= ymax; y++)
                {
                    for (int x = 0; x < ncells; x++)
                    {
                        for (int p = 0; p < cur_cells[XY(x,y)].size(); p++)
                        {
                            particle_t * pp = PP(cur_cells[XY(x,y)][p]);
                            particles.push_back(P(pp));
                        }
                    }
                }

                if (rank != 0)
                {
                    // If we are not rank 0, send our local particles to rank 0.

                    int r = MPI_Send(&particles[0], particles.size(), PARTICLE, 0, 1, world);
                    assert (r == MPI_SUCCESS);
                }
                else
                {
                    // If we are rank 0, collect everyone else's particles too, then save.
                    // For simplicity, receive in order (not necessarily efficient).
                    for (int i = 1; i < n_proc; i++)
                    {
                        MPI_Status status;
                        int cur_size = particles.size();
                        particles.resize(n);
                        assert ((n - cur_size) >= 0);
                        int r = MPI_Recv(&particles[cur_size], n - cur_size, PARTICLE, i, 1, world, &status);
                        assert (r == MPI_SUCCESS);
                        int additional_size;
                        MPI_Get_count(&status, PARTICLE, &additional_size);
                        assert ( (cur_size + additional_size) <= n );
                        particles.resize(cur_size + additional_size);
                    }
                    assert (particles.size() == n);
                    save( fsave, n, &particles[0] );

#ifdef MOVIE_MODE
                    // Easier for making movies...
                    char buf[100];
                    sprintf(buf, "coords.%03i", step/SAVEFREQ);
                    FILE * fp = fopen(buf, "w");
                    for( int i = 0; i < n; i++ )
                      fprintf( fp, "%g %g\n", particles[i].x, particles[i].y );
                    fclose(fp);
#endif
                }
            }
        }


        if( find_option( argc, argv, "-no" ) == -1 )
        {
            MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,world);
            MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,world);
            MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,world);

            #ifdef DEBUGGING
            int total_count = 0;
            MPI_Reduce(&count,&total_count,1,MPI_INT,MPI_SUM,0,world);
            if (rank == 0) assert(total_count == n);
            #endif

            if (rank == 0)
            {
                // Computing statistical data
                if (rnavg)
                {
                    absavg +=  rdavg/rnavg;
                    nabsavg++;
                }
                if (rdmin < absmin) absmin = rdmin;
            }
        }





        if (has_uhalo && !uhalo_sent)
        {
            MPI_Status status;
            //printf("%i tx uhalo %p\n", rank, &tx_uhalo_req);
            if (step > 0) MPI_Wait(&tx_uhalo_req, &status);
            int r = MPI_Isend(&tx_uhalo[0], tx_uhalo.size(), PARTICLE, rank - 1, 0, world, &tx_uhalo_req);
            assert (r == MPI_SUCCESS);
        }

        if (has_dhalo)
        {
            MPI_Status status;
            //printf("%i tx uhalo %p\n", rank, &tx_dhalo_req);
            if (step > 0) MPI_Wait(&tx_dhalo_req, &status);
            int r = MPI_Isend(&tx_dhalo[0], tx_dhalo.size(), PARTICLE, rank + 1, 0, world, &tx_dhalo_req);
            assert (r == MPI_SUCCESS);
        }

        // Clear cells
        for (int y = yminh; y <= ymaxh; y++)
        {
            for (int x = 0; x < ncells; x++)
            {
                cur_cells[XY(x,y)].resize(0);
            }
        }

        #ifdef DEBUGGING
        MPI_Barrier(world);

        for (int i = 0; i < n_proc; i++)
        {
            if (i == rank)
            {
                printf("%-4i %-4i going_up:%-4i up_halo:%-4i going_down:%-4i down_halo:%-4i old_count:%-4i count:%-4i\n", rank, step, going_up, up_halo, going_down, down_halo, count, count_particles(nxt_cells, ymin, ymax));
                fflush(stdout);
            }
            MPI_Barrier(world);
        }
        if (rank == n_proc-1) printf("\n");
        fflush(stdout);
        #endif

    }
    simulation_time = read_timer( ) - simulation_time;

    #ifdef DEBUGGING
    fflush(stdout);
    MPI_Barrier(world);
    #endif

    if (rank == 0)
    {
        printf( "n = %d, simulation time = %g seconds", n, simulation_time);

        if( find_option( argc, argv, "-no" ) == -1 )
        {
            if (nabsavg) absavg /= nabsavg;
            //
            //  -The minimum distance absmin between 2 particles during the run of the simulation
            //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
            //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
            //
            //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
            //
            printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
            if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
            if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
        }
        printf("\n");

        // Printing summary data
        if( fsum) fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }

    //  release resources
    if ( fsum ) fclose( fsum );
    if( fsave ) fclose( fsave );

    MPI_Finalize( );

    return 0;
}
