#ifndef __GNUC__
#error Please compile with GCC!
#endif
#ifdef __INTEL_COMPILER
#error Please compile with GCC!
#endif
const char* dgemm_desc = "Explicit SIMD blocked dgemm.";

#include <x86intrin.h>
#include <string.h>
#include <stdbool.h>

#define BS 16
#define BSS (BS*BS) // Block size squared (length of a linearized block)

#define REG_WIDTH 4 // Number of doubles per SIMD register

#define HAND_UNROLLED // Use hand-unrolled loops (else use real loops)

#define MAX_PREALLOCED (1024+32) // Biggest size to preallocate memory for

//#define HAVE_FMA 1 // Define if processor has SIMD FMA

#define STRIPES (BS/REG_WIDTH) // Because the SIMD instructions fit two doubles

#define ALIGNMENT 64


// Access A/B/C in original format
#define OA(_i,_ii,_j,_jj) _A[ ((_i*BS)+(_ii)) * _n + ((_j*BS)+(_jj)) ] // Transposed
#define OB(_i,_ii,_j,_jj) _B[ ((_i*BS)+(_ii)) + _n * ((_j*BS)+(_jj)) ]
#define OC(_i,_ii,_j,_jj)  C[ ((_i*BS)+(_ii)) + nn * ((_j*BS)+(_jj)) ]


// Access using block-linear addressing
#define AA(_i,_ii,_j,_jj) (A[ ((_j)+((_i)*nb))*BSS + (_jj) + (_ii)*BS ])
#define BB(_i,_ii,_j,_jj) (B[ ((_j)+((_i)*nb))*BSS + (_jj) + (_ii)*BS ])


// Preallocated work buffers
static double SA[MAX_PREALLOCED*MAX_PREALLOCED] __attribute__ ((aligned(ALIGNMENT)));
static double SB[MAX_PREALLOCED*MAX_PREALLOCED] __attribute__ ((aligned(ALIGNMENT)));
static double SC[MAX_PREALLOCED*MAX_PREALLOCED] __attribute__ ((aligned(ALIGNMENT)));

// Dynamically allocated work buffers
static double * restrict DA __attribute__ ((aligned(ALIGNMENT)));
static double * restrict DB __attribute__ ((aligned(ALIGNMENT)));
static double * restrict DC __attribute__ ((aligned(ALIGNMENT)));
static int dynamic_size = 0; // "n" for dynamic buffers

// Pointers to either preallocated or dynamic work buffers
static double * restrict A __attribute__ ((aligned(ALIGNMENT)));
static double * restrict B __attribute__ ((aligned(ALIGNMENT)));
static double * restrict C __attribute__ ((aligned(ALIGNMENT)));


// This copies a single block if you pass BS,BS.
// If you pass less than BS, it'll copy less than a full block.
#define INTRA_BLOCK_COPY(_max_ii, _max_jj)      \
        for (int ii = 0; ii < _max_ii; ++ii)    \
        {                                       \
          for (int jj = 0; jj < _max_jj; ++jj)  \
          {                                     \
            AA(i,ii,j,jj) = OA(i,ii,j,jj);      \
            BB(i,ii,j,jj) = OB(i,ii,j,jj);      \
          }                                     \
        }                                       \


void prep (int _n, double* _A, double* _B, double* restrict _C)
{
  int nbb = _n / BS; // Size in terms of blocks
  int nb = (_n+BS-1) / BS; // Size in terms of blocks
  //int nn = nb * BS;
  int remain = _n % BS;

  if (!remain) // Evenly divisible case.  Easy.
  {
    for (int i = 0; i < nbb; ++i)
      for (int j = 0; j < nbb; ++j)
        INTRA_BLOCK_COPY(BS, BS);
  }
  else // Need to special-case the edges. :(
  {
    for (int i = 0; i < nbb; ++i)
    {
      for (int j = 0; j < nbb; ++j)
        INTRA_BLOCK_COPY(BS, BS);
      int j = nbb;
      INTRA_BLOCK_COPY(BS, remain);
    }
    int i = nbb;
    for (int j = 0; j < nbb; ++j)
      INTRA_BLOCK_COPY(remain, BS);
    int j = nbb;
    INTRA_BLOCK_COPY(remain,remain);
  }
}


static bool alloc_buffers (int n)
{
  // Allocate buffers of doubles of size n*n
  int s = n * n * sizeof(double);
  if (DA) _mm_free(A);
  if (DB) _mm_free(B);
  if (DC) _mm_free(C);
  DA = _mm_malloc(s, ALIGNMENT);
  DB = _mm_malloc(s, ALIGNMENT);
  DC = _mm_malloc(s, ALIGNMENT);
  if (!DA || !DB || !DC)
  {
    if (DA) _mm_free(A);
    if (DB) _mm_free(B);
    if (DC) _mm_free(C);
    dynamic_size = 0;
    return false;
  }
  dynamic_size = n;
  return true;
}


// In general, you probably want to use the greatest width that your
// platform supports!
#if REG_WIDTH == 2
  #define REG __m128d
  #define LOAD(_d,_s) _d = _mm_load_pd(_s)
  #define SET(_d,_v) _d = _mm_set1_pd(_v)
  #ifdef HAVE_FMA
    #define FMADD(_d, _s1, _s2) _d = _mm_fmadd_pd(_s1, _s2, _d)
  #else
    #define FMADD(_d, _s1, _s2) { _s1 = _mm_mul_pd(_s1, _s2); _d = _mm_add_pd(_d, _s1);}
  #endif
  #define STORE(_v, _d) _mm_store_pd(_d, _v)
  #define ADD(_d, _s1, _s2) _d = _mm_add_pd(_s1, _s2)
#elif REG_WIDTH == 4
  #define REG __m256d
  #define LOAD(_d,_s) _d = _mm256_load_pd(_s)
  #define SET(_d,_v) _d = _mm256_set1_pd(_v)
  #ifdef HAVE_FMA
    #define FMADD(_d, _s1, _s2) _d = _mm256_fmadd_pd(_s1, _s2, _d)
  #else
    #define FMADD(_d, _s1, _s2) { _s1 = _mm256_mul_pd(_s1, _s2); _d = _mm256_add_pd(_d, _s1);}
  #endif
  #define STORE(_v, _d) _mm256_store_pd(_d, _v)
  #define ADD(_d, _s1, _s2) _d = _mm256_add_pd(_s1, _s2)
#elif REG_WIDTH == 8
  #define REG __m512d
  #define LOAD(_d,_s) _d = _mm512_load_pd(_s)
  #define SET(_d,_v) _d = _mm512_set1_pd(_v)
  #ifdef HAVE_FMA
    #define FMADD(_d, _s1, _s2) _d = _mm512_fmadd_pd(_s1, _s2, _d)
  #else
    #define FMADD(_d, _s1, _s2) { _s1 = _mm512_mul_pd(_s1, _s2); _d = _mm512_add_pd(_d, _s1);}
  #endif
  #define STORE(_v, _d) _mm512_store_pd(_d, _v)
  #define ADD(_d, _s1, _s2) _d = _mm512_add_pd(_s1, _s2)
#else
  #error Unsupported register size.
#endif


void square_dgemm (int _n, double* _A, double* _B, double* restrict _C)
{
  int nb = (_n+BS-1) / BS; // Size in terms of blocks
  int nn = nb * BS;

  if (nn <= MAX_PREALLOCED)
  {
    A = SA;
    B = SB;
    C = SC;
  }
  else
  {
    // Use dynamic buffers
    if (nn > dynamic_size)
    {
      // We need to (re)allocate them.
      if (!alloc_buffers(nn + nn / 2))
      {
        // We tried to allocate them bigger than we needed, but it
        // failed.  Try exactly as big as we need and cross fingers.
        if (!alloc_buffers(nn))
        {
          return;
        }
      }
    }
    A = DA;
    B = DB;
    C = DC;
  }


  if (nn != _n)
  {
    for (int i = 0; i < _n; ++i)
    {
      memcpy(C+nn*i, _C+_n*i, _n*sizeof(double));
      memset(C+nn*i+_n, 0, (nn-_n)*sizeof(double));
    }
    for (int i = _n; i < nn; ++i)
      memset(C+nn*i, 0, nn*sizeof(double));
  }
  else
  {
    memcpy(C,_C,nn*nn*sizeof(double));
  }

  prep(_n, _A, _B, _C);

  REG a, c[STRIPES], b;

  for (int k = 0; k < nb; k += 1)
    for (int j = 0; j < nb; j += 1)
    {
      for (int i = 0; i < nb; i += 1)
      {

        for (int jj = 0; jj < BS; ++jj)
        {
          #ifndef HAND_UNROLLED
          for (int s = 0; s < STRIPES; ++s)
            SET(c[s],0);
          #else
            SET(c[0],0);
            SET(c[1],0);
            SET(c[2],0);
            SET(c[3],0);
          #endif

          //TODO: It might help to manually unroll this
          for (int kk = 0; kk < BS; ++kk)
          {
            SET(b, BB(k,kk,j,jj));

            #ifndef HAND_UNROLLED
              for (int s = 0; s < STRIPES; ++s)
              {
                LOAD(a, &AA(k,kk,i,s*REG_WIDTH));
                FMADD(c[s], a, b);
              }
            #else
              #if STRIPES != 4
                #error Hard coded for four stripes!
              #endif

              LOAD(a, &AA(k,kk,i,0*REG_WIDTH));
              FMADD(c[0], a, b);
              LOAD(a, &AA(k,kk,i,1*REG_WIDTH));
              FMADD(c[1], a, b);
              LOAD(a, &AA(k,kk,i,2*REG_WIDTH));
              FMADD(c[2], a, b);
              LOAD(a, &AA(k,kk,i,3*REG_WIDTH));
              FMADD(c[3], a, b);
            #endif
          }

          #ifndef HAND_UNROLLED
          for (int s = 0; s < STRIPES; ++s)
          {
            LOAD(a, &OC(i,s*REG_WIDTH,j,jj));
            ADD(c[s], a, c[s]);
            STORE(c[s], &OC(i,s*REG_WIDTH,j,jj));

          }
          #else
            LOAD(a, &OC(i,0*REG_WIDTH,j,jj));
            ADD(c[0], a, c[0]);
            STORE(c[0], &OC(i,0*REG_WIDTH,j,jj));
            LOAD(a, &OC(i,1*REG_WIDTH,j,jj));
            ADD(c[1], a, c[1]);
            STORE(c[1], &OC(i,1*REG_WIDTH,j,jj));
            LOAD(a, &OC(i,2*REG_WIDTH,j,jj));
            ADD(c[2], a, c[2]);
            STORE(c[2], &OC(i,2*REG_WIDTH,j,jj));
            LOAD(a, &OC(i,3*REG_WIDTH,j,jj));
            ADD(c[3], a, c[3]);
            STORE(c[3], &OC(i,3*REG_WIDTH,j,jj));
          #endif
        }
      }
    }

  if (nn != _n)
  {
    for (int i = 0; i < _n; i++)
    {
      memcpy(_C+_n*i, C+nn*i, _n*sizeof(double));
    }
  }
  else
  {
    memcpy(_C,C,nn*nn*sizeof(double));
  }
}
