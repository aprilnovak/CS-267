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

#define BS 32
#define BSS (BS*BS) // Block size squared (length of a linearized block)

#define REG_WIDTH 4 // Number of doubles per SIMD register

#define HAND_UNROLLED // Use hand-unrolled loops (else use real loops)

#define MAX_PREALLOCED (1024+32) // Biggest size to preallocate memory for

//#define HAVE_FMA 1 // Define if processor has SIMD FMA

#define STRIPES (BS/REG_WIDTH) // Because the SIMD instructions fit two doubles

#define ALIGNMENT 64


// To use for unrolled stripe loops, have BS/REG_WIDTH entries
// (It sure would be nice to have C++ templates...)
#if STRIPES == 4
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3);
#elif STRIPES == 6
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5);
#elif STRIPES == 7
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6);
#elif STRIPES == 8
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7);
#elif STRIPES == 9
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8);
#elif STRIPES == 10
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9);
#elif STRIPES == 12
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11);
#elif STRIPES == 13
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12);
#elif STRIPES == 14
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13);
#elif STRIPES == 16
  #define FOREACH_STRIPE(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15);
#else
  #error Unsupported number of stripes
#endif

#if BS == 24
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23);
#elif BS == 28
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23); \
    _op(24); \
    _op(25); \
    _op(26); \
    _op(27);
#elif BS == 32
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23); \
    _op(24); \
    _op(25); \
    _op(26); \
    _op(27); \
    _op(28); \
    _op(29); \
    _op(30); \
    _op(31);
#elif BS == 36
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23); \
    _op(24); \
    _op(25); \
    _op(26); \
    _op(27); \
    _op(28); \
    _op(29); \
    _op(30); \
    _op(31); \
    _op(32); \
    _op(33); \
    _op(34); \
    _op(35);
#elif BS == 40
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23); \
    _op(24); \
    _op(25); \
    _op(26); \
    _op(27); \
    _op(28); \
    _op(29); \
    _op(30); \
    _op(31); \
    _op(32); \
    _op(33); \
    _op(34); \
    _op(35); \
    _op(36); \
    _op(37); \
    _op(38); \
    _op(39);
#elif BS == 48
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23); \
    _op(24); \
    _op(25); \
    _op(26); \
    _op(27); \
    _op(28); \
    _op(29); \
    _op(30); \
    _op(31); \
    _op(32); \
    _op(33); \
    _op(34); \
    _op(35); \
    _op(36); \
    _op(37); \
    _op(38); \
    _op(39); \
    _op(40); \
    _op(41); \
    _op(42); \
    _op(43); \
    _op(44); \
    _op(45); \
    _op(46); \
    _op(47);
#elif BS == 52
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23); \
    _op(24); \
    _op(25); \
    _op(26); \
    _op(27); \
    _op(28); \
    _op(29); \
    _op(30); \
    _op(31); \
    _op(32); \
    _op(33); \
    _op(34); \
    _op(35); \
    _op(36); \
    _op(37); \
    _op(38); \
    _op(39); \
    _op(40); \
    _op(41); \
    _op(42); \
    _op(43); \
    _op(44); \
    _op(45); \
    _op(46); \
    _op(47); \
    _op(48); \
    _op(49); \
    _op(50); \
    _op(51);
#elif BS == 56
  #define FOREACH_IN_BLOCK(_op) \
    _op(0); \
    _op(1); \
    _op(2); \
    _op(3); \
    _op(4); \
    _op(5); \
    _op(6); \
    _op(7); \
    _op(8); \
    _op(9); \
    _op(10); \
    _op(11); \
    _op(12); \
    _op(13); \
    _op(14); \
    _op(15); \
    _op(16); \
    _op(17); \
    _op(18); \
    _op(19); \
    _op(20); \
    _op(21); \
    _op(22); \
    _op(23); \
    _op(24); \
    _op(25); \
    _op(26); \
    _op(27); \
    _op(28); \
    _op(29); \
    _op(30); \
    _op(31); \
    _op(32); \
    _op(33); \
    _op(34); \
    _op(35); \
    _op(36); \
    _op(37); \
    _op(38); \
    _op(39); \
    _op(40); \
    _op(41); \
    _op(42); \
    _op(43); \
    _op(44); \
    _op(45); \
    _op(46); \
    _op(47); \
    _op(48); \
    _op(49); \
    _op(50); \
    _op(51); \
    _op(52); \
    _op(53); \
    _op(54); \
    _op(55);
#else
  #error Unsupported block size
#endif


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

static double tmp[BS*BS] __attribute__ ((aligned(ALIGNMENT)));


#undef HAND_UNROLLED
//TODO: Refactor this with the non-partial code (as a macro)
static void partial (int _n, int i, int j, int iBS, int jBS)
{
  int nb = (_n+BS-1) / BS; // Size in terms of blocks
  int nn = nb * BS;
  int iSTRIPES = (iBS+REG_WIDTH-1) / REG_WIDTH;
  REG a, c[STRIPES], b;

      memset(tmp, 0, sizeof(double)*BS*BS);
      for (int k = 0; k < nb; k += 1)
      {

        for (int jj = 0; jj < jBS; ++jj)
        {
          #define JJ_LOOP_SETUP(_s) SET(c[_s],0);
          #ifndef HAND_UNROLLED
            for (int s = 0; s < iSTRIPES; ++s)
              JJ_LOOP_SETUP(s);
          #else
            FOREACH_STRIPE(JJ_LOOP_SETUP);
          #endif

          //TODO: It might help to manually unroll this
          for (int kk = 0; kk < BS; ++kk)
          {
            SET(b, BB(k,kk,j,jj));

            #define KK_LOOP_BODY(_s) { \
              LOAD(a, &AA(k,kk,i,_s*REG_WIDTH)); \
              FMADD(c[_s], a, b); }
            #ifndef HAND_UNROLLED
              for (int s = 0; s < iSTRIPES; ++s)
                KK_LOOP_BODY(s);
            #else
              FOREACH_STRIPE(KK_LOOP_BODY);
            #endif
          }

          #define JJ_INCREMENT_TEMPORARY(_s) { \
            LOAD(a, tmp+_s*REG_WIDTH+jj*BS); \
            ADD(c[_s], a, c[_s]); \
            STORE(c[_s], tmp+_s*REG_WIDTH+jj*BS); }
          #ifndef HAND_UNROLLED
            for (int s = 0; s < iSTRIPES; ++s)
              JJ_INCREMENT_TEMPORARY(s);
          #else
            FOREACH_STRIPE(JJ_INCREMENT_TEMPORARY);
          #endif
        } // jj

      } // k

      #define COPY_TEMP_TO_C_ROW(_jj) \
        memcpy( &OC(i,0,j,_jj), tmp+(_jj)*BS, sizeof(double)*iBS );
      #ifndef HAND_UNROLLED
        for (int jj = 0; jj < jBS; ++jj)
          COPY_TEMP_TO_C_ROW(jj);
      #else
        FOREACH_IN_BLOCK(COPY_TEMP_TO_C_ROW);
      #endif
}
#define HAND_UNROLLED 1



void square_dgemm (int _n, double* _A, double* _B, double* restrict _C)
{
  int nb = (_n+BS-1) / BS; // Size in terms of blocks
  int nn = nb * BS;

  int nbb = _n / BS; // Size in terms of blocks
  int remain = _n % BS;

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
      //memset(C+nn*i+_n, 0, (nn-_n)*sizeof(double));
    }
    //for (int i = _n; i < nn; ++i)
    //  memset(C+nn*i, 0, nn*sizeof(double));
  }
  else
  {
    memcpy(C,_C,nn*nn*sizeof(double));
  }

  prep(_n, _A, _B, _C);

  REG a, c[STRIPES], b;

  for (int i = 0; i < nbb; i += 1)
  {
    for (int j = 0; j < nbb; j += 1)
    {
      memset(tmp, 0, sizeof(double)*BS*BS);
      for (int k = 0; k < nb; k += 1)
      {

        for (int jj = 0; jj < BS; ++jj)
        {
          #ifndef HAND_UNROLLED
          for (int s = 0; s < STRIPES; ++s)
            SET(c[s],0);
          #else
            #define JJ_LOOP_SETUP(_s) SET(c[_s],0);
            FOREACH_STRIPE(JJ_LOOP_SETUP);
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
              #define KK_LOOP_BODY(_s) \
                LOAD(a, &AA(k,kk,i,_s*REG_WIDTH)); \
                FMADD(c[_s], a, b);
              FOREACH_STRIPE(KK_LOOP_BODY);
            #endif
          }

          #ifndef HAND_UNROLLED
          for (int s = 0; s < STRIPES; ++s)
          {
            LOAD(a, tmp+s*REG_WIDTH+jj*BS);
            ADD(c[s], a, c[s]);
            STORE(c[s], tmp+s*REG_WIDTH+jj*BS);

          }
          #else
            #define JJ_INCREMENT_TEMPORARY(_s) \
              LOAD(a, tmp+_s*REG_WIDTH+jj*BS); \
              ADD(c[_s], a, c[_s]); \
              STORE(c[_s], tmp+_s*REG_WIDTH+jj*BS);
            FOREACH_STRIPE(JJ_INCREMENT_TEMPORARY);
          #endif
        } // jj
      } // k

      #define COPY_TEMP_TO_C_ROW(_jj) \
        memcpy( &OC(i,0,j,_jj), tmp+(_jj)*BS, sizeof(double)*BS );
      #ifndef HAND_UNROLLED
        for (int jj = 0; jj < BS; ++jj)
          COPY_TEMP_TO_C_ROW(jj);
      #else
        FOREACH_IN_BLOCK(COPY_TEMP_TO_C_ROW);
      #endif

    } // j
    int j = nbb;
    if (remain) partial(_n, i, j, BS, remain);
  } // i
  if (remain)
  {
    int i = nbb;
    for (int j = 0; j < nbb; ++j)
      partial(_n, i, j, remain, BS);
    int j = nbb;
    partial(_n, i, j, remain, remain);
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
