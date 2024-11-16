#ifndef KLU_H
#define KLU_H

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <getopt.h>
#include <stdint.h>
#include <stdbool.h>
#include <errno.h>
#include <inttypes.h>
#include <stddef.h>
#include <limits.h>
#include <stdarg.h>
#include <ctype.h>
#include "sched.h"

// #include "btf.h"
// #include "btf_internal.h"

// #include "klu.h" 重复定义
#endif

#define KLU_symbolic klu_symbolic
#define KLU_numeric klu_numeric
#define KLU_common klu_common
#define KLU_analyze klu_analyze
#define KLU_analyze_given klu_analyze_given
#define KLU_free_symbolic klu_free_symbolic
#define KLU_free_numeric klu_free_numeric
#define KLU_defaults klu_defaults
#define KLU_factor klu_factor
#define KLU_refactor klu_refactor
#define KLU_solve klu_solve
#define Entry double

#define EPS 1e-9
#define DRAND() ((double)rand() / ((double)RAND_MAX / 2.0F)) - 1.0F;

#define KLU_OK 0
#define EMPTY (-1)
#define KLU_OUT_OF_MEMORY (-2)
#define KLU_INVALID (-3) 
#define KLU_TOO_LARGE (-4) 
#define COLAMD_STATS 20
#define Int int32_t
#define MAX(a,b) (((a) > (b)) ?  (a) : (b))
#define MIN(a,b) (((a) < (b)) ?  (a) : (b))
#define UNASSIGNED (-1) 
#define UNVISITED (-2) 
#define BTF_FLIP(j) (-(j)-2)
#define BTF_ISFLIPPED(j) ((j) < -1)
#define BTF_UNFLIP(j) ((BTF_ISFLIPPED (j)) ? BTF_FLIP (j) : (j))
#define KLU_SINGULAR (1)  
#define REAL(c)                     (c)
#define SCALAR_IS_NAN(x)        ((x) != (x))
#define SCALAR_IS_ZERO(x)       ((x) == 0.)
#define SCALAR_IS_NONZERO(x)    ((x) != 0.)
#define SCALAR_IS_LTZERO(x)     ((x) < 0.)
#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO (x)) ? -(x) : (x))
#define ABS(s,a)                    { (s) = SCALAR_ABS (a) ; }
#define CLEAR(c)                    { (c) = 0. ; }
#define SCALE_DIV_ASSIGN(a,c,s)     { a = c / s ; }
#define SCALAR_IS_ZERO(x)       ((x) == 0.)
#define IS_ZERO(a)                  SCALAR_IS_ZERO (a)
#define BYTES(type,n) (sizeof (type) * (n))
#define CEILING(b,u)  (((b)+(u)-1) / (u))
#define UNITS(type,n) (CEILING (BYTES (type,n), sizeof (Unit)))
#define DUNITS(type,n) (ceil (BYTES (type, (double) n) / sizeof (Unit)))  //ceil向上取整
#define INT_OVERFLOW(x) ((!((x) * (1.0+1e-8) <= (double) Int_MAX)) \
                        || SCALAR_IS_NAN (x))
#define MULT_SUB(c,a,b)             { (c) -= (a) * (b) ; }
#define GET_POINTER(LU, Xip, Xlen, Xi, Xx, k, xlen) \
{ \
    Unit *xp = LU + Xip [k] ; \
    xlen = Xlen [k] ; \
    Xi = (Int *) xp ; \
    Xx = (Entry *) (xp + UNITS (Int, xlen)) ; \
}
#define DIV(c,a,b)                  { (c) = (a) / (b) ; }
#define SCALE_DIV(c,s)              { (c) /= (s) ; }

#define SCHED_FIFO      1

typedef signed char int8_t;
typedef unsigned char   uint8_t;
typedef short  int16_t;
typedef unsigned short  uint16_t;
typedef int  int32_t;
typedef unsigned   uint32_t;
typedef double Unit ;


typedef struct klu_common_struct
{

    /* ---------------------------------------------------------------------- */
    /* parameters */
    /* ---------------------------------------------------------------------- */

    double tol ;            /* pivot tolerance for diagonal preference */
    double memgrow ;        /* realloc memory growth size for LU factors */
    double initmem_amd ;    /* init. memory size with AMD: c*nnz(L) + n */
    double initmem ;        /* init. memory size: c*nnz(A) + n */
    double maxwork ;        /* maxwork for BTF, <= 0 if no limit */

    int btf ;               /* use BTF pre-ordering, or not */
    int ordering ;          /* 0: AMD, 1: COLAMD, 2: user P and Q,
                             * 3: user function */
    int scale ;             /* row scaling: -1: none (and no error check),
                             * 0: none, 1: sum, 2: max */

    /* pointer to user ordering function */
    int32_t (*user_order) (int32_t, int32_t *, int32_t *, int32_t *,
        struct klu_common_struct *) ;

    /* pointer to user data, passed unchanged as the last parameter to the
     * user ordering function (optional, the user function need not use this
     * information). */
    void *user_data ;

    int halt_if_singular ;      /* how to handle a singular matrix:
        * FALSE: keep going.  Return a Numeric object with a zero U(k,k).  A
        *   divide-by-zero may occur when computing L(:,k).  The Numeric object
        *   can be passed to klu_solve (a divide-by-zero will occur).  It can
        *   also be safely passed to klu_refactor.
        * TRUE: stop quickly.  klu_factor will free the partially-constructed
        *   Numeric object.  klu_refactor will not free it, but will leave the
        *   numerical values only partially defined.  This is the default. */

    /* ---------------------------------------------------------------------- */
    /* statistics */
    /* ---------------------------------------------------------------------- */

    int status ;                /* KLU_OK if OK, < 0 if error */
    int nrealloc ;              /* # of reallocations of L and U */

    int32_t structural_rank ;       /* 0 to n-1 if the matrix is structurally rank
        * deficient (as determined by maxtrans).  -1 if not computed.  n if the
        * matrix has full structural rank.  This is computed by klu_analyze
        * if a BTF preordering is requested. */

    int32_t numerical_rank ;        /* First k for which a zero U(k,k) was found,
        * if the matrix was singular (in the range 0 to n-1).  n if the matrix
        * has full rank. This is not a true rank-estimation.  It just reports
        * where the first zero pivot was found.  -1 if not computed.
        * Computed by klu_factor and klu_refactor. */

    int32_t singular_col ;          /* n if the matrix is not singular.  If in the
        * range 0 to n-1, this is the column index of the original matrix A that
        * corresponds to the column of U that contains a zero diagonal entry.
        * -1 if not computed.  Computed by klu_factor and klu_refactor. */

    int32_t noffdiag ;      /* # of off-diagonal pivots, -1 if not computed */

    double flops ;      /* actual factorization flop count, from klu_flops */
    double rcond ;      /* crude reciprocal condition est., from klu_rcond */
    double condest ;    /* accurate condition est., from klu_condest */
    double rgrowth ;    /* reciprocal pivot rgrowth, from klu_rgrowth */
    double work ;       /* actual work done in BTF, in klu_analyze */

    size_t memusage ;   /* current memory usage, in bytes */
    size_t mempeak ;    /* peak memory usage, in bytes */

} klu_common ;


typedef struct
{
    /* A (P,Q) is in upper block triangular form.  The kth block goes from
     * row/col index R [k] to R [k+1]-1.  The estimated number of nonzeros
     * in the L factor of the kth block is Lnz [k]. 
     */

    /* only computed if the AMD ordering is chosen: */
    double symmetry ;   /* symmetry of largest block */
    double est_flops ;  /* est. factorization flop count */
    double lnz, unz ;   /* estimated nz in L and U, including diagonals */
    double *Lnz ;       /* size n, but only Lnz [0..nblocks-1] is used */

    /* computed for all orderings: */
    int32_t
        n,              /* input matrix A is n-by-n */
        nz,             /* # entries in input matrix */
        *P,             /* size n */
        *Q,             /* size n */
        *R,             /* size n+1, but only R [0..nblocks] is used */
        nzoff,          /* nz in off-diagonal blocks */
        nblocks,        /* number of blocks */
        maxblock,       /* size of largest block */
        ordering,       /* ordering used (0:AMD, 1:COLAMD, 2:given, ... */
        do_btf ;        /* whether or not BTF preordering was requested */

    /* only computed if BTF preordering requested */
    int32_t structural_rank ;   /* 0 to n-1 if the matrix is structurally rank
                        * deficient.  -1 if not computed.  n if the matrix has
                        * full structural rank */

} klu_symbolic ;

typedef struct
{
    /* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

    int32_t n ;             /* A is n-by-n */
    int32_t nblocks ;       /* number of diagonal blocks */
    int32_t lnz ;           /* actual nz in L, including diagonal */
    int32_t unz ;           /* actual nz in U, including diagonal */
    int32_t max_lnz_block ; /* max actual nz in L in any one block, incl. diag */
    int32_t max_unz_block ; /* max actual nz in U in any one block, incl. diag */
    int32_t *Pnum ;         /* size n. final pivot permutation */
    int32_t *Pinv ;         /* size n. inverse of final pivot permutation */

    /* LU factors of each block */
    int32_t *Lip ;          /* size n. pointers into LUbx[block] for L */
    int32_t *Uip ;          /* size n. pointers into LUbx[block] for U */
    int32_t *Llen ;         /* size n. Llen [k] = # of entries in kth column of L */
    int32_t *Ulen ;         /* size n. Ulen [k] = # of entries in kth column of U */
    void **LUbx ;       /* L and U indices and entries (excl. diagonal of U) */
    size_t *LUsize ;    /* size of each LUbx [block], in sizeof (Unit) */
    void *Udiag ;       /* diagonal of U */

    /* scale factors; can be NULL if no scaling */
    double *Rs ;        /* size n. Rs [i] is scale factor for row i */

    /* permanent workspace for factorization and solve */
    size_t worksize ;   /* size (in bytes) of Work */ 
    void *Work ;        /* workspace */
    void *Xwork ;       /* alias into Numeric->Work */
    int32_t *Iwork ;        /* alias into Numeric->Work */

    /* off-diagonal entries in a conventional compressed-column sparse matrix */
    int32_t *Offp ;         /* size n+1, column pointers */
    int32_t *Offi ;         /* size nzoff, row indices */
    void *Offx ;        /* size nzoff, numerical values */
    int32_t nzoff ;

}  klu_numeric ;

KLU_symbolic *KLU_alloc_symbolic (Int n, Int *Ap, Int *Ai, KLU_common *Common) ;

void *KLU_malloc        /* returns pointer to the newly malloc'd block */
(
    /* ---- input ---- */
    size_t n,           /* number of items */
    size_t size,        /* size of each item */
    /* --------------- */
    KLU_common *Common
);

void *KLU_free          /* always returns NULL */
(
    /* ---- in/out --- */
    void *p,            /* block of memory to free */
    /* ---- input --- */
    size_t n,           /* size of block to free, in # of items */
    size_t size,        /* size of each item */
    /* --------------- */
    KLU_common *Common
);

void *SuiteSparse_malloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to malloc */
    size_t size_of_item     /* sizeof each item */
);

void *SuiteSparse_free      /* always returns NULL */
(
    void *p                 /* block to free */
);

int klu_defaults
(
    klu_common *Common
) ;

klu_symbolic *klu_analyze
(
    /* inputs, not modified */
    int32_t n,              /* A is n-by-n */
    int32_t Ap [ ],         /* size n+1, column pointers */
    int32_t Ai [ ],         /* size nz, row indices */
    klu_common *Common
) ;


#ifndef SUITESPARSE_CONFIG_H
#define SUITESPARSE_CONFIG_H

//------------------------------------------------------------------------------
// SuiteSparse-wide ANSI C11 #include files
//------------------------------------------------------------------------------

#undef  SuiteSparse_long
#undef  SuiteSparse_long_max
#undef  SuiteSparse_long_idd
#undef  SuiteSparse_long_id

#define SuiteSparse_long int64_t
#define SuiteSparse_long_max INT64_MAX
#define SuiteSparse_long_idd PRId64
#define SuiteSparse_long_id "%" SuiteSparse_long_idd

//------------------------------------------------------------------------------
// OpenMP
//------------------------------------------------------------------------------

#if defined ( _OPENMP )

    #include <omp.h>
    #define SUITESPARSE_OPENMP_MAX_THREADS       omp_get_max_threads ( )
    #define SUITESPARSE_OPENMP_GET_NUM_THREADS   omp_get_num_threads ( )
    #define SUITESPARSE_OPENMP_GET_WTIME         omp_get_wtime ( )
    #define SUITESPARSE_OPENMP_GET_THREAD_ID     omp_get_thread_num ( )

#else

    // OpenMP not available
    #define SUITESPARSE_OPENMP_MAX_THREADS       (1)
    #define SUITESPARSE_OPENMP_GET_NUM_THREADS   (1)
    #define SUITESPARSE_OPENMP_GET_WTIME         (0)
    #define SUITESPARSE_OPENMP_GET_THREAD_ID     (0)

#endif

//------------------------------------------------------------------------------
// MATLAB/Octave
//------------------------------------------------------------------------------

#if defined ( MATLAB_MEX_FILE )
#include "mex.h"
#include "matrix.h"
#endif

//------------------------------------------------------------------------------
// string and token handling macros
//------------------------------------------------------------------------------

// SUITESPARSE_STR: convert the content of x into a string "x"
#define SUITESPARSE_XSTR(x) SUITESPARSE_STR(x)
#define SUITESPARSE_STR(x) #x

// SUITESPARSE_CAT(x,y): concatenate two tokens
#define SUITESPARSE_CAT2(x,y) x ## y
#define SUITESPARSE_CAT(x,y) SUITESPARSE_CAT2(x,y)

//------------------------------------------------------------------------------
// determine which compiler is in use
//------------------------------------------------------------------------------

#define SUITESPARSE_COMPILER_NVCC    0
#define SUITESPARSE_COMPILER_ICX     0
#define SUITESPARSE_COMPILER_ICC     0
#define SUITESPARSE_COMPILER_CLANG   0
#define SUITESPARSE_COMPILER_GCC     0
#define SUITESPARSE_COMPILER_MSC     0
#define SUITESPARSE_COMPILER_XLC     0

#if defined ( __NVCC__ )

    // NVIDIA nvcc compiler
    #undef  SUITESPARSE_COMPILER_NVCC
    #define SUITESPARSE_COMPILER_NVCC    1

    #define SUITESPARSE_COMPILER_MAJOR __CUDACC_VER_MAJOR__
    #define SUITESPARSE_COMPILER_MINOR __CUDACC_VER_MINOR__
    #define SUITESPARSE_COMPILER_SUB   __CUDACC_VER_BUILD__
    #define SUITESPARSE_COMPILER_NAME  "nvcc"

#elif defined ( __INTEL_CLANG_COMPILER )

    // Intel icx compiler, 2022.0.0 based on clang/llvm 14.0.0
    #undef  SUITESPARSE_COMPILER_ICX
    #define SUITESPARSE_COMPILER_ICX     1

    #define SUITESPARSE_COMPILER_MAJOR __INTEL_CLANG_COMPILER
    #define SUITESPARSE_COMPILER_MINOR 0
    #define SUITESPARSE_COMPILER_SUB   0
    #define SUITESPARSE_COMPILER_NAME  __VERSION__

#elif defined ( __INTEL_COMPILER )

    // Intel icc compiler: 2021.5.0 uses "gcc 7.5 mode"
    #undef  SUITESPARSE_COMPILER_ICC
    #define SUITESPARSE_COMPILER_ICC     1

    #define SUITESPARSE_COMPILER_MAJOR __INTEL_COMPILER
    #define SUITESPARSE_COMPILER_MINOR __INTEL_COMPILER_UPDATE
    #define SUITESPARSE_COMPILER_SUB   0
    #define SUITESPARSE_COMPILER_NAME  __VERSION__

#elif defined ( __clang__ )

    // clang
    #undef  SUITESPARSE_COMPILER_CLANG
    #define SUITESPARSE_COMPILER_CLANG   1

    #define SUITESPARSE_COMPILER_MAJOR __clang_major__
    #define SUITESPARSE_COMPILER_MINOR __clang_minor__
    #define SUITESPARSE_COMPILER_SUB   __clang_patchlevel__
    #define SUITESPARSE_COMPILER_NAME  "clang " __clang_version__

#elif defined ( __xlC__ )

    // xlc
    #undef  SUITESPARSE_COMPILER_XLC
    #define SUITESPARSE_COMPILER_XLC     1

    #define SUITESPARSE_COMPILER_MAJOR ( __xlC__ / 256 )
    #define SUITESPARSE_COMPILER_MINOR \
        ( __xlC__ - 256 * SUITESPARSE_COMPILER_MAJOR)
    #define SUITESPARSE_COMPILER_SUB   0
    #define SUITESPARSE_COMPILER_NAME  "IBM xlc " SUITESPARSE_XSTR (__xlC__)

#elif defined ( __GNUC__ )

    // gcc
    #undef  SUITESPARSE_COMPILER_GCC
    #define SUITESPARSE_COMPILER_GCC     1

    #define SUITESPARSE_COMPILER_MAJOR __GNUC__
    #define SUITESPARSE_COMPILER_MINOR __GNUC_MINOR__
    #define SUITESPARSE_COMPILER_SUB   __GNUC_PATCHLEVEL__
    #define SUITESPARSE_COMPILER_NAME  "GNU gcc " \
        SUITESPARSE_XSTR (__GNUC__) "." \
        SUITESPARSE_XSTR (__GNUC_MINOR__) "." \
        SUITESPARSE_XSTR (__GNUC_PATCHLEVEL__)

#elif defined ( _MSC_VER )

    // Microsoft Visual Studio (cl compiler)
    #undef  SUITESPARSE_COMPILER_MSC
    #define SUITESPARSE_COMPILER_MSC     1

    #define SUITESPARSE_COMPILER_MAJOR ( _MSC_VER / 100 )
    #define SUITESPARSE_COMPILER_MINOR \
        ( _MSC_VER - 100 * SUITESPARSE_COMPILER_MAJOR)
    #define SUITESPARSE_COMPILER_SUB   0
    #define SUITESPARSE_COMPILER_NAME \
        "Microsoft Visual Studio " SUITESPARSE_XSTR (_MSC_VER)

#else

    // other compiler
    #define SUITESPARSE_COMPILER_MAJOR 0
    #define SUITESPARSE_COMPILER_MINOR 0
    #define SUITESPARSE_COMPILER_SUB   0
    #define SUITESPARSE_COMPILER_NAME  "other C compiler"

#endif

//------------------------------------------------------------------------------
// malloc.h: required include file for Microsoft Visual Studio
//------------------------------------------------------------------------------

#if SUITESPARSE_COMPILER_MSC
    #include <malloc.h>
#endif

// this was formerly "extern", or "__declspec ..." for Windows.
#define SUITESPARSE_PUBLIC

//------------------------------------------------------------------------------
// determine the ANSI C version
//------------------------------------------------------------------------------

#ifdef __STDC_VERSION__
// ANSI C17: 201710L
// ANSI C11: 201112L
// ANSI C99: 199901L
// ANSI C95: 199409L
#define SUITESPARSE_STDC_VERSION __STDC_VERSION__
#else
// assume ANSI C90 / C89
#define SUITESPARSE_STDC_VERSION 199001L
#endif

//------------------------------------------------------------------------------
// handle the restrict keyword
//------------------------------------------------------------------------------

#if defined ( __cplusplus )

    // C++ does not have the "restrict" keyword
    #define SUITESPARSE_RESTRICT

#elif SUITESPARSE_COMPILER_MSC

    // MS Visual Studio
    #define SUITESPARSE_RESTRICT __restrict

#elif SUITESPARSE_COMPILER_NVCC

    // NVIDIA nvcc
    #define SUITESPARSE_RESTRICT __restrict__

#elif SUITESPARSE_STDC_VERSION >= 199901L

    // ANSI C99 or later
    #define SUITESPARSE_RESTRICT restrict

#else

    // ANSI C95 and earlier: no restrict keyword
    #define SUITESPARSE_RESTRICT

#endif

#ifdef __cplusplus
extern "C"
{
#endif


void *(*SuiteSparse_config_malloc_func_get (void)) (size_t);
void *(*SuiteSparse_config_calloc_func_get (void)) (size_t, size_t);
void *(*SuiteSparse_config_realloc_func_get (void)) (void *, size_t);
void (*SuiteSparse_config_free_func_get (void)) (void *);
int (*SuiteSparse_config_printf_func_get (void)) (const char *, ...);
double (*SuiteSparse_config_hypot_func_get (void)) (double, double);
int (*SuiteSparse_config_divcomplex_func_get (void)) (double, double, double, double, double *, double *);

// The SuiteSparse_config_*_set methods modify the contents of the struct:
void SuiteSparse_config_malloc_func_set (void *(*malloc_func) (size_t));
void SuiteSparse_config_calloc_func_set (void *(*calloc_func) (size_t, size_t));
void SuiteSparse_config_realloc_func_set (void *(*realloc_func) (void *, size_t));
void SuiteSparse_config_free_func_set (void (*free_func) (void *));
void SuiteSparse_config_printf_func_set (int (*printf_func) (const char *, ...));
void SuiteSparse_config_hypot_func_set (double (*hypot_func) (double, double));
void SuiteSparse_config_divcomplex_func_set (int (*divcomplex_func) (double, double, double, double, double *, double *));

// The SuiteSparse_config_*_func methods are wrappers that call the function
// pointers in the struct.  Note that there is no wrapper for the printf_func.
// See the SUITESPARSE_PRINTF macro instead.
void *SuiteSparse_config_malloc (size_t s) ;
void *SuiteSparse_config_calloc (size_t n, size_t s) ;
void *SuiteSparse_config_realloc (void *, size_t s) ;
void SuiteSparse_config_free (void *) ;
double SuiteSparse_config_hypot (double x, double y) ;
int SuiteSparse_config_divcomplex
(
    double xr, double xi, double yr, double yi, double *zr, double *zi
) ;

void SuiteSparse_start ( void ) ;   // called to start SuiteSparse

void SuiteSparse_finish ( void ) ;  // called to finish SuiteSparse

void *SuiteSparse_calloc    // pointer to allocated block of memory
(
    size_t nitems,          // number of items to calloc (>=1 is enforced)
    size_t size_of_item     // sizeof each item
) ;

void *SuiteSparse_realloc   // pointer to reallocated block of memory, or
                            ///to original block if the realloc failed.
(
    size_t nitems_new,      // new number of items in the object
    size_t nitems_old,      // old number of items in the object
    size_t size_of_item,    // sizeof each item
    void *p,                // old object to reallocate
    int *ok                 // 1 if successful, 0 otherwise
) ;

void SuiteSparse_tic    // start the timer
(
    double tic [2]      // output, contents undefined on input
) ;

double SuiteSparse_toc  // return time in seconds since last tic
(
    double tic [2]      // input: from last call to SuiteSparse_tic
) ;

double SuiteSparse_time  // returns current wall clock time in seconds
(
    void
) ;

// returns sqrt (x^2 + y^2), computed reliably
double SuiteSparse_hypot (double x, double y) ;

// complex division of c = a/b
int SuiteSparse_divcomplex
(
    double ar, double ai,       // real and imaginary parts of a
    double br, double bi,       // real and imaginary parts of b
    double *cr, double *ci      // real and imaginary parts of c
) ;

// determine which timer to use, if any
#ifndef NTIMER
    // SuiteSparse_config itself can be compiled without OpenMP,
    // but other packages can themselves use OpenMP.  In this case,
    // those packages should use omp_get_wtime() directly.  This can
    // be done via the SUITESPARSE_TIME macro, defined below:
    #define SUITESPARSE_TIMER_ENABLED
    #define SUITESPARSE_HAVE_CLOCK_GETTIME
    #define SUITESPARSE_CONFIG_TIMER omp_get_wtime
    #if defined ( SUITESPARSE_TIMER_ENABLED )
        #if defined ( _OPENMP )
            // Avoid indirection through the library if the compilation unit
            // including this header happens to use OpenMP.
            #define SUITESPARSE_TIME (omp_get_wtime ( ))
        #else
            #define SUITESPARSE_TIME (SuiteSparse_time ( ))
        #endif
    #else
        // No timer is available
        #define SUITESPARSE_TIME (0)
    #endif
#else
    // The SuiteSparse_config timer is explictly disabled;
    // use the OpenMP timer omp_get_wtime if available.
    #undef SUITESPARSE_TIMER_ENABLED
    #undef SUITESPARSE_HAVE_CLOCK_GETTIME
    #undef SUITESPARSE_CONFIG_TIMER
    #if defined ( _OPENMP )
        #define SUITESPARSE_CONFIG_TIMER omp_get_wtime
        #define SUITESPARSE_TIME (omp_get_wtime ( ))
    #else
        #define SUITESPARSE_CONFIG_TIMER none
        #define SUITESPARSE_TIME (0)
    #endif
#endif

// SuiteSparse printf macro
#define SUITESPARSE_PRINTF(params)                          \
{                                                           \
    int (*printf_func) (const char *, ...) ;                \
    printf_func = SuiteSparse_config_printf_func_get ( ) ;  \
    if (printf_func != NULL)                                \
    {                                                       \
        (void) (printf_func) params ;                       \
    }                                                       \
}


int SuiteSparse_version     // returns SUITESPARSE_VERSION
(
    // output, not defined on input.  Not used if NULL.  Returns
    // the three version codes in version [0..2]:
    // version [0] is SUITESPARSE_MAIN_VERSION
    // version [1] is SUITESPARSE_SUB_VERSION
    // version [2] is SUITESPARSE_SUBSUB_VERSION
    int version [3]
) ;

#define SUITESPARSE_HAS_VERSION_FUNCTION

#define SUITESPARSE_DATE "Aug 20, 2024"
#define SUITESPARSE_MAIN_VERSION    7
#define SUITESPARSE_SUB_VERSION     8
#define SUITESPARSE_SUBSUB_VERSION  2

// version format x.y
#define SUITESPARSE_VER_CODE(main,sub) ((main) * 1000 + (sub))
#define SUITESPARSE_VERSION SUITESPARSE_VER_CODE(7, 8)

// version format x.y.z
#define SUITESPARSE__VERCODE(main,sub,patch) \
    (((main)*1000ULL + (sub))*1000ULL + (patch))
#define SUITESPARSE__VERSION SUITESPARSE__VERCODE(7,8,2)

#if defined ( BLAS_NO_UNDERSCORE )

    // no name mangling, use lower case
    #define SUITESPARSE_FORTRAN(name,NAME)  name
    #define SUITESPARSE__FORTRAN(name,NAME) name

#elif defined ( BLAS_UNDERSCORE )

    // append an underscore, use lower case
    #define SUITESPARSE_FORTRAN(name,NAME)  name ## _
    #define SUITESPARSE__FORTRAN(name,NAME) name ## _

#else

    // let CMake decide how C calls Fortran
    #define SUITESPARSE_FORTRAN(name,NAME) name##_
    #define SUITESPARSE__FORTRAN(name,NAME) name##_

#endif


#if defined ( BLAS64 )

    // override the BLAS found by CMake, and force a 64-bit interface
    #define SUITESPARSE_BLAS_INT int64_t

#elif defined ( BLAS32 )

    // override the BLAS found by CMake, and force a 32-bit interface
    #define SUITESPARSE_BLAS_INT int32_t

#else

    // let CMake determine the size of the integer in the Fortran BLAS
    #define SUITESPARSE_BLAS_INT int32_t

#endif

// SUITESPARSE_TO_BLAS_INT: convert an integer k to a BLAS integer K and set ok
// to false if the conversion changes its value.  This is implemented as a
// macro so that can work with any type of the integer k.
#define SUITESPARSE_TO_BLAS_INT(K,k,ok)         \
    SUITESPARSE_BLAS_INT K = (k) ;              \
    ok = ok && ((sizeof (K) >= sizeof (k)) || ((int64_t)(K) == (int64_t)(k))) ;

#if defined ( BLAS64__SUFFIX )

    // The suffix includes an undersore (such as "_64"), so the Fortran name
    // must be processed with the SUITESPARSE__FORTRAN macro.
    #define SUITESPARSE_G(name,NAME) SUITESPARSE__FORTRAN(name,NAME)
    #define SUITESPARSE_F(name,NAME)                            \
        SUITESPARSE_G (SUITESPARSE_CAT (name, BLAS64__SUFFIX),  \
                       SUITESPARSE_CAT (NAME, BLAS64__SUFFIX))
    #define SUITESPARSE_BLAS(name,NAME) SUITESPARSE_F(name,NAME)

#elif defined ( BLAS64_SUFFIX )

    // The suffix does not include an undersore, and neither do the original
    // names of the BLAS and LAPACK routines.  Thus, the Fortran name must be
    // processed with the SUITESPARSE_FORTRAN macro.
    #define SUITESPARSE_G(name,NAME) SUITESPARSE_FORTRAN(name,NAME)
    #define SUITESPARSE_F(name,NAME)                            \
        SUITESPARSE_G (SUITESPARSE_CAT (name, BLAS64_SUFFIX),  \
                       SUITESPARSE_CAT (NAME, BLAS64_SUFFIX))
    #define SUITESPARSE_BLAS(name,NAME) SUITESPARSE_F(name,NAME)

#else

    // No suffix is need, so the final Fortran name includes no suffix.
    #define SUITESPARSE_BLAS(name,NAME) SUITESPARSE_FORTRAN(name,NAME)

#endif

//------------------------------------------------------------------------------
// C names of Fortan BLAS and LAPACK functions used by SuiteSparse
//------------------------------------------------------------------------------

// double
#define SUITESPARSE_BLAS_DTRSV      SUITESPARSE_BLAS ( dtrsv  , DTRSV  )
#define SUITESPARSE_BLAS_DGEMV      SUITESPARSE_BLAS ( dgemv  , DGEMV  )
#define SUITESPARSE_BLAS_DTRSM      SUITESPARSE_BLAS ( dtrsm  , DTRSM  )
#define SUITESPARSE_BLAS_DGEMM      SUITESPARSE_BLAS ( dgemm  , DGEMM  )
#define SUITESPARSE_BLAS_DSYRK      SUITESPARSE_BLAS ( dsyrk  , DSYRK  )
#define SUITESPARSE_BLAS_DGER       SUITESPARSE_BLAS ( dger   , DGER   )
#define SUITESPARSE_BLAS_DSCAL      SUITESPARSE_BLAS ( dscal  , DSCAL  )
#define SUITESPARSE_BLAS_DNRM2      SUITESPARSE_BLAS ( dnrm2  , DNRM2  )

#define SUITESPARSE_LAPACK_DPOTRF   SUITESPARSE_BLAS ( dpotrf , DPOTRF )
#define SUITESPARSE_LAPACK_DLARF    SUITESPARSE_BLAS ( dlarf  , DLARF  )
#define SUITESPARSE_LAPACK_DLARFG   SUITESPARSE_BLAS ( dlarfg , DLARFG )
#define SUITESPARSE_LAPACK_DLARFT   SUITESPARSE_BLAS ( dlarft , DLARFT )
#define SUITESPARSE_LAPACK_DLARFB   SUITESPARSE_BLAS ( dlarfb , DLARFB )

// double complex
#define SUITESPARSE_BLAS_ZTRSV      SUITESPARSE_BLAS ( ztrsv  , ZTRSV  )
#define SUITESPARSE_BLAS_ZGEMV      SUITESPARSE_BLAS ( zgemv  , ZGEMV  )
#define SUITESPARSE_BLAS_ZTRSM      SUITESPARSE_BLAS ( ztrsm  , ZTRSM  )
#define SUITESPARSE_BLAS_ZGEMM      SUITESPARSE_BLAS ( zgemm  , ZGEMM  )
#define SUITESPARSE_BLAS_ZHERK      SUITESPARSE_BLAS ( zherk  , ZHERK  )
#define SUITESPARSE_BLAS_ZGERU      SUITESPARSE_BLAS ( zgeru  , ZGERU  )
#define SUITESPARSE_BLAS_ZSCAL      SUITESPARSE_BLAS ( zscal  , ZSCAL  )
#define SUITESPARSE_BLAS_DZNRM2     SUITESPARSE_BLAS ( dznrm2 , DZNRM2 )

#define SUITESPARSE_LAPACK_ZPOTRF   SUITESPARSE_BLAS ( zpotrf , ZPOTRF )
#define SUITESPARSE_LAPACK_ZLARF    SUITESPARSE_BLAS ( zlarf  , ZLARF  )
#define SUITESPARSE_LAPACK_ZLARFG   SUITESPARSE_BLAS ( zlarfg , ZLARFG )
#define SUITESPARSE_LAPACK_ZLARFT   SUITESPARSE_BLAS ( zlarft , ZLARFT )
#define SUITESPARSE_LAPACK_ZLARFB   SUITESPARSE_BLAS ( zlarfb , ZLARFB )

// single
#define SUITESPARSE_BLAS_STRSV      SUITESPARSE_BLAS ( strsv  , STRSV  )
#define SUITESPARSE_BLAS_SGEMV      SUITESPARSE_BLAS ( sgemv  , SGEMV  )
#define SUITESPARSE_BLAS_STRSM      SUITESPARSE_BLAS ( strsm  , STRSM  )
#define SUITESPARSE_BLAS_SGEMM      SUITESPARSE_BLAS ( sgemm  , SGEMM  )
#define SUITESPARSE_BLAS_SSYRK      SUITESPARSE_BLAS ( ssyrk  , SSYRK  )
#define SUITESPARSE_BLAS_SGER       SUITESPARSE_BLAS ( sger   , SGER   )
#define SUITESPARSE_BLAS_SSCAL      SUITESPARSE_BLAS ( sscal  , SSCAL  )
#define SUITESPARSE_BLAS_SNRM2      SUITESPARSE_BLAS ( snrm2  , SNRM2  )

#define SUITESPARSE_LAPACK_SPOTRF   SUITESPARSE_BLAS ( spotrf , SPOTRF )
#define SUITESPARSE_LAPACK_SLARF    SUITESPARSE_BLAS ( slarf  , SLARF  )
#define SUITESPARSE_LAPACK_SLARFG   SUITESPARSE_BLAS ( slarfg , SLARFG )
#define SUITESPARSE_LAPACK_SLARFT   SUITESPARSE_BLAS ( slarft , SLARFT )
#define SUITESPARSE_LAPACK_SLARFB   SUITESPARSE_BLAS ( slarfb , SLARFB )

// single complex
#define SUITESPARSE_BLAS_CTRSV      SUITESPARSE_BLAS ( ctrsv  , CTRSV  )
#define SUITESPARSE_BLAS_CGEMV      SUITESPARSE_BLAS ( cgemv  , CGEMV  )
#define SUITESPARSE_BLAS_CTRSM      SUITESPARSE_BLAS ( ctrsm  , CTRSM  )
#define SUITESPARSE_BLAS_CGEMM      SUITESPARSE_BLAS ( cgemm  , CGEMM  )
#define SUITESPARSE_BLAS_CHERK      SUITESPARSE_BLAS ( cherk  , CHERK  )
#define SUITESPARSE_BLAS_CGERU      SUITESPARSE_BLAS ( cgeru  , CGERU  )
#define SUITESPARSE_BLAS_CSCAL      SUITESPARSE_BLAS ( cscal  , CSCAL  )
#define SUITESPARSE_BLAS_SCNRM2     SUITESPARSE_BLAS ( scnrm2 , SCNRM2 )

#define SUITESPARSE_LAPACK_CPOTRF   SUITESPARSE_BLAS ( cpotrf , CPOTRF )
#define SUITESPARSE_LAPACK_CLARF    SUITESPARSE_BLAS ( clarf  , CLARF  )
#define SUITESPARSE_LAPACK_CLARFG   SUITESPARSE_BLAS ( clarfg , CLARFG )
#define SUITESPARSE_LAPACK_CLARFT   SUITESPARSE_BLAS ( clarft , CLARFT )
#define SUITESPARSE_LAPACK_CLARFB   SUITESPARSE_BLAS ( clarfb , CLARFB )

//------------------------------------------------------------------------------
// SuiteSparse_BLAS_library: return name of BLAS library found
//------------------------------------------------------------------------------

// Returns the name of the BLAS library found by SuiteSparse_config

const char *SuiteSparse_BLAS_library ( void ) ;

//------------------------------------------------------------------------------
// SuiteSparse_BLAS_integer_size: return sizeof (SUITESPARSE_BLAS_INT)
//------------------------------------------------------------------------------

size_t SuiteSparse_BLAS_integer_size ( void ) ;

#ifdef __cplusplus
}
#endif
#endif



#ifndef AMD_H
#define AMD_H

/* make it easy for C++ programs to include AMD */
#ifdef __cplusplus
extern "C" {
#endif

int amd_order  /* returns AMD_OK, AMD_OK_BUT_JUMBLED,
                                    * AMD_INVALID, or AMD_OUT_OF_MEMORY */
(
    int32_t n,                     /* A is n-by-n.  n must be >= 0. */
    const int32_t Ap [ ],          /* column pointers for A, of size n+1 */
    const int32_t Ai [ ],          /* row indices of A, of size nz = Ap [n] */
    int32_t P [ ],                 /* output permutation, of size n */
    double Control [ ],     /* input Control settings, of size AMD_CONTROL */
    double Info [ ]         /* output Info statistics, of size AMD_INFO */
) ;

int amd_l_order  /* see above for description */
(
    int64_t n,
    const int64_t Ap [ ],
    const int64_t Ai [ ],
    int64_t P [ ],
    double Control [ ],
    double Info [ ]
) ;

void amd_2
(
    int32_t n,
    int32_t Pe [ ],
    int32_t Iw [ ],
    int32_t Len [ ],
    int32_t iwlen,
    int32_t pfree,
    int32_t Nv [ ],
    int32_t Next [ ], 
    int32_t Last [ ],
    int32_t Head [ ],
    int32_t Elen [ ],
    int32_t Degree [ ],
    int32_t W [ ],
    double Control [ ],
    double Info [ ]
) ;

void amd_l2
(
    int64_t n,
    int64_t Pe [ ],
    int64_t Iw [ ],
    int64_t Len [ ],
    int64_t iwlen,
    int64_t pfree,
    int64_t Nv [ ],
    int64_t Next [ ], 
    int64_t Last [ ],
    int64_t Head [ ],
    int64_t Elen [ ],
    int64_t Degree [ ],
    int64_t W [ ],
    double Control [ ],
    double Info [ ]
) ;

/* ------------------------------------------------------------------------- */
/* amd_valid */
/* ------------------------------------------------------------------------- */


int amd_valid
(
    int32_t n_row,                 /* # of rows */
    int32_t n_col,                 /* # of columns */
    const int32_t Ap [ ],          /* column pointers, of size n_col+1 */
    const int32_t Ai [ ]           /* row indices, of size Ap [n_col] */
) ;

int amd_l_valid
(
    int64_t n_row,
    int64_t n_col,
    const int64_t Ap [ ],
    const int64_t Ai [ ]
) ;

/* ------------------------------------------------------------------------- */
/* AMD Control and Info arrays */
/* ------------------------------------------------------------------------- */

/* amd_defaults:  sets the default control settings */
void amd_defaults   (double Control [ ]) ;
void amd_l_defaults (double Control [ ]) ;

/* amd_control: prints the control settings */
void amd_control    (double Control [ ]) ;
void amd_l_control  (double Control [ ]) ;

/* amd_info: prints the statistics */
void amd_info       (double Info [ ]) ;
void amd_l_info     (double Info [ ]) ;

// amd_version: return AMD version.  The version array is returned with
// version [0..2] = {AMD_MAIN_VERSION, AMD_SUB_VERSION, AMD_SUBSUB_VERSION}
void amd_version (int version [3]) ;

#ifdef __cplusplus
}
#endif

#define AMD_CONTROL 5          /* size of Control array */
#define AMD_INFO 20            /* size of Info array */

/* contents of Control */
#define AMD_DENSE 0            /* "dense" if degree > Control [0] * sqrt (n) */
#define AMD_AGGRESSIVE 1    /* do aggressive absorption if Control [1] != 0 */

/* default Control settings */
#define AMD_DEFAULT_DENSE 10.0          /* default "dense" degree 10*sqrt(n) */
#define AMD_DEFAULT_AGGRESSIVE 1    /* do aggressive absorption by default */

/* contents of Info */
#define AMD_STATUS 0           /* return value of amd_order and amd_l_order */
#define AMD_N 1                /* A is n-by-n */
#define AMD_NZ 2      /* number of nonzeros in A */ 
#define AMD_SYMMETRY 3         /* symmetry of pattern (1 is sym., 0 is unsym.) */
#define AMD_NZDIAG 4           /* # of entries on diagonal */
#define AMD_NZ_A_PLUS_AT 5  /* nz in A+A' */
#define AMD_NDENSE 6           /* number of "dense" rows/columns in A */
#define AMD_MEMORY 7           /* amount of memory used by AMD */
#define AMD_NCMPA 8            /* number of garbage collections in AMD */
#define AMD_LNZ 9     /* approx. nz in L, excluding the diagonal */
#define AMD_NDIV 10            /* number of fl. point divides for LU and LDL' */
#define AMD_NMULTSUBS_LDL 11 /* number of fl. point (*,-) pairs for LDL' */
#define AMD_NMULTSUBS_LU 12  /* number of fl. point (*,-) pairs for LU */
#define AMD_DMAX 13             /* max nz. in any column of L, incl. diagonal */


#define AMD_OK 0           /* success */
#define AMD_OUT_OF_MEMORY -1        
#define AMD_INVALID -2             
#define AMD_OK_BUT_JUMBLED 1        



#define AMD_DATE "June 20, 2024"
#define AMD_MAIN_VERSION   3
#define AMD_SUB_VERSION    3
#define AMD_SUBSUB_VERSION 3

#define AMD_VERSION_CODE(main,sub) SUITESPARSE_VER_CODE(main,sub)
#define AMD_VERSION AMD_VERSION_CODE(3,3)

#define AMD__VERSION SUITESPARSE__VERCODE(3,3,3)
#if !defined (SUITESPARSE__VERSION) || \
    (SUITESPARSE__VERSION < SUITESPARSE__VERCODE(7,8,0))
#error "AMD 3.3.3 requires SuiteSparse_config 7.8.0 or later"
#endif

#endif


#ifndef NDEBUG
#define NDEBUG
#endif

// To enable debugging, uncomment the following line:
// #undef NDEBUG

/* ------------------------------------------------------------------------- */
/* basic definitions */
/* ------------------------------------------------------------------------- */

#ifdef FLIP
#undef FLIP
#endif

#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

#ifdef EMPTY
#undef EMPTY
#endif

#define PRIVATE static

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < EMPTY) ? FLIP (i) : (i))

/* for integer MAX/MIN, or for doubles when we don't care how NaN's behave: */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* logical expression of p implies q: */
#define IMPLIES(p,q) (!(p) || (q))

/* Note that the IBM RS 6000 xlc predefines TRUE and FALSE in <types.h>. */
/* The Compaq Alpha also predefines TRUE and FALSE. */
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#define TRUE (1)
#define FALSE (0)
#define EMPTY (-1)

/* largest value of size_t */
#ifndef SIZE_T_MAX
#ifdef SIZE_MAX
/* C99 only */
#define SIZE_T_MAX SIZE_MAX
#else
#define SIZE_T_MAX ((size_t) (-1))
#endif
#endif

/* ------------------------------------------------------------------------- */
/* integer type for AMD: int32_t or int64_t */
/* ------------------------------------------------------------------------- */

#if defined (DLONG) || defined (ZLONG)

#define Int int64_t
#define UInt uint64_t
#define ID  "%" PRId64
#define Int_MAX INT64_MAX

#define AMD_order amd_l_order
#define AMD_defaults amd_l_defaults
#define AMD_control amd_l_control
#define AMD_info amd_l_info
#define AMD_1 amd_l1
#define AMD_2 amd_l2
#define AMD_valid amd_l_valid
#define AMD_aat amd_l_aat
#define AMD_postorder amd_l_postorder
#define AMD_post_tree amd_l_post_tree
#define AMD_dump amd_l_dump
#define AMD_debug amd_l_debug
#define AMD_debug_init amd_l_debug_init
#define AMD_preprocess amd_l_preprocess

#else

#define Int int32_t
#define UInt uint32_t
#define ID "%d"
#define Int_MAX INT32_MAX

#define AMD_defaults amd_defaults
#define AMD_control amd_control
// #define AMD_info amd_info
// #define AMD_1 amd_1
// #define AMD_2 amd_2
// #define AMD_valid amd_valid
// #define AMD_aat amd_aat
// #define AMD_postorder amd_postorder
// #define AMD_post_tree amd_post_tree
#define AMD_dump amd_dump
#define AMD_debug amd_debug
#define AMD_debug_init amd_debug_init
// #define AMD_preprocess amd_preprocess

#endif

/* ------------------------------------------------------------------------- */
/* AMD routine definitions (not user-callable) */
/* ------------------------------------------------------------------------- */

size_t AMD_aat
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],
    Int Tp [ ],
    double Info [ ]
) ;

void AMD_1
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    Int Pinv [ ],
    Int Len [ ],
    Int slen,
    Int S [ ],
    double Control [ ],
    double Info [ ]
) ;

void AMD_postorder
(
    Int nn,
    Int Parent [ ],
    Int Npiv [ ],
    Int Fsize [ ],
    Int Order [ ],
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
) ;

Int AMD_post_tree
(
    Int root,
    Int k,
    Int Child [ ],
    const Int Sibling [ ],
    Int Order [ ],
    Int Stack [ ]
#ifndef NDEBUG
    , Int nn
#endif
) ;

void AMD_preprocess
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Rp [ ],
    Int Ri [ ],
    Int W [ ],
    Int Flag [ ]
) ;

/* ------------------------------------------------------------------------- */
/* debugging definitions */
/* ------------------------------------------------------------------------- */

/* no debugging */
#define ASSERT(expression)
#define AMD_DEBUG0(params)
#define AMD_DEBUG1(params)
#define AMD_DEBUG2(params)
#define AMD_DEBUG3(params)
#define AMD_DEBUG4(params)




/* make it easy for C++ programs to include BTF */
#ifdef __cplusplus
extern "C" {
#endif

int32_t btf_maxtrans    /* returns # of columns matched */
(
    /* --- input, not modified: --- */
    int32_t nrow,   /* A is nrow-by-ncol in compressed column form */
    int32_t ncol,
    int32_t Ap [ ], /* size ncol+1 */
    int32_t Ai [ ], /* size nz = Ap [ncol] */
    double maxwork, /* maximum amount of work to do is maxwork*nnz(A); no limit
                     * if <= 0 */

    /* --- output, not defined on input --- */
    double *work,   /* work = -1 if maxwork > 0 and the total work performed
                     * reached the maximum of maxwork*nnz(A).
                     * Otherwise, work = the total work performed. */

    int32_t Match [ ], /* size nrow. Match [i] = j if column j matched to row i
                     * (see above for the singular-matrix case) */

    /* --- workspace, not defined on input or output --- */
    int32_t Work [ ]    /* size 5*ncol */
) ;

/* int64_t integer version */
int64_t btf_l_maxtrans (int64_t, int64_t,
    int64_t *, int64_t *, double, double *,
    int64_t *, int64_t *) ;


/* ========================================================================== */
/* === BTF_STRONGCOMP ======================================================= */
/* ========================================================================== */

int32_t btf_strongcomp  /* return # of strongly connected components */
(
    /* input, not modified: */
    int32_t n,      /* A is n-by-n in compressed column form */
    int32_t Ap [ ], /* size n+1 */
    int32_t Ai [ ], /* size nz = Ap [n] */

    /* optional input, modified (if present) on output: */
    int32_t Q [ ],  /* size n, input column permutation */

    /* output, not defined on input */
    int32_t P [ ],  /* size n.  P [k] = j if row and column j are kth row/col
                     * in permuted matrix. */

    int32_t R [ ],  /* size n+1.  block b is in rows/cols R[b] ... R[b+1]-1 */

    /* workspace, not defined on input or output */
    int32_t Work [ ]    /* size 4n */
) ;

int64_t btf_l_strongcomp (int64_t, int64_t *,
    int64_t *, int64_t *, int64_t *,
    int64_t *, int64_t *) ;


/* ========================================================================== */
/* === BTF_ORDER ============================================================ */
/* ========================================================================== */


int32_t btf_order       /* returns number of blocks found */
(
    /* --- input, not modified: --- */
    int32_t n,      /* A is n-by-n in compressed column form */
    int32_t Ap [ ], /* size n+1 */
    int32_t Ai [ ], /* size nz = Ap [n] */
    double maxwork, /* do at most maxwork*nnz(A) work in the maximum
                     * transversal; no limit if <= 0 */

    /* --- output, not defined on input --- */
    double *work,   /* return value from btf_maxtrans */
    int32_t P [ ],  /* size n, row permutation */
    int32_t Q [ ],  /* size n, column permutation */
    int32_t R [ ],  /* size n+1.  block b is in rows/cols R[b] ... R[b+1]-1 */
    int32_t *nmatch, /* # nonzeros on diagonal of P*A*Q */

    /* --- workspace, not defined on input or output --- */
    int32_t Work [ ] /* size 5n */
) ;

int64_t btf_l_order (int64_t, int64_t *, int64_t *, double , double *,
    int64_t *, int64_t *, int64_t *, int64_t *, int64_t *) ;

//------------------------------------------------------------------------------
// btf_version: return BTF version
//------------------------------------------------------------------------------

void btf_version (int version [3]) ;

#ifdef __cplusplus
}
#endif



#define BTF_FLIP(j) (-(j)-2)
#define BTF_ISFLIPPED(j) ((j) < -1)
#define BTF_UNFLIP(j) ((BTF_ISFLIPPED (j)) ? BTF_FLIP (j) : (j))

/* ========================================================================== */
/* === BTF version ========================================================== */
/* ========================================================================== */

#define BTF_DATE "Mar 22, 2024"
#define BTF_MAIN_VERSION   2
#define BTF_SUB_VERSION    3
#define BTF_SUBSUB_VERSION 2

#define BTF_VERSION_CODE(main,sub) SUITESPARSE_VER_CODE(main,sub)
#define BTF_VERSION BTF_VERSION_CODE(2,3)

#define BTF__VERSION SUITESPARSE__VERCODE(2,3,2)
#if !defined (SUITESPARSE__VERSION) || \
    (SUITESPARSE__VERSION < SUITESPARSE__VERCODE(7,7,0))
#error "BTF 2.3.2 requires SuiteSparse_config 7.7.0 or later"
#endif

#ifndef _BTF_INTERNAL_H
#define _BTF_INTERNAL_H

/* Not to be included in any user program. */

#ifdef DLONG
#define Int int64_t
#define Int_id "%" PRId64
#define BTF(name) btf_l_ ## name
#else
#define Int int32_t
#define Int_id "%d"
#define BTF(name) btf_ ## name
#endif

/* ========================================================================== */
/* make sure debugging and printing is turned off */

#ifndef NDEBUG
#define NDEBUG
#endif
#ifndef NPRINT
#define NPRINT
#endif

/* To enable debugging and assertions, uncomment this line: 
 #undef NDEBUG
*/
/* To enable diagnostic printing, uncomment this line: 
 #undef NPRINT
*/

/* ========================================================================== */



#undef TRUE
#undef FALSE
#undef PRINTF
#undef MIN

#ifndef NPRINT
#define PRINTF(s) { printf s ; } ;
#else
#define PRINTF(s)
#endif

#define TRUE 1
#define FALSE 0
#define EMPTY (-1)
#define MIN(a,b) (((a) < (b)) ?  (a) : (b))

#endif
