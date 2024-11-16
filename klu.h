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
#include <sched.h>
#include "amd_internal.h"
#include "amd.h"
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


// static void dfs
// (
//     /* inputs, not modified on output: */
//     Int j,              /* start the DFS at node j */
//     Int Ap [ ],         /* size n+1, column pointers for the matrix A */
//     Int Ai [ ],         /* row indices, size nz = Ap [n] */
//     Int Q [ ],          /* input column permutation */

//     /* inputs, modified on output (each array is of size n): */
//     Int Time [ ],       /* Time [j] = "time" that node j was first visited */
//     Int Flag [ ],       /* Flag [j]: see above */
//     Int Low [ ],        /* Low [j]: see definition below */
//     Int *p_nblocks,     /* number of blocks (aka strongly-connected-comp.)*/
//     Int *p_timestamp,   /* current "time" */

//     /* workspace, not defined on input or output: */
//     Int Cstack [ ],     /* size n, output stack to hold nodes of components */
//     Int Jstack [ ],     /* size n, stack for the variable j */
//     Int Pstack [ ]      /* size n, stack for the variable p */
// );