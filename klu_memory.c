#include "klu.h"


KLU_symbolic *KLU_alloc_symbolic
(
    Int n,
    Int *Ap,
    Int *Ai,
    KLU_common *Common
)
{
    KLU_symbolic *Symbolic ;
    Int *P, *Q, *R ;
    double *Lnz ;
    Int nz, i, j, p, pend ;

    if (Common == NULL)
    {
        return (NULL) ;
    }
    Common->status = KLU_OK ;

    /* A is n-by-n, with n > 0.  Ap [0] = 0 and nz = Ap [n] >= 0 required.
     * Ap [j] <= Ap [j+1] must hold for all j = 0 to n-1.  Row indices in Ai
     * must be in the range 0 to n-1, and no duplicate entries can be present.
     * The list of row indices in each column of A need not be sorted.
     */

    if (n <= 0 || Ap == NULL || Ai == NULL)
    {
        /* Ap and Ai must be present, and n must be > 0 */
        Common->status = KLU_INVALID ;
        return (NULL) ;
    }

    nz = Ap [n] ;
    if (Ap [0] != 0 || nz < 0)
    {
        /* nz must be >= 0 and Ap [0] must equal zero */
        Common->status = KLU_INVALID ;
        return (NULL) ;
    }

    for (j = 0 ; j < n ; j++)
    {
        if (Ap [j] > Ap [j+1])
        {
            /* column pointers must be non-decreasing */
            Common->status = KLU_INVALID ;
            return (NULL) ;
        }
    }
    P = KLU_malloc (n, sizeof (Int), Common) ;
    if (Common->status < KLU_OK)
    {
        /* out of memory */
        Common->status = KLU_OUT_OF_MEMORY ;
        return (NULL) ;
    }
    for (i = 0 ; i < n ; i++)
    {
        P [i] = EMPTY ;
    }
    for (j = 0 ; j < n ; j++)
    {
        pend = Ap [j+1] ;
        for (p = Ap [j] ; p < pend ; p++)
        {
            i = Ai [p] ;
            if (i < 0 || i >= n || P [i] == j)
            {
                /* row index out of range, or duplicate entry */
                KLU_free (P, n, sizeof (Int), Common) ;
                Common->status = KLU_INVALID ;
                return (NULL) ;
            }
            /* flag row i as appearing in column j */
            P [i] = j ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Symbolic = KLU_malloc (1, sizeof (KLU_symbolic), Common) ;
    if (Common->status < KLU_OK)
    {
        /* out of memory */
        KLU_free (P, n, sizeof (Int), Common) ;
        Common->status = KLU_OUT_OF_MEMORY ;
        return (NULL) ;
    }

    Q = KLU_malloc (n, sizeof (Int), Common) ;
    R = KLU_malloc (n+1, sizeof (Int), Common) ;
    Lnz = KLU_malloc (n, sizeof (double), Common) ;

    Symbolic->n = n ;
    Symbolic->nz = nz ;
    Symbolic->P = P ;
    Symbolic->Q = Q ;
    Symbolic->R = R ;
    Symbolic->Lnz = Lnz ;

    // if (Common->status < KLU_OK)  status:0    KLU_OK:0
    // {
    //     /* out of memory */
    //     KLU_free_symbolic (&Symbolic, Common) ;
    //     Common->status = KLU_OUT_OF_MEMORY ;
    //     return (NULL) ;
    // }

    return (Symbolic) ;
}


void *KLU_malloc        /* returns pointer to the newly malloc'd block */
(
    /* ---- input ---- */
    size_t n,           /* number of items */
    size_t size,        /* size of each item */
    /* --------------- */
    KLU_common *Common
)
{
    void *p ;

    if (Common == NULL)
    {
        p = NULL ;
    }
    else if (size == 0)
    {
        /* size must be > 0 */
        Common->status = KLU_INVALID ;
        p = NULL ;
    }
    else if (sizeof (size_t) > sizeof (Int) && n >= 2147483647)
    {
        /* object is too big to allocate; p[i] where i is an Int will not
         * be enough. */
        Common->status = KLU_TOO_LARGE ;
        p = NULL ;
    }
    else
    {
        /* call malloc, or its equivalent */
        p = SuiteSparse_malloc (n, size) ;
        if (p == NULL)
        {
            /* failure: out of memory */
            Common->status = KLU_OUT_OF_MEMORY ;
        }
        else
        {
            Common->memusage += (MAX (1,n) * size) ;
            Common->mempeak = MAX (Common->mempeak, Common->memusage) ;
        }
    }
    return (p) ;
}

void *SuiteSparse_malloc    /* pointer to allocated block of memory */
(
    size_t nitems,          /* number of items to malloc */
    size_t size_of_item     /* sizeof each item */
)
{
    void *p ;
    size_t size ;
    if (nitems < 1) nitems = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems * size_of_item  ;

    if (size != ((double) nitems) * size_of_item)
    {
        /* size_t overflow */
        p = NULL ;
    }
    else
    {
 //  p = (void *) (SuiteSparse_config.malloc_func) (size) ; 原这样写是方便在不同的编译环境适配
        p =  malloc (size) ;
    }
    return (p) ;
}


void *KLU_free          /* always returns NULL */
(
    /* ---- in/out --- */
    void *p,            /* block of memory to free */
    /* ---- input --- */
    size_t n,           /* size of block to free, in # of items */
    size_t size,        /* size of each item */
    /* --------------- */
    KLU_common *Common
)
{
    if (p != NULL && Common != NULL)
    {
        /* only free the object if the pointer is not NULL */
        /* call free, or its equivalent */
        SuiteSparse_free (p) ;
        Common->memusage -= (MAX (1,n) * size) ;
    }
    /* return NULL, and the caller should assign this to p.  This avoids
     * freeing the same pointer twice. */
    return (NULL) ;
}

void *SuiteSparse_free      /* always returns NULL */
(
    void *p                 /* block to free */
)
{
    if (p)
    {
        free(p) ;
    }
    return (NULL) ;
}