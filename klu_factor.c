
#include"klu.h"

/* ========================================================================== */
/* === KLU_factor ========================================================== */
/* ========================================================================== */
size_t KLU_add_size_t (size_t a, size_t b, Int *ok)
{
    (*ok) = (*ok) && ((a + b) >= MAX (a,b)) ;
    return ((*ok) ? (a + b) : ((size_t) -1)) ;
}

size_t KLU_mult_size_t (size_t a, size_t k, Int *ok)
{
    size_t i, s = 0 ;
    for (i = 0 ; i < k ; i++)
    {
        s = KLU_add_size_t (s, a, ok) ;
    }
    return ((*ok) ? s : ((size_t) -1)) ;
}


int KLU_free_numeric
(
    KLU_numeric **NumericHandle,
    KLU_common  *Common
)
{
    KLU_numeric *Numeric ;
    Unit **LUbx ;
    size_t *LUsize ;
    Int block, n, nzoff, nblocks ;

    if (Common == NULL)
    {
        return (FALSE) ;
    }
    if (NumericHandle == NULL || *NumericHandle == NULL)
    {
        return (TRUE) ;
    }

    Numeric = *NumericHandle ;

    n = Numeric->n ;
    nzoff = Numeric->nzoff ;
    nblocks = Numeric->nblocks ;
    LUsize = Numeric->LUsize ;

    LUbx = (Unit **) Numeric->LUbx ;
    if (LUbx != NULL)
    {
        for (block = 0 ; block < nblocks ; block++)
        {
            KLU_free (LUbx [block], LUsize ? LUsize [block] : 0,
                sizeof (Unit), Common) ;
        }
    }

    KLU_free (Numeric->Pnum, n, sizeof (Int), Common) ;
    KLU_free (Numeric->Offp, n+1, sizeof (Int), Common) ;
    KLU_free (Numeric->Offi, nzoff+1, sizeof (Int), Common) ;
    KLU_free (Numeric->Offx, nzoff+1, sizeof (Entry), Common) ;

    KLU_free (Numeric->Lip,  n, sizeof (Int), Common) ;
    KLU_free (Numeric->Llen, n, sizeof (Int), Common) ;
    KLU_free (Numeric->Uip,  n, sizeof (Int), Common) ;
    KLU_free (Numeric->Ulen, n, sizeof (Int), Common) ;

    KLU_free (Numeric->LUsize, nblocks, sizeof (size_t), Common) ;

    KLU_free (Numeric->LUbx, nblocks, sizeof (Unit *), Common) ;

    KLU_free (Numeric->Udiag, n, sizeof (Entry), Common) ;

    KLU_free (Numeric->Rs,   n, sizeof (double), Common) ;
    KLU_free (Numeric->Pinv, n, sizeof (Int), Common) ;

    KLU_free (Numeric->Work, Numeric->worksize, 1, Common) ;

    KLU_free (Numeric, 1, sizeof (KLU_numeric), Common) ;

    *NumericHandle = NULL ;
    return (TRUE) ;
}

void *KLU_realloc       /* returns pointer to reallocated block */
(
    /* ---- input ---- */
    size_t nnew,        /* requested # of items in reallocated block */
    size_t nold,        /* old # of items */
    size_t size,        /* size of each item */
    /* ---- in/out --- */
    void *p,            /* block of memory to realloc */
    /* --------------- */
    KLU_common *Common
)
{
    void *pnew ;
    int ok = TRUE ;

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
    else if (p == NULL)
    {
        /* A fresh object is being allocated. */
        p = KLU_malloc (nnew, size, Common) ;
    }
    else if (sizeof (size_t) > sizeof (Int) && nnew >= Int_MAX)
    {
        /* failure: nnew is too big.  Do not change p */
        Common->status = KLU_TOO_LARGE ;
    }
    else
    {
        /* The object exists, and is changing to some other nonzero size. */
        /* call realloc, or its equivalent */
        pnew = SuiteSparse_realloc (nnew, nold, size, p, &ok) ;
        if (ok)
        {
            /* success: return the new p and change the size of the block */
            Common->memusage += ((nnew-nold) * size) ;
            Common->mempeak = MAX (Common->mempeak, Common->memusage) ;
            p = pnew ;
        }
        else
        {
            /* Do not change p, since it still points to allocated memory */
            Common->status = KLU_OUT_OF_MEMORY ;
        }
    }
    return (p) ;
}


static Int dfss
(
    /* input, not modified on output: */
    Int j,              /* node at which to start the DFS */
    Int k,              /* mark value, for the Flag array */
    Int Pinv [ ],       /* Pinv [i] = k if row i is kth pivot row, or EMPTY if
                         * row i is not yet pivotal.  */
    Int Llen [ ],       /* size n, Llen [k] = # nonzeros in column k of L */
    Int Lip [ ],        /* size n, Lip [k] is position in LU of column k of L */

    /* workspace, not defined on input or output */
    Int Stack [ ],      /* size n */

    /* input/output: */
    Int Flag [ ],       /* Flag [i] == k means i is marked */
    Int Lpend [ ],      /* for symmetric pruning */
    Int top,            /* top of stack on input*/
    Unit LU [],
    Int *Lik,           /* Li row index array of the kth column */
    Int *plength,

    /* other, not defined on input or output */
    Int Ap_pos [ ]      /* keeps track of position in adj list during DFS */
)
{
    Int i, pos, jnew, head, l_length ;
    Int *Li ;

    l_length = *plength ;

    head = 0 ;
    Stack [0] = j ;
    ASSERT (Flag [j] != k) ;

    while (head >= 0)
    {
        j = Stack [head] ;
        jnew = Pinv [j] ;
        ASSERT (jnew >= 0 && jnew < k) ;        /* j is pivotal */

        if (Flag [j] != k)          /* a node is not yet visited */
        {
            /* first time that j has been visited */
            Flag [j] = k ;
            /* set Ap_pos [head] to one past the last entry in col j to scan */
            Ap_pos [head] =
                (Lpend [jnew] == EMPTY) ?  Llen [jnew] : Lpend [jnew] ;
        }

        /* add the adjacent nodes to the recursive stack by iterating through
         * until finding another non-visited pivotal node */
        Li = (Int *) (LU + Lip [jnew]) ;
        for (pos = --Ap_pos [head] ; pos >= 0 ; --pos)
        {
            i = Li [pos] ;
            if (Flag [i] != k)
            {
                /* node i is not yet visited */
                if (Pinv [i] >= 0)
                {
                    /* keep track of where we left off in the scan of the
                     * adjacency list of node j so we can restart j where we
                     * left off. */
                    Ap_pos [head] = pos ;

                    /* node i is pivotal; push it onto the recursive stack
                     * and immediately break so we can recurse on node i. */
                    Stack [++head] = i ;
                    break ;
                }
                else
                {
                    /* node i is not pivotal (no outgoing edges). */
                    /* Flag as visited and store directly into L,
                     * and continue with current node j. */
                    Flag [i] = k ;
                    Lik [l_length] = i ;
                    l_length++ ;
                }
            }
        }

        if (pos == -1)
        {
            /* if all adjacent nodes of j are already visited, pop j from
             * recursive stack and push j onto output stack */
            head-- ;
            Stack[--top] = j ;
        }
    }

    *plength = l_length ;
    return (top) ;
}

//确定下三角矩阵的非零模式
static Int lsolve_symbolic
(
    /* input, not modified on output: */
    Int n,              /* L is n-by-n, where n >= 0 */
    Int k,              /* also used as the mark value, for the Flag array */
    Int Ap [ ],
    Int Ai [ ],
    Int Q [ ],
    Int Pinv [ ],       /* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
                         * is not yet pivotal.  */

    /* workspace, not defined on input or output */
    Int Stack [ ],      /* size n */

    /* workspace, defined on input and output */
    Int Flag [ ],       /* size n.  Initially, all of Flag [0..n-1] < k.  After
                         * lsolve_symbolic is done, Flag [i] == k if i is in
                         * the pattern of the output, and Flag [0..n-1] <= k. */

    /* other */
    Int Lpend [ ],      /* for symmetric pruning */
    Int Ap_pos [ ],     /* workspace used in dfs */

    Unit LU [ ],        /* LU factors (pattern and values) */
    Int lup,            /* pointer to free space in LU */
    Int Llen [ ],       /* size n, Llen [k] = # nonzeros in column k of L */
    Int Lip [ ],        /* size n, Lip [k] is position in LU of column k of L */

    /* ---- the following are only used in the BTF case --- */

    Int k1,             /* the block of A is from k1 to k2-1 */
    Int PSinv [ ]       /* inverse of P from symbolic factorization */
)
{
    Int *Lik ;
    Int i, p, pend, oldcol, kglobal, top, l_length ;

    top = n ;
    l_length = 0 ;
    Lik = (Int *) (LU + lup);

    /* ---------------------------------------------------------------------- */
    /* BTF factorization of A (k1:k2-1, k1:k2-1) */
    /* ---------------------------------------------------------------------- */

    kglobal = k + k1 ;  /* column k of the block is col kglobal of A */
    oldcol = Q [kglobal] ;      /* Q must be present for BTF case */
    pend = Ap [oldcol+1] ;
    for (p = Ap [oldcol] ; p < pend ; p++)
    {
        i = PSinv [Ai [p]] - k1 ;
        if (i < 0) continue ;   /* skip entry outside the block */

        /* (i,k) is an entry in the block.  start a DFS at node i */
        if (Flag [i] != k)
        {
            if (Pinv [i] >= 0)
            {
                //这是 static int dfs 而不是 void dfs 这里把int dfs弄成int dfs
                top = dfss (i, k, Pinv, Llen, Lip, Stack, Flag,
                           Lpend, top, LU, Lik, &l_length, Ap_pos) ;
            }
            else
            {
                /* i is not pivotal, and not flagged. Flag and put in L */
                Flag [i] = k ;
                Lik [l_length] = i ;
                l_length++;
            }
        }
    }

    /* If Llen [k] is zero, the matrix is structurally singular */
    Llen [k] = l_length ;
    return (top) ;
}

void *SuiteSparse_realloc   /* pointer to reallocated block of memory, or
                               to original block if the realloc failed. */
(
    size_t nitems_new,      /* new number of items in the object */
    size_t nitems_old,      /* old number of items in the object */
    size_t size_of_item,    /* sizeof each item */
    void *p,                /* old object to reallocate */
    int *ok                 /* 1 if successful, 0 otherwise */
)
{
    size_t size ;
    if (nitems_old < 1) nitems_old = 1 ;
    if (nitems_new < 1) nitems_new = 1 ;
    if (size_of_item < 1) size_of_item = 1 ;
    size = nitems_new * size_of_item  ;

    if (size != ((double) nitems_new) * size_of_item)
    {
        /* size_t overflow */
        (*ok) = 0 ;
    }
    else if (p == NULL)
    {
        /* a fresh object is being allocated */
        p = SuiteSparse_malloc (nitems_new, size_of_item) ;
        (*ok) = (p != NULL) ;
    }
    else if (nitems_old == nitems_new)
    {
        /* the object does not change; do nothing */
        (*ok) = 1 ;
    }
    else
    {
        /* change the size of the object from nitems_old to nitems_new */
        void *pnew ;
        pnew = realloc (p, size) ;
        if (pnew == NULL)
        {
            if (nitems_new < nitems_old)
            {
                /* the attempt to reduce the size of the block failed, but
                   the old block is unchanged.  So pretend to succeed. */
                (*ok) = 1 ;
            }
            else
            {
                /* out of memory */
                (*ok) = 0 ;
            }
        }
        else
        {
            /* success */
            p = pnew ;
            (*ok) = 1 ;
        }
    }
    return (p) ;
}


static void construct_column
(
    /* inputs, not modified on output */
    Int k,          /* the column of A (or the column of the block) to get */
    Int Ap [ ],
    Int Ai [ ],
    Entry Ax [ ],
    Int Q [ ],      /* column pre-ordering */

    /* zero on input, modified on output */
    Entry X [ ],

    /* ---- the following are only used in the BTF case --- */

    /* inputs, not modified on output */
    Int k1,         /* the block of A is from k1 to k2-1 */
    Int PSinv [ ],  /* inverse of P from symbolic factorization */
    double Rs [ ],  /* scale factors for A */
    Int scale,      /* 0: no scaling, nonzero: scale the rows with Rs */

    /* inputs, modified on output */
    Int Offp [ ],   /* off-diagonal matrix (modified by this routine) */
    Int Offi [ ],
    Entry Offx [ ]
)
{
    Entry aik ;
    Int i, p, pend, oldcol, kglobal, poff, oldrow ;

    /* ---------------------------------------------------------------------- */
    /* Scale and scatter the column into X. */
    /* ---------------------------------------------------------------------- */

    kglobal = k + k1 ;          /* column k of the block is col kglobal of A */
    poff = Offp [kglobal] ;     /* start of off-diagonal column */
    oldcol = Q [kglobal] ;
    pend = Ap [oldcol+1] ;

    if (scale <= 0)
    {
        /* no scaling */
        for (p = Ap [oldcol] ; p < pend ; p++)
        {
            oldrow = Ai [p] ;
            i = PSinv [oldrow] - k1 ;
            aik = Ax [p] ;
            if (i < 0)
            {
                /* this is an entry in the off-diagonal part */
                Offi [poff] = oldrow ;
                Offx [poff] = aik ;
                poff++ ;
            }
            else
            {
                /* (i,k) is an entry in the block.  scatter into X */
                X [i] = aik ;
            }
        }
    }
    else
    {
        /* row scaling */
        for (p = Ap [oldcol] ; p < pend ; p++)
        {
            oldrow = Ai [p] ;
            i = PSinv [oldrow] - k1 ;
            aik = Ax [p] ;
            SCALE_DIV (aik, Rs [oldrow]) ;
            if (i < 0)
            {
                /* this is an entry in the off-diagonal part */
                Offi [poff] = oldrow ;
                Offx [poff] = aik ;
                poff++ ;
            }
            else
            {
                /* (i,k) is an entry in the block.  scatter into X */
                X [i] = aik ;
            }
        }
    }

    Offp [kglobal+1] = poff ;   /* start of the next col of off-diag part */
}

static void lsolve_numeric
(
    /* input, not modified on output: */
    Int Pinv [ ],       /* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
                         * is not yet pivotal.  */
    Unit *LU,           /* LU factors (pattern and values) */
    Int Stack [ ],      /* stack for dfs */
    Int Lip [ ],        /* size n, Lip [k] is position in LU of column k of L */
    Int top,            /* top of stack on input */
    Int n,              /* A is n-by-n */
    Int Llen [ ],       /* size n, Llen [k] = # nonzeros in column k of L */

    /* output, must be zero on input: */
    Entry X [ ] /* size n, initially zero.  On output,
                 * X [Ui [up1..up-1]] and X [Li [lp1..lp-1]]
                 * contains the solution. */

)
{
    Entry xj ;
    Entry *Lx ;
    Int *Li ;
    Int p, s, j, jnew, len ;

    /* solve Lx=b */
    for (s = top ; s < n ; s++)
    {
        /* forward solve with column j of L */
        j = Stack [s] ;
        jnew = Pinv [j] ;
        ASSERT (jnew >= 0) ;
        xj = X [j] ;
        GET_POINTER (LU, Lip, Llen, Li, Lx, jnew, len) ;
        ASSERT (Lip [jnew] <= Lip [jnew+1]) ;
        for (p = 0 ; p < len ; p++)
        {
            /*X [Li [p]] -= Lx [p] * xj ; */
            MULT_SUB (X [Li [p]], Lx [p], xj) ;
        }
    }
}

static Int lpivot
(
    Int diagrow,
    Int *p_pivrow,
    Entry *p_pivot,
    double *p_abs_pivot,
    double tol,
    Entry X [ ],
    Unit *LU,           /* LU factors (pattern and values) */
    Int Lip [ ],
    Int Llen [ ],
    Int k,
    Int n,

    Int Pinv [ ],       /* Pinv [i] = k if row i is kth pivot row, or EMPTY if
                         * row i is not yet pivotal.  */

    Int *p_firstrow,
    KLU_common *Common
)
{
    Entry x, pivot, *Lx ;
    double abs_pivot, xabs ;
    Int p, i, ppivrow, pdiag, pivrow, *Li, last_row_index, firstrow, len ;

    pivrow = EMPTY ;
    if (Llen [k] == 0)
    {
        /* matrix is structurally singular */
        if (Common->halt_if_singular)
        {
            return (FALSE) ;
        }
        for (firstrow = *p_firstrow ; firstrow < n ; firstrow++)
        {
            if (Pinv [firstrow] < 0)
            {
                /* found the lowest-numbered non-pivotal row.  Pick it. */
                pivrow = firstrow ;
                break ;
            }
        }
        ASSERT (pivrow >= 0 && pivrow < n) ;
        CLEAR (pivot) ;
        *p_pivrow = pivrow ;
        *p_pivot = pivot ;
        *p_abs_pivot = 0 ;
        *p_firstrow = firstrow ;
        return (FALSE) ;
    }

    pdiag = EMPTY ;
    ppivrow = EMPTY ;
    abs_pivot = EMPTY ;
    i = Llen [k] - 1 ;
    GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
    last_row_index = Li [i] ;

    /* decrement the length by 1 */
    Llen [k] = i ;
    GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;

    /* look in Li [0 ..Llen [k] - 1 ] for a pivot row */
    for (p = 0 ; p < len ; p++)
    {
        /* gather the entry from X and store in L */
        i = Li [p] ;
        x = X [i] ;
        CLEAR (X [i]) ;

        Lx [p] = x ;
        /* xabs = ABS (x) ; */
        ABS (xabs, x) ;

        /* find the diagonal */
        if (i == diagrow)
        {
            pdiag = p ;
        }

        /* find the partial-pivoting choice */
        if (xabs > abs_pivot)
        {
            abs_pivot = xabs ;
            ppivrow = p ;
        }
    }

    /* xabs = ABS (X [last_row_index]) ;*/
    ABS (xabs, X [last_row_index]) ;
    if (xabs > abs_pivot)
    {
        abs_pivot = xabs ;
        ppivrow = EMPTY ;
    }

    /* compare the diagonal with the largest entry */
    if (last_row_index == diagrow)
    {
        if (xabs >= tol * abs_pivot)
        {
            abs_pivot = xabs ;
            ppivrow = EMPTY ;
        }
    }
    else if (pdiag != EMPTY)
    {
        /* xabs = ABS (Lx [pdiag]) ;*/
        ABS (xabs, Lx [pdiag]) ;
        if (xabs >= tol * abs_pivot)
        {
            /* the diagonal is large enough */
            abs_pivot = xabs ;
            ppivrow = pdiag ;
        }
    }

    if (ppivrow != EMPTY)
    {
        pivrow = Li [ppivrow] ;
        pivot  = Lx [ppivrow] ;
        /* overwrite the ppivrow values with last index values */
        Li [ppivrow] = last_row_index ;
        Lx [ppivrow] = X [last_row_index] ;
    }
    else
    {
        pivrow = last_row_index ;
        pivot = X [last_row_index] ;
    }
    CLEAR (X [last_row_index]) ;

    *p_pivrow = pivrow ;
    *p_pivot = pivot ;
    *p_abs_pivot = abs_pivot ;
    ASSERT (pivrow >= 0 && pivrow < n) ;

    if (IS_ZERO (pivot) && Common->halt_if_singular)
    {
        /* numerically singular case */
        return (FALSE) ;
    }

    /* divide L by the pivot value */
    for (p = 0 ; p < Llen [k] ; p++)
    {
        /* Lx [p] /= pivot ; */
        DIV (Lx [p], Lx [p], pivot) ;
    }

    return (TRUE) ;
}

static void prune
(
    /* input/output: */
    Int Lpend [ ],      /* Lpend [j] marks symmetric pruning point for L(:,j) */

    /* input: */
    Int Pinv [ ],       /* Pinv [i] = k if row i is kth pivot row, or EMPTY if
                         * row i is not yet pivotal.  */
    Int k,              /* prune using column k of U */
    Int pivrow,         /* current pivot row */

    /* input/output: */
    Unit *LU,           /* LU factors (pattern and values) */

    /* input */
    Int Uip [ ],        /* size n, column pointers for U */
    Int Lip [ ],        /* size n, column pointers for L */
    Int Ulen [ ],       /* size n, column length of U */
    Int Llen [ ]        /* size n, column length of L */
)
{
    Entry x ;
    Entry *Lx, *Ux ;
    Int *Li, *Ui ;
    Int p, i, j, p2, phead, ptail, llen, ulen ;

    /* check to see if any column of L can be pruned */
    /* Ux is set but not used.  This OK. */
    GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
    for (p = 0 ; p < ulen ; p++)
    {
        j = Ui [p] ;
        ASSERT (j < k) ;
        if (Lpend [j] == EMPTY)
        {
            /* scan column j of L for the pivot row */
            GET_POINTER (LU, Lip, Llen, Li, Lx, j, llen) ;
            for (p2 = 0 ; p2 < llen ; p2++)
            {
                if (pivrow == Li [p2])
                {
                    /* found it!  This column can be pruned */
                    /* partition column j of L.  The unit diagonal of L
                     * is not stored in the column of L. */
                    phead = 0 ;
                    ptail = Llen [j] ;
                    while (phead < ptail)
                    {
                        i = Li [phead] ;
                        if (Pinv [i] >= 0)
                        {
                            /* leave at the head */
                            phead++ ;
                        }
                        else
                        {
                            /* swap with the tail */
                            ptail-- ;
                            Li [phead] = Li [ptail] ;
                            Li [ptail] = i ;
                            x = Lx [phead] ;
                            Lx [phead] = Lx [ptail] ;
                            Lx [ptail] = x ;
                        }
                    }

                    /* set Lpend to one past the last entry in the
                     * first part of the column of L.  Entries in
                     * Li [0 ... Lpend [j]-1] are the only part of
                     * column j of L that needs to be scanned in the DFS.
                     * Lpend [j] was EMPTY; setting it >= 0 also flags
                     * column j as pruned. */
                    Lpend [j] = ptail ;

                    break ;
                }
            }
        }
    }
}

//指定缩放方式，0 表示不进行缩放，1 表示按行求和缩放，2 表示按行最大值缩放
int KLU_scale           /* return TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    int scale,          /* 0: none, 1: sum, 2: max */
    Int n,
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    double Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ],      /* size n, can be NULL if scale <= 0 */
    /* workspace, not defined on input or output */
    Int W [ ],          /* size n, can be NULL */
    /* --------------- */
    KLU_common *Common
)
{
    double a ;
    Entry *Az ;
    Int row, col, p, pend, check_duplicates ;

    /* ---------------------------------------------------------------------- */
    /* 检查各种输入的情况*/
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE) ;
    }
    Common->status = KLU_OK ;

    if (scale < 0)
    {
        /* return without checking anything and without computing the
         * scale factors */
        return (TRUE) ;
    }

    Az = (Entry *) Ax ;

    if (n <= 0 || Ap == NULL || Ai == NULL || Az == NULL ||
        (scale > 0 && Rs == NULL))
    {
        /* Ap, Ai, Ax and Rs must be present, and n must be > 0 */
        Common->status = KLU_INVALID ;
        return (FALSE) ;
    }
    if (Ap [0] != 0 || Ap [n] < 0)
    {
        /* nz = Ap [n] must be >= 0 and Ap [0] must equal zero */
        Common->status = KLU_INVALID ;
        return (FALSE) ;
    }
    for (col = 0 ; col < n ; col++)
    {
        if (Ap [col] > Ap [col+1])
        {
            /* column pointers must be non-decreasing */
            Common->status = KLU_INVALID ;
            return (FALSE) ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* scale  进行缩放*/
    /* ---------------------------------------------------------------------- */

    if (scale > 0)
    {
        /* 初始化存放缩放因子数组 为0 */
        for (row = 0 ; row < n ; row++)
        {
            Rs [row] = 0 ;
        }
    }

    /* 检查工作空间 */
    check_duplicates = (W != (Int *) NULL) ;
    if (check_duplicates)
    {
        for (row = 0 ; row < n ; row++)
        {
            W [row] = EMPTY ;
        }
    }

	//双循环
    for (col = 0 ; col < n ; col++)
    {
        pend = Ap [col+1] ;
        for (p = Ap [col] ; p < pend ; p++)
        {
            row = Ai [p] ;
            if (row < 0 || row >= n)
            {
                /* row index out of range, or duplicate entry */
                Common->status = KLU_INVALID ;
                return (FALSE) ;
            }
			//检查重复项
            if (check_duplicates)
            {
                if (W [row] == col)
                {
                    /* duplicate entry */
                    Common->status = KLU_INVALID ;
                    return (FALSE) ;
                }
                /* flag row i as appearing in column col */
                W [row] = col ;
            }
            /* a = ABS (Az [p]) ;*/
            ABS (a, Az [p]) ;
            if (scale == 1) //为1 按行求和缩放
            {
                /* accumulate the abs. row sum */
                Rs [row] += a ;
            }
            else if (scale > 1) // >1 按最大值缩放
            {
                /* find the max abs. value in the row */
                Rs [row] = MAX (Rs [row], a) ;
            }
        }
    }

    if (scale > 0)
    {
        /* do not scale empty rows */
        for (row = 0 ; row < n ; row++)
        {
        

            if (Rs [row] == 0.0)
            {
                Rs [row] = 1.0 ;
            }
        }
    }

    return (TRUE) ;
}

size_t KLU_kernel   /* final size of LU on output */
(
    /* input, not modified */
    Int n,          /* A is n-by-n */
    Int Ap [ ],     /* size n+1, column pointers for A */
    Int Ai [ ],     /* size nz = Ap [n], row indices for A */
    Entry Ax [ ],   /* size nz, values of A */
    Int Q [ ],      /* size n, optional input permutation */
    size_t lusize,  /* initial size of LU on input */

    /* output, not defined on input */
    Int Pinv [ ],   /* size n, inverse row permutation, where Pinv [i] = k if
                     * row i is the kth pivot row */
    Int P [ ],      /* size n, row permutation, where P [k] = i if row i is the
                     * kth pivot row. */
    Unit **p_LU,        /* LU array, size lusize on input */
    Entry Udiag [ ],    /* size n, diagonal of U */
    Int Llen [ ],       /* size n, column length of L */
    Int Ulen [ ],       /* size n, column length of U */
    Int Lip [ ],        /* size n, column pointers for L */
    Int Uip [ ],        /* size n, column pointers for U */
    Int *lnz,           /* size of L*/
    Int *unz,           /* size of U*/
    /* workspace, not defined on input */
    Entry X [ ],    /* size n, undefined on input, zero on output */

    /* workspace, not defined on input or output */
    Int Stack [ ],  /* size n */
    Int Flag [ ],   /* size n */
    Int Ap_pos [ ],     /* size n */

    /* other workspace: */
    Int Lpend [ ],                  /* size n workspace, for pruning only */

    /* inputs, not modified on output */
    Int k1,             /* the block of A is from k1 to k2-1 */
    Int PSinv [ ],      /* inverse of P from symbolic factorization */
    double Rs [ ],      /* scale factors for A */

    /* inputs, modified on output */
    Int Offp [ ],   /* off-diagonal matrix (modified by this routine) */
    Int Offi [ ],
    Entry Offx [ ],
    /* --------------- */
    KLU_common *Common
)
{
    Entry pivot ;
    double abs_pivot, xsize, nunits, tol, memgrow ;
    Entry *Ux ;
    Int *Li, *Ui ;
    Unit *LU ;          /* LU factors (pattern and values) */
    Int k, p, i, j, pivrow = 0, kbar, diagrow, firstrow, lup, top, scale, len ;
    size_t newlusize ;

    ASSERT (Common != NULL) ;
    scale = Common->scale ;
    tol = Common->tol ;
    memgrow = Common->memgrow ;
    *lnz = 0 ;
    *unz = 0 ;
    CLEAR (pivot) ;

    /* ---------------------------------------------------------------------- */
    /* get initial Li, Lx, Ui, and Ux */
    /* ---------------------------------------------------------------------- */

    ASSERT (lusize > 0) ;
    LU = *p_LU ;

    /* ---------------------------------------------------------------------- */
    /* 初始化 */
    /* ---------------------------------------------------------------------- */

    firstrow = 0 ;
    lup = 0 ;

    for (k = 0 ; k < n ; k++)
    {
        /* X [k] = 0 ; */
        CLEAR (X [k]) ;
        Flag [k] = EMPTY ;
        Lpend [k] = EMPTY ;     /* flag k as not pruned */
    }

    /* ---------------------------------------------------------------------- */
    /* 标记所有行为非枢纽，并确定初始对角映射 */
    /* ---------------------------------------------------------------------- */

    /* PSinv does the symmetric permutation, so don't do it here */
    for (k = 0 ; k < n ; k++)
    {
        P [k] = k ;
        Pinv [k] = FLIP (k) ;   /* mark all rows as non-pivotal */
    }
    /*初始化非对角矩阵的构造*/
    Offp [0] = 0 ;

    /* P [k] = row means that UNFLIP (Pinv [row]) = k, and visa versa.
     * If row is pivotal, then Pinv [row] >= 0.  A row is initially "flipped"
     * (Pinv [k] < EMPTY), and then marked "unflipped" when it becomes
     * pivotal. */

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n ; k++)
    {

        /* ------------------------------------------------------------------ */
        /*确定内存是否足够*/
        /* ------------------------------------------------------------------ */

        /* (n - k) entries for L and k entries for U */
		//首先计算当前状态下下三角矩阵L和上三角矩阵U所需的内存单元数量nunits，根据剩余的行数和当前列的索引计算的整数单元和Entry类型单元的总和。
        nunits = DUNITS (Int, n - k) + DUNITS (Int, k) +
                 DUNITS (Entry, n - k) + DUNITS (Entry, k) ;

        xsize = ((double) lup) + nunits ; //计算已经使用的内存大小 如果超过了则进行内存增长
		
		//如果内存增长后仍然超出整数表示范围设置状态为KLU_TOO_LARGE并返回当前的内存大小lusize。否则，重新分配内存给LU，更新相关状态和指针
        if (xsize > (double) lusize)
        {
            /* check here how much to grow */
            xsize = (memgrow * ((double) lusize) + 4*n + 1) ;
            if (INT_OVERFLOW (xsize))
            {
                Common->status = KLU_TOO_LARGE ;
                return (lusize) ;
            }
            newlusize = memgrow * lusize + 2*n + 1 ;
            /* Future work: retry mechanism in case of malloc failure */
            LU = KLU_realloc (newlusize, lusize, sizeof (Unit), LU, Common) ;
            Common->nrealloc++ ;
            *p_LU = LU ;
            if (Common->status == KLU_OUT_OF_MEMORY)
            {
                return (lusize) ;
            }
            lusize = newlusize ;
        }

        /* ------------------------------------------------------------------ */
        /* 处理下三角和上三角矩阵的当前列 */
        /* ------------------------------------------------------------------ */

        Lip [k] = lup ;

        /* ------------------------------------------------------------------ */
        /* 计算L和U的第k列的非零模式 */
        /* ------------------------------------------------------------------ */

        top = lsolve_symbolic (n, k, Ap, Ai, Q, Pinv, Stack, Flag,
                    Lpend, Ap_pos, LU, lup, Llen, Lip, k1, PSinv) ;

        /* ------------------------------------------------------------------ */
        /* 构建要分解的矩阵的当前列，并将其分散到工作数组X中*/
        /* ------------------------------------------------------------------ */

        construct_column (k, Ap, Ai, Ax, Q, X,
            k1, PSinv, Rs, scale, Offp, Offi, Offx) ;

        /* ------------------------------------------------------------------ */
        /* 进行数值求解，计算L\A(:,k)，即下三角矩阵作用于当前列的结果 (s = L \ A (:,k)) */
        /* ------------------------------------------------------------------ */

        lsolve_numeric (Pinv, LU, Stack, Lip, top, n, Llen, X) ;

        /* ------------------------------------------------------------------ */
        /* 部分选主元 */
        /* ------------------------------------------------------------------ */

        /* 确定当前的对角行 */
        diagrow = P [k] ;   /* might already be pivotal */

        /* 调用lpivot函数寻找主元并进行列缩放。如果找不到合适的主元，设置矩阵为奇异状态，并根据需要决定是否停止因子分解。*/
        if (!lpivot (diagrow, &pivrow, &pivot, &abs_pivot, tol, X, LU, Lip,
                    Llen, k, n, Pinv, &firstrow, Common))
        {
            /* matrix is structurally or numerically singular */
            Common->status = KLU_SINGULAR ;
            if (Common->numerical_rank == EMPTY)
            {
                Common->numerical_rank = k+k1 ;
                Common->singular_col = Q [k+k1] ;
            }
            if (Common->halt_if_singular)
            {
                /* do not continue the factorization */
                return (lusize) ;
            }
        }


        /* 在确定了有效的主元行后，进行了一系列操作来更新上三角矩阵的相关信息、记录主元置换以及处理可能的非对角主元情况。 */
        ASSERT (pivrow >= 0 && pivrow < n) ;
        ASSERT (Pinv [pivrow] < 0) ;


		/*准确地设置了上三角矩阵的指针、提取了上三角矩阵的数据，并更新了相关的全局指针和存储了上三角矩阵的对角元素值，*/
        /* set the Uip pointer */
        Uip [k] = Lip [k] + UNITS (Int, Llen [k]) + UNITS (Entry, Llen [k]) ;

        /* move the lup pointer to the position where indices of U
         * should be stored */
        lup += UNITS (Int, Llen [k]) + UNITS (Entry, Llen [k]) ;

        Ulen [k] = n - top ;

        /* extract Stack [top..n-1] to Ui and the values to Ux and clear X */
        GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, len) ;
        for (p = top, i = 0 ; p < n ; p++, i++)
        {
            j = Stack [p] ;
            Ui [i] = Pinv [j] ;
            Ux [i] = X [j] ;
            CLEAR (X [j]) ;
        }

        /* position the lu index at the starting point for next column */
        lup += UNITS (Int, Ulen [k]) + UNITS (Entry, Ulen [k]) ;

        /* U(k,k) = pivot */
        Udiag [k] = pivot ;

        /* ------------------------------------------------------------------ */
        /*记录主元置换 */
        /* ------------------------------------------------------------------ */

        ASSERT (UNFLIP (Pinv [diagrow]) < n) ;
        ASSERT (P [UNFLIP (Pinv [diagrow])] == diagrow) ;

        if (pivrow != diagrow)
        {
            /* an off-diagonal pivot has been chosen */
            Common->noffdiag++ ;
            if (Pinv [diagrow] < 0)
            {
                /* the former diagonal row index, diagrow, has not yet been
                 * chosen as a pivot row.  Log this diagrow as the "diagonal"
                 * entry in the column kbar for which the chosen pivot row,
                 * pivrow, was originally logged as the "diagonal" */
                kbar = FLIP (Pinv [pivrow]) ;
                P [kbar] = diagrow ;
                Pinv [diagrow] = FLIP (kbar) ;
            }
        }
        P [k] = pivrow ;
        Pinv [pivrow] = k ;


        /* ------------------------------------------------------------------ */
        /* symmetric pruning */
        /* ------------------------------------------------------------------ */
//删减L中的列，以减少后续深度优先搜索的工作量
        prune (Lpend, Pinv, k, pivrow, LU, Uip, Lip, Ulen, Llen) ;

        *lnz += Llen [k] + 1 ; /* 1 added to lnz for diagonal */
        *unz += Ulen [k] + 1 ; /* 1 added to unz for diagonal */
    }

    /* ---------------------------------------------------------------------- */
   /*最终确定L和U的列指针，并将L放在枢纽顺序*/
    /* ---------------------------------------------------------------------- */

    for (p = 0 ; p < n ; p++)
    {
        Li = (Int *) (LU + Lip [p]) ;
        for (i = 0 ; i < Llen [p] ; i++)
        {
            Li [i] = Pinv [Li [i]] ;
        }
    }


    /* ---------------------------------------------------------------------- */
    /*将LU因子缩小到所需的大小*/
    /* ---------------------------------------------------------------------- */

    newlusize = lup ;
    ASSERT ((size_t) newlusize <= lusize) ;

    /* this cannot fail, since the block is descreasing in size */
    LU = KLU_realloc (newlusize, lusize, sizeof (Unit), LU, Common) ;
    *p_LU = LU ;
    return (newlusize) ;
}


size_t KLU_kernel_factor            /* 0 if failure, size of LU if OK */
(
    /* inputs, not modified */
    Int n,          /* A is n-by-n. n must be > 0. */
    Int Ap [ ],     /* size n+1, column pointers for A */
    Int Ai [ ],     /* size nz = Ap [n], row indices for A */
    Entry Ax [ ],   /* size nz, values of A */
    Int Q [ ],      /* size n, optional column permutation */
    double Lsize,   /* estimate of number of nonzeros in L */

    /* outputs, not defined on input */
    Unit **p_LU,        /* row indices and values of L and U */
    Entry Udiag [ ],    /* size n, diagonal of U */
    Int Llen [ ],       /* size n, column length of L */
    Int Ulen [ ],       /* size n, column length of U */
    Int Lip [ ],        /* size n, column pointers for L */
    Int Uip [ ],        /* size n, column pointers for U */
    Int P [ ],          /* row permutation, size n */
    Int *lnz,           /* size of L */
    Int *unz,           /* size of U */

    /* workspace, undefined on input */
    Entry *X,       /* size n double's, zero on output */
    Int *Work,      /* size 5n Int's */

    /* inputs, not modified on output */
    Int k1,             /* the block of A is from k1 to k2-1 */
    Int PSinv [ ],      /* inverse of P from symbolic factorization */
    double Rs [ ],      /* scale factors for A */

    /* inputs, modified on output */
    Int Offp [ ],   /* off-diagonal matrix (modified by this routine) */
    Int Offi [ ],
    Entry Offx [ ],
    /* --------------- */
    KLU_common *Common
)
{
    double maxlnz, dunits ;
    Unit *LU ;
    Int *Pinv, *Lpend, *Stack, *Flag, *Ap_pos, *W ;
    Int lsize, usize, anz, ok ;
    size_t lusize ;
    ASSERT (Common != NULL) ;

    /* ---------------------------------------------------------------------- */
    /* 获取控制参数，或使用默认值 参数初始化并检查是否溢出 */
    /* ---------------------------------------------------------------------- */

    n = MAX (1, n) ;
    anz = Ap [n+k1] - Ap [k1] ;

    if (Lsize <= 0)
    {
        Lsize = -Lsize ;
        Lsize = MAX (Lsize, 1.0) ;
        lsize = Lsize * anz + n ;
    }
    else
    {
        lsize = Lsize ;
    }

    usize = lsize ;

    lsize  = MAX (n+1, lsize) ;
    usize  = MAX (n+1, usize) ;

    maxlnz = (((double) n) * ((double) n) + ((double) n)) / 2. ;
    maxlnz = MIN (maxlnz, ((double) Int_MAX)) ;
    lsize  = MIN (maxlnz, lsize) ;
    usize  = MIN (maxlnz, usize) ;

    /* ---------------------------------------------------------------------- */
    /* 分配工作空间和输出空间的指针 */
    /* ---------------------------------------------------------------------- */

    /* return arguments are not yet assigned */
    *p_LU = (Unit *) NULL ;

    /* these computations are safe from size_t overflow */
    W = Work ;
    Pinv = (Int *) W ;      W += n ;
    Stack = (Int *) W ;     W += n ;
    Flag = (Int *) W ;      W += n ;
    Lpend = (Int *) W ;     W += n ;
    Ap_pos = (Int *) W ;    W += n ;

    dunits = DUNITS (Int, lsize) + DUNITS (Entry, lsize) +
             DUNITS (Int, usize) + DUNITS (Entry, usize) ;
    lusize = (size_t) dunits ;
    ok = !INT_OVERFLOW (dunits) ; 
    LU = ok ? KLU_malloc (lusize, sizeof (Unit), Common) : NULL ; //分配一个内存空间
    if (LU == NULL)
    {
        /* out of memory, or problem too large */
        Common->status = KLU_OUT_OF_MEMORY ;
        lusize = 0 ;
        return (lusize) ;
    }

    /* ---------------------------------------------------------------------- */
    /* 因子分解 */
    /* ---------------------------------------------------------------------- */

    /* with pruning, and non-recursive depth-first-search */
    lusize = KLU_kernel (n, Ap, Ai, Ax, Q, lusize,
            Pinv, P, &LU, Udiag, Llen, Ulen, Lip, Uip, lnz, unz,
            X, Stack, Flag, Ap_pos, Lpend,
            k1, PSinv, Rs, Offp, Offi, Offx, Common) ;

    /* ---------------------------------------------------------------------- */
    /* 返回结果 */
    /* ---------------------------------------------------------------------- */

    if (Common->status < KLU_OK)
    {
        LU = KLU_free (LU, lusize, sizeof (Unit), Common) ;
        lusize = 0 ;
    }
    *p_LU = LU ; //分配的内存地址
    return (lusize) ;
}

static void factor2
(
    /* inputs, not modified */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Entry Ax [ ],
    KLU_symbolic *Symbolic,

    /* inputs, modified on output: */
    KLU_numeric *Numeric,
    KLU_common *Common
)
{
    double lsize ;
    double *Lnz, *Rs ;
    Int *P, *Q, *R, *Pnum, *Offp, *Offi, *Pblock, *Pinv, *Iwork,
        *Lip, *Uip, *Llen, *Ulen ;
    Entry *Offx, *X, s, *Udiag ;
    Unit **LUbx ;
    Int k1, k2, nk, k, block, oldcol, pend, oldrow, n, lnz, unz, p, newrow,
        nblocks, poff, nzoff, lnz_block, unz_block, scale, max_lnz_block,
        max_unz_block ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* get the contents of the Symbolic object */
    n = Symbolic->n ;
    P = Symbolic->P ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    Lnz = Symbolic->Lnz ;
    nblocks = Symbolic->nblocks ;
    nzoff = Symbolic->nzoff ;

    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = (Entry *) Numeric->Offx ;

    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    Llen = Numeric->Llen ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;
    Udiag = Numeric->Udiag ;

    Rs = Numeric->Rs ;
    Pinv = Numeric->Pinv ;
    X = (Entry *) Numeric->Xwork ;              /* X is of size n */
    Iwork = Numeric->Iwork ;                    /* 5*maxblock for KLU_factor */
                                                /* 1*maxblock for Pblock */
    Pblock = Iwork + 5*((size_t) Symbolic->maxblock) ;
    Common->nrealloc = 0 ;
    scale = Common->scale ;
    max_lnz_block = 1 ;
    max_unz_block = 1 ;

  /*从符号分析中计算P的倒数。将更新为因式分解时*变成数值分解的逆已经完成，用于KLU_refactor */

    for (k = 0 ; k < n ; k++)
    {
        ASSERT (P [k] >= 0 && P [k] < n) ;
        Pinv [P [k]] = k ;
    }
    lnz = 0 ;
    unz = 0 ;
    Common->noffdiag = 0 ;
    Offp [0] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* optionally check input matrix and compute scale factors */
    /* ---------------------------------------------------------------------- */

    if (scale >= 0)
    {
        /*Pnum被用作工作空间，并且尺度因子在因子分解过程中的状态和变化。
        它强调了尺度因子在不同阶段与矩阵的行和最终的行置换之间的关系 */
        KLU_scale (scale, n, Ap, Ai, (double *) Ax, Rs, Pnum, Common) ;
        if (Common->status < KLU_OK)
        {
            /* matrix is invalid */
            return ;
        }
    }


    /* ---------------------------------------------------------------------- */
    /* factor each block using klu */
    /* ---------------------------------------------------------------------- */

    for (block = 0 ; block < nblocks ; block++) //遍历所有块
    {

        /* ------------------------------------------------------------------ */
        /* the block is from rows/columns k1 to k2-1 */
        /* ------------------------------------------------------------------ */

        k1 = R [block] ;
        k2 = R [block+1] ;
        nk = k2 - k1 ; //计算块大小

		//表示当前块是单个元素的块
        if (nk == 1)
        {

            /* -------------------------------------------------------------- */
            /* singleton case */
            /* -------------------------------------------------------------- */

            poff = Offp [k1] ;
            oldcol = Q [k1] ;
            pend = Ap [oldcol+1] ;
            CLEAR (s) ;

            if (scale <= 0) //此时没有进行缩放 则直接复制元素 
            {
                /* no scaling */
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    oldrow = Ai [p] ;
                    newrow = Pinv [oldrow] ;
                    if (newrow < k1)
                    {
                        Offi [poff] = oldrow ;
                        Offx [poff] = Ax [p] ;
                        poff++ ;
                    }
                    else
                    {
                        ASSERT (newrow == k1) ;
                        s = Ax [p] ;
                    }
                }
            }
            else //进行缩放
            {
         /* 正在进行行缩放操作。行缩放通常是为了改善数值稳定性或在特定算法中调整矩阵的性质。 */
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    oldrow = Ai [p] ;
                    newrow = Pinv [oldrow] ;
                    if (newrow < k1)
                    {
                        Offi [poff] = oldrow ;
                        /* Offx [poff] = Ax [p] / Rs [oldrow] ; */
                        SCALE_DIV_ASSIGN (Offx [poff], Ax [p], Rs [oldrow]) ;//除以对应缩放因子
                        poff++ ;
                    }
                    else
                    {
                        ASSERT (newrow == k1) ;
                        SCALE_DIV_ASSIGN (s, Ax [p], Rs [oldrow]) ;
                    }
                }
            }

            Udiag [k1] = s ; //指向U的对角线元素的指针

            if (IS_ZERO (s))
            {
                /* singular singleton */
                Common->status = KLU_SINGULAR ;
                Common->numerical_rank = k1 ;
                Common->singular_col = oldcol ;
                if (Common->halt_if_singular)
                {
                    return ;
                }
            }

            Offp [k1+1] = poff ;
            Pnum [k1] = P [k1] ;
            lnz++ ;
            unz++ ;

        }


        else //处理多个块
        {

            /* -------------------------------------------------------------- */
            /* construct and factorize the kth block */
            /* -------------------------------------------------------------- */

            if (Lnz [block] < 0)
            {
                /* COLAMD was used - no estimate of fill-in */
                /* use 10 times the nnz in A, plus n */
                lsize = -(Common->initmem) ;
            }
            else
            {
                lsize = Common->initmem_amd * Lnz [block] + nk ;
            }

            /* allocates 1 arrays: LUbx [block]  函数全部包含 LUsize被正确分配*/
			//对给定的矩阵进行因子分解
            Numeric->LUsize [block] = KLU_kernel_factor (nk, Ap, Ai, Ax, Q,
                    lsize, &LUbx [block], Udiag + k1, Llen + k1, Ulen + k1,
                    Lip + k1, Uip + k1, Pblock, &lnz_block, &unz_block,
                    X, Iwork, k1, Pinv, Rs, Offp, Offi, Offx, Common) ;

            if (Common->status < KLU_OK ||
               (Common->status == KLU_SINGULAR && Common->halt_if_singular))
            {
                /* out of memory, invalid inputs, or singular */
                return ;
            }

            /* -------------------------------------------------------------- */
            /* get statistics */
            /* -------------------------------------------------------------- */

            lnz += lnz_block ;
            unz += unz_block ;
            max_lnz_block = MAX (max_lnz_block, lnz_block) ;
            max_unz_block = MAX (max_unz_block, unz_block) ;

            if (Lnz [block] == EMPTY)
            {
                /* revise estimate for subsequent factorization */
                Lnz [block] = MAX (lnz_block, unz_block) ;
            }

            /* -------------------------------------------------------------- */
            /*结合klu行排序和符号预排序*/
            /* -------------------------------------------------------------- */
            for (k = 0 ; k < nk ; k++)
            {
                ASSERT (k + k1 < n) ;
                ASSERT (Pblock [k] + k1 < n) ;
                Pnum [k + k1] = P [Pblock [k] + k1] ;
            }

            /*不再需要本地主行排列Pblock */
        }
    }
    ASSERT (nzoff == Offp [n]) ;

    Numeric->lnz = lnz ;
    Numeric->unz = unz ;
    Numeric->max_lnz_block = max_lnz_block ;
    Numeric->max_unz_block = max_unz_block ;

    /* 计算Pnum的逆 */
    for (k = 0 ; k < n ; k++)
    {
        ASSERT (Pnum [k] >= 0 && Pnum [k] < n) ;
        Pinv [Pnum [k]] = k ;
    }

    /* 根据关键行顺序排列比例因子Rs */
    if (scale > 0)
    {
        for (k = 0 ; k < n ; k++)
        {
            REAL (X [k]) = Rs [Pnum [k]] ;
        }
        for (k = 0 ; k < n ; k++)
        {
            Rs [k] = REAL (X [k]) ;
        }
    }

    /* 将主行排列应用于非对角线项 */
    for (p = 0 ; p < nzoff ; p++)
    {
        ASSERT (Offi [p] >= 0 && Offi [p] < n) ;
        Offi [p] = Pinv [Offi [p]] ;
    }
}


KLU_numeric *KLU_factor         /* returns NULL if error, or a valid
                                   KLU_numeric object if successful */
(
    /* --- inputs --- */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    double Ax [ ],
    KLU_symbolic *Symbolic,
    /* -------------- */
    KLU_common *Common
)
{
    Int n, nzoff, nblocks, maxblock, k, ok = TRUE ;
    KLU_numeric *Numeric ;
    size_t n1, nzoff1, s, b6, n3 ;

    Common->status = KLU_OK ;
    Common->numerical_rank = EMPTY ;
    Common->singular_col = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    /* check for a valid Symbolic object */
    if (Symbolic == NULL)
    {
        Common->status = KLU_INVALID ;
        return (NULL) ;
    }

    n = Symbolic->n ;
    nzoff = Symbolic->nzoff ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;

    /* ---------------------------------------------------------------------- */
    /* get control parameters and make sure they are in the proper range */
    /* ---------------------------------------------------------------------- */

    Common->initmem_amd = MAX (1.0, Common->initmem_amd) ;
    Common->initmem = MAX (1.0, Common->initmem) ;
    Common->tol = MIN (Common->tol, 1.0) ;
    Common->tol = MAX (0.0, Common->tol) ;
    Common->memgrow = MAX (1.0, Common->memgrow) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the Numeric object  */
    /* ---------------------------------------------------------------------- */

    /* this will not cause size_t overflow (already checked by KLU_symbolic) */
    n1 = ((size_t) n) + 1 ;
    nzoff1 = ((size_t) nzoff) + 1 ;

    Numeric = KLU_malloc (1, sizeof (KLU_numeric), Common) ;
    if (Common->status < KLU_OK)
    {
        /* out of memory */
        Common->status = KLU_OUT_OF_MEMORY ;
        return (NULL) ;
    }
    Numeric->n = n ;
    Numeric->nblocks = nblocks ;
    Numeric->nzoff = nzoff ;
    Numeric->Pnum = KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Offp = KLU_malloc (n1, sizeof (Int), Common) ;
    Numeric->Offi = KLU_malloc (nzoff1, sizeof (Int), Common) ;
    Numeric->Offx = KLU_malloc (nzoff1, sizeof (Entry), Common) ;

    Numeric->Lip  = KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Uip  = KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Llen = KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Ulen = KLU_malloc (n, sizeof (Int), Common) ;

    Numeric->LUsize = KLU_malloc (nblocks, sizeof (size_t), Common) ;

    Numeric->LUbx = KLU_malloc (nblocks, sizeof (Unit *), Common) ;
    if (Numeric->LUbx != NULL)
    {
        for (k = 0 ; k < nblocks ; k++)
        {
            Numeric->LUbx [k] = NULL ;
        }
    }

    Numeric->Udiag = KLU_malloc (n, sizeof (Entry), Common) ;

    if (Common->scale > 0)
    {
        Numeric->Rs = KLU_malloc (n, sizeof (double), Common) ;
    }
    else
    {
        /* no scaling */
        Numeric->Rs = NULL ;
    }

    Numeric->Pinv = KLU_malloc (n, sizeof (Int), Common) ;

    /* allocate permanent workspace for factorization and solve.  Note that the
     * solver will use an Xwork of size 4n, whereas the factorization codes use
     * an Xwork of size n and integer space (Iwork) of size 6n. KLU_condest
     * uses an Xwork of size 2n.  Total size is:
     *
     *    n*sizeof(Entry) + max (6*maxblock*sizeof(Int), 3*n*sizeof(Entry))
     */
    s = KLU_mult_size_t (n, sizeof (Entry), &ok) ;
    n3 = KLU_mult_size_t (n, 3 * sizeof (Entry), &ok) ;
    b6 = KLU_mult_size_t (maxblock, 6 * sizeof (Int), &ok) ;
    Numeric->worksize = KLU_add_size_t (s, MAX (n3, b6), &ok) ;
    Numeric->Work = KLU_malloc (Numeric->worksize, 1, Common) ;
    Numeric->Xwork = Numeric->Work ;
    Numeric->Iwork = (Int *) ((Entry *) Numeric->Xwork + n) ;
    if (!ok || Common->status < KLU_OK)
    {
        /* out of memory or problem too large */
        Common->status = ok ? KLU_OUT_OF_MEMORY : KLU_TOO_LARGE ;
        KLU_free_numeric (&Numeric, Common) ;
        return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize the blocks */
    /* ---------------------------------------------------------------------- */

    factor2 (Ap, Ai, (Entry *) Ax, Symbolic, Numeric, Common) ;

    /* ---------------------------------------------------------------------- */
    /* return or free the Numeric object */
    /* ---------------------------------------------------------------------- */

    if (Common->status < KLU_OK)
    {
        /* out of memory or inputs invalid */
        KLU_free_numeric (&Numeric, Common) ;
    }
    else if (Common->status == KLU_SINGULAR)
    {
        if (Common->halt_if_singular)
        {
            /* Matrix is singular, and the Numeric object is only partially
             * defined because we halted early.  This is the default case for
             * a singular matrix. */
            KLU_free_numeric (&Numeric, Common) ;
        }
    }
    else if (Common->status == KLU_OK)
    {
        /* successful non-singular factorization */
        Common->numerical_rank = n ;
        Common->singular_col = n ;
    }
    return (Numeric) ;
}
