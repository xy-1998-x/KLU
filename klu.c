#include "klu.h"
#include "stdbool.h"
#include "sched.h"

//pid进程ID，policy调度策略，prio优先级
bool __rsn_sched_setscheduler(pid_t pid, int policy, int prio)
{
	struct sched_param param;
	memset(&param, 0x00, sizeof(param));
	param.sched_priority = prio;    //完成对调度参数结构体中优先级成员的设置
	
    //调用系统函数sched_setscheduler  尝试将指定进程（由pid指定）的调度策略设置为policy，并使用设置好的调度参数&param
	if (sched_setscheduler(pid, policy, &param) != 0)
	{
		return false;
	}
	return true;
}

//对上面函数的包装 将id直接设定为获取当前进程ID
bool rsn_sched_setscheduler(int policy, int prio)
{
    //这段代码使用getpid()这个函数来获取当前的进程ID
	return __rsn_sched_setscheduler(getpid(), policy, prio);
}

/* ========================================================================== */
/* === KLU_defaults ========================================================= */
/* ========================================================================== */
int KLU_defaults
(
    KLU_common *Common
)
{
    if (Common == NULL)
    {
        return (FALSE) ;
    }

    /* parameters */
    Common->tol = 0.001 ;       /* pivot tolerance for diagonal */
    Common->memgrow = 1.2;      /* realloc size ratio increase for LU factors */
    Common->initmem_amd = 1.2 ; /* init. mem with AMD:  c*nnz(L) + n */
    Common->initmem = 10 ;      /* init. mem otherwise: c*nnz(A) + n */
    Common->btf = TRUE ;        /* use BTF pre-ordering, or not */
    Common->maxwork = 0 ;       /* no limit to work done by btf_order */
    Common->ordering = 0 ;      /* 0: AMD, 1: COLAMD, 2: user-provided P and Q,
                                 * 3: user-provided function */
    Common->scale = 2 ;         /* scale: -1: none, and do not check for errors
                                 * in the input matrix in KLU_refactor.
                                 * 0: none, but check for errors,
                                 * 1: sum, 2: max */
    Common->halt_if_singular = TRUE ;   /* quick halt if matrix is singular */

    /* user ordering function and optional argument */
    Common->user_order = NULL ;
    Common->user_data = NULL ;

    /* statistics */
    Common->status = KLU_OK ;
    Common->nrealloc = 0 ;
    Common->structural_rank = EMPTY ;
    Common->numerical_rank = EMPTY ;
    Common->noffdiag = EMPTY ;
    Common->flops = EMPTY ;
    Common->rcond = EMPTY ;
    Common->condest = EMPTY ;
    Common->rgrowth = EMPTY ;
    Common->work = 0 ;          /* work done by btf_order */

    Common->memusage = 0 ;
    Common->mempeak = 0 ;

    return (TRUE) ;
}


/* ========================================================================== */
/* === KLU_analyze ========================================================== */
/* ========================================================================== */
static Int analyze_worker       /* returns KLU_OK or < 0 if error */
(
    /* inputs, not modified */
    Int n,              /* A is n-by-n */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Int nblocks,        /* # of blocks */
    Int Pbtf [ ],       /* BTF row permutation */
    Int Qbtf [ ],       /* BTF col permutation */
    Int R [ ],          /* size n+1, but only Rbtf [0..nblocks] is used */
    Int ordering,       /* what ordering to use (0, 1, or 3 for this routine) */

    /* output only, not defined on input */
    Int P [ ],          /* size n */
    Int Q [ ],          /* size n */
    double Lnz [ ],     /* size n, but only Lnz [0..nblocks-1] is used */

    /* workspace, not defined on input or output */
    Int Pblk [ ],       /* size maxblock */
    Int Cp [ ],         /* size maxblock+1 */
    Int Ci [ ],         /* size MAX (nz+1, Cilen) */
    Int Cilen,          /* nz+1, or COLAMD_recommend(nz,n,n) for COLAMD */
    Int Pinv [ ],       /* size maxblock */

    /* input/output */
    KLU_symbolic *Symbolic,
    KLU_common *Common
)
{
    double amd_Info [AMD_INFO], lnz, lnz1, flops, flops1 ;
    Int k1, k2, nk, k, block, oldcol, pend, newcol, result, pc, p, newrow,
        maxnz, nzoff, cstats [COLAMD_STATS], ok, err = KLU_INVALID ;


    for (k = 0 ; k < n ; k++)
    {
        ASSERT (Pbtf [k] >= 0 && Pbtf [k] < n) ;
        Pinv [Pbtf [k]] = k ;
    }
    nzoff = 0 ;
    lnz = 0 ;
    maxnz = 0 ;
    flops = 0 ;
    Symbolic->symmetry = EMPTY ;        /* only computed by AMD */

    /* ---------------------------------------------------------------------- */
    /* 排序每个块*/
    /* ---------------------------------------------------------------------- */

    for (block = 0 ; block < nblocks ; block++)
    {

        /* ------------------------------------------------------------------ */
        /* the block is from rows/columns k1 to k2-1 */
        /* ------------------------------------------------------------------ */

        k1 = R [block] ;
        k2 = R [block+1] ;
        nk = k2 - k1 ;

        /* ------------------------------------------------------------------ */
        /* construct the kth block, C */
        /* ------------------------------------------------------------------ */

        Lnz [block] = EMPTY ;
        pc = 0 ;
        for (k = k1 ; k < k2 ; k++)
        {
            newcol = k-k1 ;
            Cp [newcol] = pc ;
            oldcol = Qbtf [k] ;
            pend = Ap [oldcol+1] ;
            for (p = Ap [oldcol] ; p < pend ; p++)
            {
                newrow = Pinv [Ai [p]] ;
                if (newrow < k1)
                {
                    nzoff++ ;
                }
                else
                {
                    /* (newrow,newcol) is an entry in the block */
                    ASSERT (newrow < k2) ;
                    newrow -= k1 ;
                    Ci [pc++] = newrow ;
                }
            }
        }
        Cp [nk] = pc ;
        maxnz = MAX (maxnz, pc) ;

        /* ------------------------------------------------------------------ */
        /* order the block C */
        /* ------------------------------------------------------------------ */
 
        if (nk <= 3)
        {

            /* -------------------------------------------------------------- */
            /* use natural ordering for tiny blocks (3-by-3 or less) */
            /* -------------------------------------------------------------- */

            for (k = 0 ; k < nk ; k++)
            {
                Pblk [k] = k ;
            }
            lnz1 = nk * (nk + 1) / 2 ;
            flops1 = nk * (nk - 1) / 2 + (nk-1)*nk*(2*nk-1) / 6 ;
            ok = TRUE ;

        }
        else if (ordering == 0)
        {

            /* -------------------------------------------------------------- */
            /* order the block with AMD (C+C') */
            /* -------------------------------------------------------------- */

            result = AMD_order (nk, Cp, Ci, Pblk, NULL, amd_Info) ;
            ok = (result >= AMD_OK) ;
            if (result == AMD_OUT_OF_MEMORY)
            {
                err = KLU_OUT_OF_MEMORY ;
            }

            /* account for memory usage in AMD */
            Common->mempeak = MAX (Common->mempeak,
                Common->memusage + amd_Info [AMD_MEMORY]) ;

            /* get the ordering statistics from AMD */
            lnz1 = (Int) (amd_Info [AMD_LNZ]) + nk ;
            flops1 = 2 * amd_Info [AMD_NMULTSUBS_LU] + amd_Info [AMD_NDIV] ;
            if (pc == maxnz)
            {
                /* get the symmetry of the biggest block */
                Symbolic->symmetry = amd_Info [AMD_SYMMETRY] ;
            }

        }
        //否则使用自定义排序
        else
        {

            /* -------------------------------------------------------------- */
            /* pass the block to the user-provided ordering function */
            /* -------------------------------------------------------------- */

            lnz1 = (Common->user_order) (nk, Cp, Ci, Pblk, Common) ;
            flops1 = EMPTY ;
            ok = (lnz1 != 0) ;
        }

        if (!ok)
        {
            return (err) ;  /* ordering method failed */
        }

        /* ------------------------------------------------------------------ */
        /* keep track of nnz(L) and flops statistics */
        /* ------------------------------------------------------------------ */

        Lnz [block] = lnz1 ;
        lnz = (lnz == EMPTY || lnz1 == EMPTY) ? EMPTY : (lnz + lnz1) ;
        flops = (flops == EMPTY || flops1 == EMPTY) ? EMPTY : (flops + flops1) ;

        /* ------------------------------------------------------------------ */
        /* combine the preordering with the BTF ordering */
        /* ------------------------------------------------------------------ */

 
        for (k = 0 ; k < nk ; k++)
        {
            ASSERT (k + k1 < n) ;
            ASSERT (Pblk [k] + k1 < n) ;
            Q [k + k1] = Qbtf [Pblk [k] + k1] ;
        }
        for (k = 0 ; k < nk ; k++)
        {
            ASSERT (k + k1 < n) ;
            ASSERT (Pblk [k] + k1 < n) ;
            P [k + k1] = Pbtf [Pblk [k] + k1] ;
        }
    }

    
    ASSERT (nzoff >= 0 && nzoff <= Ap [n]) ;

    /* return estimates of # of nonzeros in L including diagonal */
    Symbolic->lnz = lnz ;           /* EMPTY if COLAMD used */
    Symbolic->unz = lnz ;
    Symbolic->nzoff = nzoff ;
    Symbolic->est_flops = flops ;   /* EMPTY if COLAMD or user-ordering used */
    return (KLU_OK) ;
}



static KLU_symbolic *order_and_analyze  /* returns NULL if error, or a valid
                                           KLU_symbolic object if successful */
(
    /* inputs, not modified */
    Int n,              /* A is n-by-n */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    /* --------------------- */
    KLU_common *Common
)
{
    double work ;
    KLU_symbolic *Symbolic ;
    double *Lnz ;
    Int *Qbtf, *Cp, *Ci, *Pinv, *Pblk, *Pbtf, *P, *Q, *R ;
    Int nblocks, nz, block, maxblock, k1, k2, nk, do_btf, ordering, k, Cilen,
        *Work ;

    /* ---------------------------------------------------------------------- */
    /* 检查输入矩阵并分配一个Symbolic对象 */
    /* ---------------------------------------------------------------------- */

    Symbolic = KLU_alloc_symbolic (n, Ap, Ai, Common) ;
    if (Symbolic == NULL)
    {
        return (NULL) ;
    }
    P = Symbolic->P ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    Lnz = Symbolic->Lnz ;
    nz = Symbolic->nz ;

    ordering = Common->ordering ;
	if (ordering == 0 || (ordering == 3 && Common->user_order != NULL))
    {
        /* AMD or user ordering function */
        Cilen = nz+1 ;
    }

    /* ---------------------------------------------------------------------- */
    /*  BTF 操作分配所需的内存空间 */
    /* ---------------------------------------------------------------------- */

    Pbtf = KLU_malloc (n, sizeof (Int), Common) ;
    Qbtf = KLU_malloc (n, sizeof (Int), Common) ;
 

    /* ---------------------------------------------------------------------- */
    /* 获取用于 BTF和排序方法的共同参数。 */
    /* ---------------------------------------------------------------------- */

    do_btf = Common->btf ; //Common->btf = 1
    do_btf = (do_btf) ? TRUE : FALSE ;
    Symbolic->ordering = ordering ;
    Symbolic->do_btf = do_btf ;
    Symbolic->structural_rank = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* 如果被请求，则找到矩阵的块三角形式。 */
    /* ---------------------------------------------------------------------- */

    Common->work = 0 ;

    if (do_btf)
    {
        Work = KLU_malloc (5*n, sizeof (Int), Common) ;
    
		//	BTF的函数依赖都正常导入
        nblocks = BTF_order (n, Ap, Ai, Common->maxwork, &work, Pbtf, Qbtf, R,
                &(Symbolic->structural_rank), Work) ;
        Common->structural_rank = Symbolic->structural_rank ;
        Common->work += work ;

        KLU_free (Work, 5*n, sizeof (Int), Common) ;

        /* unflip Qbtf if the matrix does not have full structural rank */
        if (Symbolic->structural_rank < n)
        {
            for (k = 0 ; k < n ; k++)
            {
                Qbtf [k] = BTF_UNFLIP (Qbtf [k]) ;
            }
        }

        /* find the size of the largest block */
        maxblock = 1 ;
        for (block = 0 ; block < nblocks ; block++)
        {
            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;
            maxblock = MAX (maxblock, nk) ;
        }
    }
    else
    {
        /* BTF not requested */
        nblocks = 1 ;
        maxblock = n ;
        R [0] = 0 ;
        R [1] = n ;
        for (k = 0 ; k < n ; k++)
        {
            Pbtf [k] = k ;
            Qbtf [k] = k ;
        }
    }

    Symbolic->nblocks = nblocks ;
    Symbolic->maxblock = maxblock ;

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, for analyze_worker */
    /* ---------------------------------------------------------------------- */

    Pblk = KLU_malloc (maxblock, sizeof (Int), Common) ;
    Cp   = KLU_malloc (maxblock + 1, sizeof (Int), Common) ;
    Ci   = KLU_malloc (MAX (Cilen, nz+1), sizeof (Int), Common) ;
    Pinv = KLU_malloc (n, sizeof (Int), Common) ;

    /* ---------------------------------------------------------------------- */
    /* 对块三角形式（BTF）排序中的每个块进行排序，并且进行一种减少填充的排序 */
    /* ---------------------------------------------------------------------- */

    if (Common->status == KLU_OK)
    {
        Common->status = analyze_worker (n, Ap, Ai, nblocks, Pbtf, Qbtf, R,
            ordering, P, Q, Lnz, Pblk, Cp, Ci, Cilen, Pinv, Symbolic, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* free all workspace */
    /* ---------------------------------------------------------------------- */

    KLU_free (Pblk, maxblock, sizeof (Int), Common) ;
    KLU_free (Cp, maxblock+1, sizeof (Int), Common) ;
    KLU_free (Ci, MAX (Cilen, nz+1), sizeof (Int), Common) ;
    KLU_free (Pinv, n, sizeof (Int), Common) ;
    KLU_free (Pbtf, n, sizeof (Int), Common) ;
    KLU_free (Qbtf, n, sizeof (Int), Common) ;

    /* ---------------------------------------------------------------------- */
    /* return the symbolic object */
    /* ---------------------------------------------------------------------- */

    return (Symbolic) ;
}

KLU_symbolic *KLU_analyze      
(
    /* inputs, not modified */
    Int n,              /* A is n-by-n */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    /* -------------------- */
    KLU_common *Common
)
{

    if (Common == NULL)
    {
        return (NULL) ;
    }
    Common->status = KLU_OK ;
    Common->structural_rank = EMPTY ;
    /* ---------------------------------------------------------------------- */
    /* order and analyze */
    /* ---------------------------------------------------------------------- */

        /* order with P and Q */
       
   return (  order_and_analyze (n, Ap, Ai, Common) );
}



Int BTF_order      /* returns number of blocks found */
(
    /* input, not modified: */
    Int n,          /* A is n-by-n in compressed column form */
    Int Ap [ ],     /* size n+1 */
    Int Ai [ ],     /* size nz = Ap [n] */
    double maxwork, /* do at most maxwork*nnz(A) work in the maximum
                     * transversal; no limit if <= 0 */

    /* output, not defined on input */
    double *work,   /* work performed in maxtrans, or -1 if limit reached */
    Int P [ ],      /* size n, row permutation */
    Int Q [ ],      /* size n, column permutation */
    Int R [ ],      /* size n+1.  block b is in rows/cols R[b] ... R[b+1]-1 */
    Int *nmatch,    /* # nonzeros on diagonal of P*A*Q */

    /* workspace, not defined on input or output */
    Int Work [ ]    /* size 5n */
)
{
    Int *Flag ;
    Int nblocks, i, j, nbadcol ;

    /* ---------------------------------------------------------------------- */
    /* compute the maximum matching */
    /* ---------------------------------------------------------------------- */

    /* if maxwork > 0, then a maximum matching might not be found */

    *nmatch = BTF_maxtrans (n, n, Ap, Ai, maxwork, work, Q, Work) ;

    /* ---------------------------------------------------------------------- */
    /* complete permutation if the matrix is structurally singular */
    /* ---------------------------------------------------------------------- */

    /* Since the matrix is square, ensure BTF_UNFLIP(Q[0..n-1]) is a
     * permutation of the columns of A so that A has as many nonzeros on the
     * diagonal as possible.
     */

    if (*nmatch < n)
    {
        /* get a size-n work array */
        Flag = Work + n ;
        for (j = 0 ; j < n ; j++)
        {
            Flag [j] = 0 ;
        }

        /* flag all matched columns */
        for (i = 0 ; i < n ; i++)
        {
            j = Q [i] ;
            if (j != EMPTY)
            {
                /* row i and column j are matched to each other */
                Flag [j] = 1 ;
            }
        }

        /* make a list of all unmatched columns, in Work [0..nbadcol-1]  */
        nbadcol = 0 ;
        for (j = n-1 ; j >= 0 ; j--)
        {
            if (!Flag [j])
            {
                /* j is matched to nobody */
                Work [nbadcol++] = j ;
            }
        }
        ASSERT (*nmatch + nbadcol == n) ;

        /* make an assignment for each unmatched row */
        for (i = 0 ; i < n ; i++)
        {
            if (Q [i] == EMPTY && nbadcol > 0)
            {
                /* get an unmatched column j */
                j = Work [--nbadcol] ;
                /* assign j to row i and flag the entry by "flipping" it */
                Q [i] = BTF_FLIP (j) ;
            }
        }
    }

    /* The permutation of a square matrix can be recovered as follows: Row i is
     * matched with column j, where j = BTF_UNFLIP (Q [i]) and where j
     * will always be in the valid range 0 to n-1.  The entry A(i,j) is zero
     * if BTF_ISFLIPPED (Q [i]) is true, and nonzero otherwise.  nmatch
     * is the number of entries in the Q array that are non-negative. */

    /* ---------------------------------------------------------------------- */
    /* find the strongly connected components */
    /* ---------------------------------------------------------------------- */

    nblocks = BTF_strongcomp (n, Ap, Ai, Q, P, R, Work) ;
    return (nblocks) ;
}


static int augment
(
    Int k,              /* which stage of the main loop we're in */
    Int Ap [ ],         /* column pointers, size n+1 */
    Int Ai [ ],         /* row indices, size nz = Ap [n] */
    Int Match [ ],      /* size n,  Match [i] = j if col j matched to i */
    Int Cheap [ ],      /* rows Ai [Ap [j] .. Cheap [j]-1] alread matched */
    Int Flag [ ],       /* Flag [j] = k if j already visited this stage */
    Int Istack [ ],     /* size n.  Row index stack. */
    Int Jstack [ ],     /* size n.  Column index stack. */
    Int Pstack [ ],     /* size n.  Keeps track of position in adjacency list */
    double *work,       /* work performed by the depth-first-search */
    double maxwork      /* maximum work allowed */
)
{
    /* local variables, but "global" to all DFS levels: */
    Int found ; /* true if match found.  */
    Int head ;  /* top of stack */

    /* variables that are purely local to any one DFS level: */
    Int j2 ;    /* the next DFS goes to node j2 */
    Int pend ;  /* one past the end of the adjacency list for node j */
    Int pstart ;
    Int quick ;

    /* variables that need to be pushed then popped from the stack: */
    Int i ;     /* the row tentatively matched to i if DFS successful */
    Int j ;     /* the DFS is at the current node j */
    Int p ;     /* current index into the adj. list for node j */
    /* the variables i, j, and p are stacked in Istack, Jstack, and Pstack */

    quick = (maxwork > 0) ;

    /* start a DFS to find a match for column k */
    found = FALSE ;
    i = EMPTY ;
    head = 0 ;
    Jstack [0] = k ;
    ASSERT (Flag [k] != k) ;
	
    while (head >= 0)
    {
        j = Jstack [head] ;
        pend = Ap [j+1] ;

        if (Flag [j] != k)          /* a node is not yet visited */
        {

            /* -------------------------------------------------------------- */
            /* prework for node j */
            /* -------------------------------------------------------------- */

            /* first time that j has been visited */
            Flag [j] = k ;
            /* cheap assignment: find the next unmatched row in col j.  This
             * loop takes at most O(nnz(A)) time for the sum total of all
             * calls to augment. */
            for (p = Cheap [j] ; p < pend && !found ; p++)
            {
                i = Ai [p] ;
                found = (Match [i] == EMPTY) ;
            }
            Cheap [j] = p ;

            /* -------------------------------------------------------------- */

            /* prepare for DFS */
            if (found)
            {
                /* end of augmenting path, column j matched with row i */
                Istack [head] = i ;
                break ;
            }
            /* set Pstack [head] to the first entry in column j to scan */
            Pstack [head] = Ap [j] ;
        }

        /* ------------------------------------------------------------------ */
        /* quick return if too much work done */
        /* ------------------------------------------------------------------ */

        if (quick && *work > maxwork)
        {
            /* too much work has been performed; abort the search */
            return (EMPTY) ;
        }

        /* ------------------------------------------------------------------ */
        /* DFS for nodes adjacent to j */
        /* ------------------------------------------------------------------ */

        /* If cheap assignment not made, continue the depth-first search.  All
         * rows in column j are already matched.  Add the adjacent nodes to the
         * stack by iterating through until finding another non-visited node.
         *
         * It is the following loop that can force maxtrans to take
         * O(n*nnz(A)) time. */

        pstart = Pstack [head] ;
        for (p = pstart ; p < pend ; p++)
        {
            i = Ai [p] ;
            j2 = Match [i] ;
            ASSERT (j2 != EMPTY) ;
            if (Flag [j2] != k)
            {
                /* Node j2 is not yet visited, start a depth-first search on
                 * node j2.  Keep track of where we left off in the scan of adj
                 * list of node j so we can restart j where we left off. */
                Pstack [head] = p + 1 ;
                /* Push j2 onto the stack and immediately break so we can
                 * recurse on node j2.  Also keep track of row i which (if this
                 * search for an augmenting path works) will be matched with the
                 * current node j. */
                Istack [head] = i ;
                Jstack [++head] = j2 ;
                break ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* determine how much work was just performed */
        /* ------------------------------------------------------------------ */

        *work += (p - pstart + 1) ;

        /* ------------------------------------------------------------------ */
        /* node j is done, but the postwork is postponed - see below */
        /* ------------------------------------------------------------------ */

        if (p == pend)
        {
            /* If all adjacent nodes of j are already visited, pop j from
             * stack and continue.  We failed to find a match. */
            head-- ;
        }
    }

    /* postwork for all nodes j in the stack */
    /* unwind the path and make the corresponding matches */
    if (found)
    {
        for (p = head ; p >= 0 ; p--)
        {
            j = Jstack [p] ;
            i = Istack [p] ;

            /* -------------------------------------------------------------- */
            /* postwork for node j */
            /* -------------------------------------------------------------- */
            /* if found, match row i with column j */
            Match [i] = j ;
        }
    }
    return (found) ;
}


Int BTF_maxtrans   /* returns # of columns in the matching */
(
    /* --- input --- */
    Int nrow,       /* A is nrow-by-ncol in compressed column form */
    Int ncol,
    Int Ap [ ],     /* size ncol+1 */
    Int Ai [ ],     /* size nz = Ap [ncol] */
    double maxwork, /* do at most maxwork*nnz(A) work; no limit if <= 0.  This
                     * work limit excludes the O(nnz(A)) cheap-match phase. */

    /* --- output --- */
    double *work,   /* work = -1 if maxwork > 0 and the total work performed
                     * reached the maximum of maxwork*nnz(A)).
                     * Otherwise, work = the total work performed. */

    Int Match [ ],  /* size nrow.  Match [i] = j if column j matched to row i */

    /* --- workspace --- */
    Int Work [ ]    /* size 5*ncol */
)
{
    Int *Cheap, *Flag, *Istack, *Jstack, *Pstack ;
    Int i, j, k, nmatch, work_limit_reached ;
    int result ;

    /* ---------------------------------------------------------------------- */
    /* get workspace and initialize */
    /* ---------------------------------------------------------------------- */

    Cheap  = Work ; Work += ncol ;
    Flag   = Work ; Work += ncol ;

    /* stack for non-recursive depth-first search in augment function */
    Istack = Work ; Work += ncol ;
    Jstack = Work ; Work += ncol ;
    Pstack = Work ;

    /* in column j, rows Ai [Ap [j] .. Cheap [j]-1] are known to be matched */
    for (j = 0 ; j < ncol ; j++)
    {
        Cheap [j] = Ap [j] ;
        Flag [j] = EMPTY ; 
    }

    /* all rows and columns are currently unmatched */
    for (i = 0 ; i < nrow ; i++)
    {
        Match [i] = EMPTY ;
    }

    if (maxwork > 0)
    {
        maxwork *= Ap [ncol] ;
    }
    *work = 0 ;

    /* ---------------------------------------------------------------------- */
    /* find a matching row for each column k */
    /* ---------------------------------------------------------------------- */

    nmatch = 0 ;
    work_limit_reached = FALSE ;
    for (k = 0 ; k < ncol ; k++)
    {
        /* find an augmenting path to match some row i to column k */
        result = augment (k, Ap, Ai, Match, Cheap, Flag, Istack, Jstack, Pstack,
            work, maxwork) ;
        if (result == TRUE)
        {
            /* we found it.  Match [i] = k for some row i has been done. */
            nmatch++ ;
        }
        else if (result == EMPTY)
        {
            /* augment gave up because of too much work, and no match found */
            work_limit_reached = TRUE ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* return the Match, and the # of matches made */
    /* ---------------------------------------------------------------------- */

    /* At this point, row i is matched to j = Match [i] if j >= 0.  i is an
     * unmatched row if Match [i] == EMPTY. */

    if (work_limit_reached)
    {
        /* return -1 if the work limit of maxwork*nnz(A) was reached */
        *work = EMPTY ;
    }

    return (nmatch) ;
}


static void dfs
(
    /* inputs, not modified on output: */
    Int j,              /* start the DFS at node j */
    Int Ap [ ],         /* size n+1, column pointers for the matrix A */
    Int Ai [ ],         /* row indices, size nz = Ap [n] */
    Int Q [ ],          /* input column permutation */

    /* inputs, modified on output (each array is of size n): */
    Int Time [ ],       /* Time [j] = "time" that node j was first visited */
    Int Flag [ ],       /* Flag [j]: see above */
    Int Low [ ],        /* Low [j]: see definition below */
    Int *p_nblocks,     /* number of blocks (aka strongly-connected-comp.)*/
    Int *p_timestamp,   /* current "time" */

    /* workspace, not defined on input or output: */
    Int Cstack [ ],     /* size n, output stack to hold nodes of components */
    Int Jstack [ ],     /* size n, stack for the variable j */
    Int Pstack [ ]      /* size n, stack for the variable p */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables, and initializations */
    /* ---------------------------------------------------------------------- */

    /* local variables, but "global" to all DFS levels: */
    Int chead ;     /* top of Cstack */
    Int jhead ;     /* top of Jstack and Pstack */

    /* variables that are purely local to any one DFS level: */
    Int i ;         /* edge (j,i) considered; i can be next node to traverse */
    Int parent ;    /* parent of node j in the DFS tree */
    Int pend ;      /* one past the end of the adjacency list for node j */
    Int jj ;        /* column j of A*Q is column jj of the input matrix A */

    /* variables that need to be pushed then popped from the stack: */
    Int p ;         /* current index into the adj. list for node j */
    /* the variables j and p are stacked in Jstack and Pstack */

    /* local copies of variables in the calling routine */
    Int nblocks   = *p_nblocks ;
    Int timestamp = *p_timestamp ;

    /* ---------------------------------------------------------------------- */
    /* start a DFS at node j (same as the recursive call dfs (EMPTY, j)) */
    /* ---------------------------------------------------------------------- */

    chead = 0 ;             /* component stack is empty */
    jhead = 0 ;             /* Jstack and Pstack are empty */
    Jstack [0] = j ;        /* put the first node j on the Jstack */
    ASSERT (Flag [j] == UNVISITED) ;

    while (jhead >= 0)
    {
        j = Jstack [jhead] ;        /* grab the node j from the top of Jstack */

        /* determine which column jj of the A is column j of A*Q */
        jj = (Q == (Int *) NULL) ? (j) : (BTF_UNFLIP (Q [j])) ;
        pend = Ap [jj+1] ;          /* j's row index list ends at Ai [pend-1] */

        if (Flag [j] == UNVISITED)
        {

            /* -------------------------------------------------------------- */
            /* prework at node j */
            /* -------------------------------------------------------------- */

            /* node j is being visited for the first time */
            Cstack [++chead] = j ;          /* push j onto the stack */
            timestamp++ ;                   /* get a timestamp */
            Time [j] = timestamp ;          /* give the timestamp to node j */
            Low [j] = timestamp ;
            Flag [j] = UNASSIGNED ;         /* flag node j as visited */

            /* -------------------------------------------------------------- */
            /* set Pstack [jhead] to the first entry in column j to scan */
            /* -------------------------------------------------------------- */

            Pstack [jhead] = Ap [jj] ;
        }

        /* ------------------------------------------------------------------ */
        /* DFS rooted at node j (start it, or continue where left off) */
        /* ------------------------------------------------------------------ */

        for (p = Pstack [jhead] ; p < pend ; p++)
        {
            i = Ai [p] ;    /* examine the edge from node j to node i */
            if (Flag [i] == UNVISITED)
            {
                /* Node i has not been visited - start a DFS at node i.
                 * Keep track of where we left off in the scan of adjacency list
                 * of node j so we can restart j where we left off. */
                Pstack [jhead] = p + 1 ;
                /* Push i onto the stack and immediately break
                 * so we can recurse on node i. */
                Jstack [++jhead] = i ;
                ASSERT (Time [i] == EMPTY) ;
                ASSERT (Low [i] == EMPTY) ;
                /* break here to do what the recursive call dfs (j,i) does */
                break ;
            }
            else if (Flag [i] == UNASSIGNED)
            {
                /* Node i has been visited, but still unassigned to a block
                 * this is a back or cross edge if Time [i] < Time [j].
                 * Note that i might equal j, in which case this code does
                 * nothing. */
                ASSERT (Time [i] > 0) ;
                ASSERT (Low [i] > 0) ;
                Low [j] = MIN (Low [j], Time [i]) ;
            }
        }

        if (p == pend)
        {
            /* If all adjacent nodes of j are already visited, pop j from
             * Jstack and do the post work for node j.  This also pops p
             * from the Pstack. */
            jhead-- ;

            /* -------------------------------------------------------------- */
            /* postwork at node j */
            /* -------------------------------------------------------------- */

            /* determine if node j is the head of a component */
            if (Low [j] == Time [j])
            {
                /* pop all nodes in this SCC from Cstack */
                while (TRUE)
                {
                    ASSERT (chead >= 0) ;       /* stack not empty (j in it) */
                    i = Cstack [chead--] ;      /* pop a node from the Cstack */
                    ASSERT (i >= 0) ;
                    ASSERT (Flag [i] == UNASSIGNED) ;
                    Flag [i] = nblocks ;        /* assign i to current block */
                    if (i == j) break ;         /* current block ends at j */
                }
                nblocks++ ;     /* one more block has been found */
            }
            /* update Low [parent], if the parent exists */
            if (jhead >= 0)
            {
                parent = Jstack [jhead] ;
                Low [parent] = MIN (Low [parent], Low [j]) ;
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* cleanup: update timestamp and nblocks */
    /* ---------------------------------------------------------------------- */

    *p_timestamp = timestamp ;
    *p_nblocks   = nblocks ;
}

Int BTF_strongcomp /* return # of strongly connected components */
(
    /* input, not modified: */
    Int n,          /* A is n-by-n in compressed column form */
    Int Ap [ ],     /* size n+1 */
    Int Ai [ ],     /* size nz = Ap [n] */

    /* optional input, modified (if present) on output: */
    Int Q [ ],      /* size n, input column permutation.  The permutation Q can
                     * include a flag which indicates an unmatched row.
                     * jold = BTF_UNFLIP (Q [jnew]) is the permutation;
                     * this function ingnores these flags.  On output, it is
                     * modified according to the permutation P. */

    /* output, not defined on input: */
    Int P [ ],      /* size n.  P [k] = j if row and column j are kth row/col
                     * in permuted matrix. */
    Int R [ ],      /* size n+1.  kth block is in rows/cols R[k] ... R[k+1]-1
                     * of the permuted matrix. */

    /* workspace, not defined on input or output: */
    Int Work [ ]    /* size 4n */
)
{
    Int j, k, b ;
    Int timestamp, nblocks, *Flag, *Cstack, *Time, *Low, *Jstack, *Pstack ; 

    Time   = Work ; Work += n ;
    Flag   = Work ; Work += n ;
    Low    = P ;                /* use output array P as workspace for Low */
    Cstack = R ;                /* use output array R as workspace for Cstack */

    /* stack for non-recursive dfs */
    Jstack = Work ; Work += n ;     /* stack for j */
    Pstack = Work ;                 /* stack for p */


    for (j = 0 ; j < n ; j++)
    {
        Flag [j] = UNVISITED ;
        Low [j] = EMPTY ;
        Time [j] = EMPTY ;

    }

    timestamp = 0 ;     /* each node given a timestamp when it is visited */
    nblocks = 0 ;       /* number of blocks found so far */



    for (j = 0 ; j < n ; j++)
    {
        /* node j is unvisited or assigned to a block. Cstack is empty. */
        ASSERT (Flag [j] == UNVISITED || (Flag [j] >= 0 && Flag [j] < nblocks));
        if (Flag [j] == UNVISITED)
        {
            /* non-recursive dfs (default) */
            dfs (j, Ap, Ai, Q, Time, Flag, Low, &nblocks, &timestamp,
                    Cstack, Jstack, Pstack) ;
        }
    }
    ASSERT (timestamp == n) ;

    /* ---------------------------------------------------------------------- */
    /* construct the block boundary array, R */
    /* ---------------------------------------------------------------------- */

    for (b = 0 ; b < nblocks ; b++)
    {
        R [b] = 0 ;
    }
    for (j = 0 ; j < n ; j++)
    {
        /* node j has been assigned to block b = Flag [j] */
        ASSERT (Time [j] > 0 && Time [j] <= n) ;
        ASSERT (Low [j] > 0 && Low [j] <= n) ;
        ASSERT (Flag [j] >= 0 && Flag [j] < nblocks) ;
        R [Flag [j]]++ ;
    }
    /* R [b] is now the number of nodes in block b.  Compute cumulative sum
     * of R, using Time [0 ... nblocks-1] as workspace. */
    Time [0] = 0 ;
    for (b = 1 ; b < nblocks ; b++)
    {
        Time [b] = Time [b-1] + R [b-1] ;
    }
    for (b = 0 ; b < nblocks ; b++)
    {
        R [b] = Time [b] ;
    }
    R [nblocks] = n ;

    for (j = 0 ; j < n ; j++)
    {
        /* place column j in the permutation */
        P [Time [Flag [j]]++] = j ;
    }


    if (Q != (Int *) NULL)
    {
        
        for (k = 0 ; k < n ; k++)
        {
            Time [k] = Q [P [k]] ;
        }
        for (k = 0 ; k < n ; k++)
        {
            Q [k] = Time [k] ;
        }
    }

    return (nblocks) ;
}

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

int KLU_free_symbolic
(
    KLU_symbolic **SymbolicHandle,
    KLU_common   *Common
)
{
    KLU_symbolic *Symbolic ;
    Int n ;
    if (Common == NULL)
    {
        return (FALSE) ;
    }
    if (SymbolicHandle == NULL || *SymbolicHandle == NULL)
    {
        return (TRUE) ;
    }
    Symbolic = *SymbolicHandle ;
    n = Symbolic->n ;
    KLU_free (Symbolic->P, n, sizeof (Int), Common) ;
    KLU_free (Symbolic->Q, n, sizeof (Int), Common) ;
    KLU_free (Symbolic->R, n+1, sizeof (Int), Common) ;
    KLU_free (Symbolic->Lnz, n, sizeof (double), Common) ;
    KLU_free (Symbolic, 1, sizeof (KLU_symbolic), Common) ;
    *SymbolicHandle = NULL ;
    return (TRUE) ;
}

int AMD_order
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    double Control [ ],
    double Info [ ]
)
{
    Int *Len, *S, nz, i, *Pinv, info, status, *Rp, *Ri, *Cp, *Ci, ok ;
    size_t nzaat, slen ;
    double mem = 0 ;

    /* clear the Info array, if it exists */
    info = Info != (double *) NULL ;
    if (info)
    {
	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [AMD_N] = n ;
	Info [AMD_STATUS] = AMD_OK ;
    }

    /* make sure inputs exist and n is >= 0 */
    if (Ai == (Int *) NULL || Ap == (Int *) NULL || P == (Int *) NULL || n < 0)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;	    /* arguments are invalid */
    }

    if (n == 0)
    {
	return (AMD_OK) ;	    /* n is 0 so there's nothing to do */
    }

    nz = Ap [n] ;
    if (info)
    {
	Info [AMD_NZ] = nz ;
    }
    if (nz < 0)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;
    }

    /* check if n or nz will cause integer overflow */
    if (((size_t) n) >= Int_MAX / sizeof (Int)
     || ((size_t) nz) >= Int_MAX / sizeof (Int))
    {
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;	    /* problem too large */
    }

    /* check the input matrix:	AMD_OK, AMD_INVALID, or AMD_OK_BUT_JUMBLED */
    status = AMD_valid (n, n, Ap, Ai) ; 

    if (status == AMD_INVALID)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;	    /* matrix is invalid */
    }

    /* allocate two size-n integer workspaces */
    size_t nn = (size_t) n ;
    Len  = SuiteSparse_malloc (nn, sizeof (Int)) ;
    Pinv = SuiteSparse_malloc (nn, sizeof (Int)) ;
    mem += n ;
    mem += n ;
    if (!Len || !Pinv)
    {
	/* :: out of memory :: */
	SuiteSparse_free (Len) ;
	SuiteSparse_free (Pinv) ;
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;
    }

    if (status == AMD_OK_BUT_JUMBLED)
    {
	/* sort the input matrix and remove duplicate entries */
	AMD_DEBUG1 (("Matrix is jumbled\n")) ;
	Rp = SuiteSparse_malloc (nn+1, sizeof (Int)) ;
	Ri = SuiteSparse_malloc (nz,  sizeof (Int)) ;
	mem += (n+1) ;
	mem += MAX (nz,1) ;
	if (!Rp || !Ri)
	{
	    /* :: out of memory :: */
	    SuiteSparse_free (Rp) ;
	    SuiteSparse_free (Ri) ;
	    SuiteSparse_free (Len) ;
	    SuiteSparse_free (Pinv) ;
	    if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	    return (AMD_OUT_OF_MEMORY) ;
	}
	/* use Len and Pinv as workspace to create R = A' */
	AMD_preprocess (n, Ap, Ai, Rp, Ri, Len, Pinv) ;
	Cp = Rp ;
	Ci = Ri ;
    }
    
    else
    {
	/* order the input matrix as-is.  No need to compute R = A' first */
	Rp = NULL ;
	Ri = NULL ;
	Cp = (Int *) Ap ;
	Ci = (Int *) Ai ;
    }

    /* --------------------------------------------------------------------- */
    /* determine the symmetry and count off-diagonal nonzeros in A+A' */
    /* --------------------------------------------------------------------- */

    nzaat = AMD_aat (n, Cp, Ci, Len, P, Info) ;
    AMD_DEBUG1 (("nzaat: %g\n", (double) nzaat)) ;
    ASSERT ((MAX (nz-n, 0) <= nzaat) && (nzaat <= 2 * (size_t) nz)) ;

    /* --------------------------------------------------------------------- */
    /* allocate workspace for matrix, elbow room, and 6 size-n vectors */
    /* --------------------------------------------------------------------- */

    S = NULL ;
    slen = nzaat ;			/* space for matrix */
    ok = ((slen + nzaat/5) >= slen) ;	/* check for size_t overflow */
    slen += nzaat/5 ;			/* add elbow room */
    for (i = 0 ; ok && i < 7 ; i++)
    {
	ok = ((slen + nn) > slen) ;	/* check for size_t overflow */
	slen += nn ;			/* size-n elbow room, 6 size-n work */
    }
    mem += slen ;
    ok = ok && (slen < SIZE_T_MAX / sizeof (Int)) ; /* check for overflow */
    if (ok)
    {
	S = SuiteSparse_malloc (slen, sizeof (Int)) ;
    }
    AMD_DEBUG1 (("slen %g\n", (double) slen)) ;
    if (!S)
    {
	/* :: out of memory :: (or problem too large) */
	SuiteSparse_free (Rp) ;
	SuiteSparse_free (Ri) ;
	SuiteSparse_free (Len) ;
	SuiteSparse_free (Pinv) ;
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;
    }
    if (info)
    {
	/* memory usage, in bytes. */
	Info [AMD_MEMORY] = mem * sizeof (Int) ;
    }

    /* --------------------------------------------------------------------- */
    /* order the matrix */
    /* --------------------------------------------------------------------- */

    AMD_1 (n, Cp, Ci, P, Pinv, Len, slen, S, Control, Info) ;

    /* --------------------------------------------------------------------- */
    /* free the workspace */
    /* --------------------------------------------------------------------- */

    SuiteSparse_free (Rp) ;
    SuiteSparse_free (Ri) ;
    SuiteSparse_free (Len) ;
    SuiteSparse_free (Pinv) ;
    SuiteSparse_free (S) ;
    if (info) Info [AMD_STATUS] = status ;
    return (status) ;	    /* successful ordering */
}


int AMD_valid
(
    /* inputs, not modified on output: */
    Int n_row,		/* A is n_row-by-n_col */
    Int n_col,
    const Int Ap [ ],	/* column pointers of A, of size n_col+1 */
    const Int Ai [ ]	/* row indices of A, of size nz = Ap [n_col] */
)
{
    Int nz, j, p1, p2, ilast, i, p ;
    int result = AMD_OK ;

    if (n_row < 0 || n_col < 0 || Ap == NULL || Ai == NULL)
    {
	return (AMD_INVALID) ;
    }
    nz = Ap [n_col] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	AMD_DEBUG0 (("column 0 pointer bad or nz < 0\n")) ;
	return (AMD_INVALID) ;
    }
    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	AMD_DEBUG2 (("\nColumn: "ID" p1: "ID" p2: "ID"\n", j, p1, p2)) ;
	if (p1 > p2)
	{
	    /* column pointers must be ascending */
	    AMD_DEBUG0 (("column "ID" pointer bad\n", j)) ;
	    return (AMD_INVALID) ;
	}
	ilast = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    AMD_DEBUG3 (("row: "ID"\n", i)) ;
	    if (i < 0 || i >= n_row)
	    {
		/* row index out of range */
		AMD_DEBUG0 (("index out of range, col "ID" row "ID"\n", j, i));
		return (AMD_INVALID) ;
	    }
	    if (i <= ilast)
	    {
		/* row index unsorted, or duplicate entry present */
		AMD_DEBUG1 (("index unsorted/dupl col "ID" row "ID"\n", j, i));
		result = AMD_OK_BUT_JUMBLED ;
	    }
	    ilast = i ;
	}
    }
    return (result) ;
}


size_t AMD_aat	/* returns nz in A+A' */
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],	/* Len [j]: length of column j of A+A', excl diagonal*/
    Int Tp [ ],		/* workspace of size n */
    double Info [ ]
)
{
    Int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz ;
    double sym ;
    size_t nzaat ;

    if (Info != (double *) NULL)
    {
	/* clear the Info array, if it exists */
	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [AMD_STATUS] = AMD_OK ;
    }

    for (k = 0 ; k < n ; k++)
    {
	Len [k] = 0 ;
    }

    nzdiag = 0 ;
    nzboth = 0 ;
    nz = Ap [n] ;

    for (k = 0 ; k < n ; k++)
    {
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;
	AMD_DEBUG2 (("\nAAT Column: "ID" p1: "ID" p2: "ID"\n", k, p1, p2)) ;

	/* construct A+A' */
	for (p = p1 ; p < p2 ; )
	{
	    /* scan the upper triangular part of A */
	    j = Ai [p] ;
	    if (j < k)
	    {
		/* entry A (j,k) is in the strictly upper triangular part,
		 * add both A (j,k) and A (k,j) to the matrix A+A' */
		Len [j]++ ;
		Len [k]++ ;
		AMD_DEBUG3 (("    upper ("ID","ID") ("ID","ID")\n", j,k, k,j));
		p++ ;
	    }
	    else if (j == k)
	    {
		/* skip the diagonal */
		p++ ;
		nzdiag++ ;
		break ;
	    }
	    else /* j > k */
	    {
		/* first entry below the diagonal */
		break ;
	    }
	    /* scan lower triangular part of A, in column j until reaching
	     * row k.  Start where last scan left off. */
	    ASSERT (Tp [j] != EMPTY) ;
	    ASSERT (Ap [j] <= Tp [j] && Tp [j] <= Ap [j+1]) ;
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		if (i < k)
		{
		    /* A (i,j) is only in the lower part, not in upper.
		     * add both A (i,j) and A (j,i) to the matrix A+A' */
		    Len [i]++ ;
		    Len [j]++ ;
		    AMD_DEBUG3 (("    lower ("ID","ID") ("ID","ID")\n",
			i,j, j,i)) ;
		    pj++ ;
		}
		else if (i == k)
		{
		    /* entry A (k,j) in lower part and A (j,k) in upper */
		    pj++ ;
		    nzboth++ ;
		    break ;
		}
		else /* i > k */
		{
		    /* consider this entry later, when k advances to i */
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	/* Tp [k] points to the entry just below the diagonal in column k */
	Tp [k] = p ;
    }

    /* clean up, for remaining mismatched entries */
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    /* A (i,j) is only in the lower part, not in upper.
	     * add both A (i,j) and A (j,i) to the matrix A+A' */
	    Len [i]++ ;
	    Len [j]++ ;
	    AMD_DEBUG3 (("    lower cleanup ("ID","ID") ("ID","ID")\n",
		i,j, j,i)) ;
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the symmetry of the nonzero pattern of A */
    /* --------------------------------------------------------------------- */

    /* Given a matrix A, the symmetry of A is:
     *	B = tril (spones (A), -1) + triu (spones (A), 1) ;
     *  sym = nnz (B & B') / nnz (B) ;
     *  or 1 if nnz (B) is zero.
     */

    if (nz == nzdiag)
    {
	sym = 1 ;
    }
    else
    {
	sym = (2 * (double) nzboth) / ((double) (nz - nzdiag)) ;
    }

    nzaat = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	nzaat += Len [k] ;
    }

    AMD_DEBUG1 (("AMD nz in A+A', excluding diagonal (nzaat) = %g\n",
	(double) nzaat)) ;
    AMD_DEBUG1 (("   nzboth: "ID" nz: "ID" nzdiag: "ID" symmetry: %g\n",
		nzboth, nz, nzdiag, sym)) ;

    if (Info != (double *) NULL)
    {
	Info [AMD_STATUS] = AMD_OK ;
	Info [AMD_N] = n ;
	Info [AMD_NZ] = nz ;
	Info [AMD_SYMMETRY] = sym ;	    /* symmetry of pattern of A */
	Info [AMD_NZDIAG] = nzdiag ;	    /* nonzeros on diagonal of A */
	Info [AMD_NZ_A_PLUS_AT] = nzaat ;   /* nonzeros in A+A' */
    }

    return (nzaat) ;
}


void AMD_1
(
    Int n,		/* n > 0 */
    const Int Ap [ ],	/* input of size n+1, not modified */
    const Int Ai [ ],	/* input of size nz = Ap [n], not modified */
    Int P [ ],		/* size n output permutation */
    Int Pinv [ ],	/* size n output inverse permutation */
    Int Len [ ],	/* size n input, undefined on output */
    Int slen,		/* slen >= sum (Len [0..n-1]) + 7n,
			 * ideally slen = 1.2 * sum (Len) + 8n */
    Int S [ ],		/* size slen workspace */
    double Control [ ],	/* input array of size AMD_CONTROL */
    double Info [ ]	/* output array of size AMD_INFO */
)
{
    Int i, j, k, p, pfree, iwlen, pj, p1, p2, pj2, *Iw, *Pe, *Nv, *Head,
	*Elen, *Degree, *s, *W, *Sp, *Tp ;

    /* --------------------------------------------------------------------- */
    /* construct the matrix for AMD_2 */
    /* --------------------------------------------------------------------- */

    ASSERT (n > 0) ;

    iwlen = slen - 6*n ;
    s = S ;
    Pe = s ;	    s += n ;
    Nv = s ;	    s += n ;
    Head = s ;	    s += n ;
    Elen = s ;	    s += n ;
    Degree = s ;    s += n ;
    W = s ;	    s += n ;
    Iw = s ;	    s += iwlen ;

    ASSERT (AMD_valid (n, n, Ap, Ai) == AMD_OK) ;

    /* construct the pointers for A+A' */
    Sp = Nv ;			/* use Nv and W as workspace for Sp and Tp [ */
    Tp = W ;
    pfree = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	Pe [j] = pfree ;
	Sp [j] = pfree ;
	pfree += Len [j] ;
    }

    /* Note that this restriction on iwlen is slightly more restrictive than
     * what is strictly required in AMD_2.  AMD_2 can operate with no elbow
     * room at all, but it will be very slow.  For better performance, at
     * least size-n elbow room is enforced. */
    ASSERT (iwlen >= pfree + n) ;


    for (k = 0 ; k < n ; k++)
    {
	AMD_DEBUG1 (("Construct row/column k= "ID" of A+A'\n", k))  ;
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;

	/* construct A+A' */
	for (p = p1 ; p < p2 ; )
	{
	    /* scan the upper triangular part of A */
	    j = Ai [p] ;
	    ASSERT (j >= 0 && j < n) ;
	    if (j < k)
	    {
		/* entry A (j,k) in the strictly upper triangular part */
		ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
		ASSERT (Sp [k] < (k == n-1 ? pfree : Pe [k+1])) ;
		Iw [Sp [j]++] = k ;
		Iw [Sp [k]++] = j ;
		p++ ;
	    }
	    else if (j == k)
	    {
		/* skip the diagonal */
		p++ ;
		break ;
	    }
	    else /* j > k */
	    {
		/* first entry below the diagonal */
		break ;
	    }
	    /* scan lower triangular part of A, in column j until reaching
	     * row k.  Start where last scan left off. */
	    ASSERT (Ap [j] <= Tp [j] && Tp [j] <= Ap [j+1]) ;
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		ASSERT (i >= 0 && i < n) ;
		if (i < k)
		{
		    /* A (i,j) is only in the lower part, not in upper */
		    ASSERT (Sp [i] < (i == n-1 ? pfree : Pe [i+1])) ;
		    ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
		    Iw [Sp [i]++] = j ;
		    Iw [Sp [j]++] = i ;
		    pj++ ;
		}
		else if (i == k)
		{
		    /* entry A (k,j) in lower part and A (j,k) in upper */
		    pj++ ;
		    break ;
		}
		else /* i > k */
		{
		    /* consider this entry later, when k advances to i */
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	Tp [k] = p ;
    }

    /* clean up, for remaining mismatched entries */
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    ASSERT (i >= 0 && i < n) ;
	    /* A (i,j) is only in the lower part, not in upper */
	    ASSERT (Sp [i] < (i == n-1 ? pfree : Pe [i+1])) ;
	    ASSERT (Sp [j] < (j == n-1 ? pfree : Pe [j+1])) ;
	    Iw [Sp [i]++] = j ;
	    Iw [Sp [j]++] = i ;
	}
    }

    /* Tp and Sp no longer needed ] */

    /* --------------------------------------------------------------------- */
    /* order the matrix */
    /* --------------------------------------------------------------------- */

    AMD_2 (n, Pe, Iw, Len, iwlen, pfree,
	Nv, Pinv, P, Head, Elen, Degree, W, Control, Info) ;
}


static Int clear_flag (Int wflg, Int wbig, Int W [ ], Int n)
{
    Int x ;
    if (wflg < 2 || wflg >= wbig)
    {
	for (x = 0 ; x < n ; x++)
	{
	    if (W [x] != 0) W [x] = 1 ;
	}
	wflg = 2 ;
    }
    /*  at this point, W [0..n-1] < wflg holds */
    return (wflg) ;
}


/* ========================================================================= */
/* === AMD_2 =============================================================== */
/* ========================================================================= */

void AMD_2
(
    Int n,		/* A is n-by-n, where n > 0 */
    Int Pe [ ],		/* Pe [0..n-1]: index in Iw of row i on input */
    Int Iw [ ],		/* workspace of size iwlen. Iw [0..pfree-1]
			 * holds the matrix on input */
    Int Len [ ],	/* Len [0..n-1]: length for row/column i on input */
    Int iwlen,		/* length of Iw. iwlen >= pfree + n */
    Int pfree,		/* Iw [pfree ... iwlen-1] is empty on input */

    /* 7 size-n workspaces, not defined on input: */
    Int Nv [ ],		/* the size of each supernode on output */
    Int Next [ ],	/* the output inverse permutation */
    Int Last [ ],	/* the output permutation */
    Int Head [ ],
    Int Elen [ ],	/* the size columns of L for each supernode */
    Int Degree [ ],
    Int W [ ],

    /* control parameters and output statistics */
    double Control [ ],	/* array of size AMD_CONTROL */
    double Info [ ]	/* array of size AMD_INFO */
)
{
    Int deg, degme, dext, lemax, e, elenme, eln, i, ilast, inext, j,
	jlast, jnext, k, knt1, knt2, knt3, lenj, ln, me, mindeg, nel, nleft,
	nvi, nvj, nvpiv, slenme, wbig, we, wflg, wnvi, ok, ndense, ncmpa,
	dense, aggressive ;

    UInt hash ;	    /* unsigned, so that hash % n is well defined.*/


    double f, r, ndiv, s, nms_lu, nms_ldl, dmax, alpha, lnz, lnzme ;


    Int p, p1, p2, p3, p4, pdst, pend, pj, pme, pme1, pme2, pn, psrc ;

    ASSERT (iwlen >= pfree + n) ;
    ASSERT (n > 0) ;

    /* initialize output statistics */
    lnz = 0 ;
    ndiv = 0 ;
    nms_lu = 0 ;
    nms_ldl = 0 ;
    dmax = 1 ;
    me = EMPTY ;

    mindeg = 0 ;
    ncmpa = 0 ;
    nel = 0 ;
    lemax = 0 ;

    /* get control parameters */
    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = (Control [AMD_AGGRESSIVE] != 0) ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }
    /* Note: if alpha is NaN, this is undefined: */
    if (alpha < 0)
    {
	/* only remove completely dense rows/columns */
	dense = n-2 ;
    }
    else
    {
	dense = alpha * sqrt ((double) n) ;      //sqrt 开方函数 undefine
    }
    dense = MAX (16, dense) ;
    dense = MIN (n,  dense) ;
    AMD_DEBUG1 (("\n\nAMD (debug), alpha %g, aggr. "ID"\n",
	alpha, aggressive)) ;

    for (i = 0 ; i < n ; i++)
    {
	Last [i] = EMPTY ;
	Head [i] = EMPTY ;
	Next [i] = EMPTY ;
	/* if separate Hhead array is used for hash buckets: *
	Hhead [i] = EMPTY ;
	*/
	Nv [i] = 1 ;
	W [i] = 1 ;
	Elen [i] = 0 ;
	Degree [i] = Len [i] ;
    }

    /* initialize wflg */
    wbig = Int_MAX - n ;
    wflg = clear_flag (0, wbig, W, n) ;

    /* --------------------------------------------------------------------- */
    /* initialize degree lists and eliminate dense and empty rows */
    /* --------------------------------------------------------------------- */

    ndense = 0 ;

    for (i = 0 ; i < n ; i++)
    {
	deg = Degree [i] ;
	ASSERT (deg >= 0 && deg < n) ;
	if (deg == 0)
	{

	    /* -------------------------------------------------------------
	     * we have a variable that can be eliminated at once because
	     * there is no off-diagonal non-zero in its row.  Note that
	     * Nv [i] = 1 for an empty variable i.  It is treated just
	     * the same as an eliminated element i.
	     * ------------------------------------------------------------- */

	    Elen [i] = FLIP (1) ;
	    nel++ ;
	    Pe [i] = EMPTY ;
	    W [i] = 0 ;

	}
	else if (deg > dense)
	{

	    /* -------------------------------------------------------------
	     * Dense variables are not treated as elements, but as unordered,
	     * non-principal variables that have no parent.  They do not take
	     * part in the postorder, since Nv [i] = 0.  Note that the Fortran
	     * version does not have this option.
	     * ------------------------------------------------------------- */

	    AMD_DEBUG1 (("Dense node "ID" degree "ID"\n", i, deg)) ;
	    ndense++ ;
	    Nv [i] = 0 ;		/* do not postorder this node */
	    Elen [i] = EMPTY ;
	    nel++ ;
	    Pe [i] = EMPTY ;

	}
	else
	{

	    /* -------------------------------------------------------------
	     * place i in the degree list corresponding to its degree
	     * ------------------------------------------------------------- */

	    inext = Head [deg] ;
	    ASSERT (inext >= EMPTY && inext < n) ;
	    if (inext != EMPTY) Last [inext] = i ;
	    Next [i] = inext ;
	    Head [deg] = i ;

	}
    }

/* ========================================================================= */
/* WHILE (selecting pivots) DO */
/* ========================================================================= */

    while (nel < n)
    {

/* ========================================================================= */
/* GET PIVOT OF MINIMUM DEGREE */
/* ========================================================================= */

	/* ----------------------------------------------------------------- */
	/* find next supervariable for elimination */
	/* ----------------------------------------------------------------- */

	ASSERT (mindeg >= 0 && mindeg < n) ;
	for (deg = mindeg ; deg < n ; deg++)
	{
	    me = Head [deg] ;
	    if (me != EMPTY) break ;
	}
	mindeg = deg ;
	ASSERT (me >= 0 && me < n) ;
	AMD_DEBUG1 (("=================me: "ID"\n", me)) ;

	/* ----------------------------------------------------------------- */
	/* remove chosen variable from link list */
	/* ----------------------------------------------------------------- */

	inext = Next [me] ;
	ASSERT (inext >= EMPTY && inext < n) ;
	if (inext != EMPTY) Last [inext] = EMPTY ;
	Head [deg] = inext ;

	/* ----------------------------------------------------------------- */
	/* me represents the elimination of pivots nel to nel+Nv[me]-1. */
	/* place me itself as the first in this set. */
	/* ----------------------------------------------------------------- */

	elenme = Elen [me] ;
	nvpiv = Nv [me] ;
	ASSERT (nvpiv > 0) ;
	nel += nvpiv ;

/* ========================================================================= */
/* CONSTRUCT NEW ELEMENT */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * At this point, me is the pivotal supervariable.  It will be
	 * converted into the current element.  Scan list of the pivotal
	 * supervariable, me, setting tree pointers and constructing new list
	 * of supervariables for the new element, me.  p is a pointer to the
	 * current position in the old list.
	 * ----------------------------------------------------------------- */

	/* flag the variable "me" as being in Lme by negating Nv [me] */
	Nv [me] = -nvpiv ;
	degme = 0 ;
	ASSERT (Pe [me] >= 0 && Pe [me] < iwlen) ;

	if (elenme == 0)
	{

	    /* ------------------------------------------------------------- */
	    /* construct the new element in place */
	    /* ------------------------------------------------------------- */

	    pme1 = Pe [me] ;
	    pme2 = pme1 - 1 ;

	    for (p = pme1 ; p <= pme1 + Len [me] - 1 ; p++)
	    {
		i = Iw [p] ;
		ASSERT (i >= 0 && i < n && Nv [i] >= 0) ;
		nvi = Nv [i] ;
		if (nvi > 0)
		{

		    /* ----------------------------------------------------- */
		    /* i is a principal variable not yet placed in Lme. */
		    /* store i in new list */
		    /* ----------------------------------------------------- */

		    /* flag i as being in Lme by negating Nv [i] */
		    degme += nvi ;
		    Nv [i] = -nvi ;
		    Iw [++pme2] = i ;

		    /* ----------------------------------------------------- */
		    /* remove variable i from degree list. */
		    /* ----------------------------------------------------- */

		    ilast = Last [i] ;
		    inext = Next [i] ;
		    ASSERT (ilast >= EMPTY && ilast < n) ;
		    ASSERT (inext >= EMPTY && inext < n) ;
		    if (inext != EMPTY) Last [inext] = ilast ;
		    if (ilast != EMPTY)
		    {
			Next [ilast] = inext ;
		    }
		    else
		    {
			/* i is at the head of the degree list */
			ASSERT (Degree [i] >= 0 && Degree [i] < n) ;
			Head [Degree [i]] = inext ;
		    }
		}
	    }
	}
	else
	{

	    /* ------------------------------------------------------------- */
	    /* construct the new element in empty space, Iw [pfree ...] */
	    /* ------------------------------------------------------------- */

	    p = Pe [me] ;
	    pme1 = pfree ;
	    slenme = Len [me] - elenme ;

	    for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)
	    {

		if (knt1 > elenme)
		{
		    /* search the supervariables in me. */
		    e = me ;
		    pj = p ;
		    ln = slenme ;
		    AMD_DEBUG2 (("Search sv: "ID" "ID" "ID"\n", me,pj,ln)) ;
		}
		else
		{
		    /* search the elements in me. */
		    e = Iw [p++] ;
		    ASSERT (e >= 0 && e < n) ;
		    pj = Pe [e] ;
		    ln = Len [e] ;
		    AMD_DEBUG2 (("Search element e "ID" in me "ID"\n", e,me)) ;
		    ASSERT (Elen [e] < EMPTY && W [e] > 0 && pj >= 0) ;
		}
		ASSERT (ln >= 0 && (ln == 0 || (pj >= 0 && pj < iwlen))) ;

		/* ---------------------------------------------------------
		 * search for different supervariables and add them to the
		 * new list, compressing when necessary. this loop is
		 * executed once for each element in the list and once for
		 * all the supervariables in the list.
		 * --------------------------------------------------------- */

		for (knt2 = 1 ; knt2 <= ln ; knt2++)
		{
		    i = Iw [pj++] ;
		    ASSERT (i >= 0 && i < n && (i == me || Elen [i] >= EMPTY));
		    nvi = Nv [i] ;
		    AMD_DEBUG2 ((": "ID" "ID" "ID" "ID"\n",
				i, Elen [i], Nv [i], wflg)) ;

		    if (nvi > 0)
		    {

			/* ------------------------------------------------- */
			/* compress Iw, if necessary */
			/* ------------------------------------------------- */

			if (pfree >= iwlen)
			{

			    AMD_DEBUG1 (("GARBAGE COLLECTION\n")) ;

			    /* prepare for compressing Iw by adjusting pointers
			     * and lengths so that the lists being searched in
			     * the inner and outer loops contain only the
			     * remaining entries. */

			    Pe [me] = p ;
			    Len [me] -= knt1 ;
			    /* check if nothing left of supervariable me */
			    if (Len [me] == 0) Pe [me] = EMPTY ;
			    Pe [e] = pj ;
			    Len [e] = ln - knt2 ;
			    /* nothing left of element e */
			    if (Len [e] == 0) Pe [e] = EMPTY ;

			    ncmpa++ ;	/* one more garbage collection */

			    /* store first entry of each object in Pe */
			    /* FLIP the first entry in each object */
			    for (j = 0 ; j < n ; j++)
			    {
				pn = Pe [j] ;
				if (pn >= 0)
				{
				    ASSERT (pn >= 0 && pn < iwlen) ;
				    Pe [j] = Iw [pn] ;
				    Iw [pn] = FLIP (j) ;
				}
			    }

			    /* psrc/pdst point to source/destination */
			    psrc = 0 ;
			    pdst = 0 ;
			    pend = pme1 - 1 ;

			    while (psrc <= pend)
			    {
				/* search for next FLIP'd entry */
				j = FLIP (Iw [psrc++]) ;
				if (j >= 0)
				{
				    AMD_DEBUG2 (("Got object j: "ID"\n", j)) ;
				    Iw [pdst] = Pe [j] ;
				    Pe [j] = pdst++ ;
				    lenj = Len [j] ;
				    /* copy from source to destination */
				    for (knt3 = 0 ; knt3 <= lenj - 2 ; knt3++)
				    {
					Iw [pdst++] = Iw [psrc++] ;
				    }
				}
			    }

			    /* move the new partially-constructed element */
			    p1 = pdst ;
			    for (psrc = pme1 ; psrc <= pfree-1 ; psrc++)
			    {
				Iw [pdst++] = Iw [psrc] ;
			    }
			    pme1 = p1 ;
			    pfree = pdst ;
			    pj = Pe [e] ;
			    p = Pe [me] ;

			}

			/* ------------------------------------------------- */
			/* i is a principal variable not yet placed in Lme */
			/* store i in new list */
			/* ------------------------------------------------- */

			/* flag i as being in Lme by negating Nv [i] */
			degme += nvi ;
			Nv [i] = -nvi ;
			Iw [pfree++] = i ;
			AMD_DEBUG2 (("     s: "ID"     nv "ID"\n", i, Nv [i]));

			/* ------------------------------------------------- */
			/* remove variable i from degree link list */
			/* ------------------------------------------------- */

			ilast = Last [i] ;
			inext = Next [i] ;
			ASSERT (ilast >= EMPTY && ilast < n) ;
			ASSERT (inext >= EMPTY && inext < n) ;
			if (inext != EMPTY) Last [inext] = ilast ;
			if (ilast != EMPTY)
			{
			    Next [ilast] = inext ;
			}
			else
			{
			    /* i is at the head of the degree list */
			    ASSERT (Degree [i] >= 0 && Degree [i] < n) ;
			    Head [Degree [i]] = inext ;
			}
		    }
		}

		if (e != me)
		{
		    /* set tree pointer and flag to indicate element e is
		     * absorbed into new element me (the parent of e is me) */
		    AMD_DEBUG1 ((" Element "ID" => "ID"\n", e, me)) ;
		    Pe [e] = FLIP (me) ;
		    W [e] = 0 ;
		}
	    }

	    pme2 = pfree - 1 ;
	}

	/* ----------------------------------------------------------------- */
	/* me has now been converted into an element in Iw [pme1..pme2] */
	/* ----------------------------------------------------------------- */

	/* degme holds the external degree of new element */
	Degree [me] = degme ;
	Pe [me] = pme1 ;
	Len [me] = pme2 - pme1 + 1 ;
	ASSERT (Pe [me] >= 0 && Pe [me] < iwlen) ;

	Elen [me] = FLIP (nvpiv + degme) ;
	/* FLIP (Elen (me)) is now the degree of pivot (including
	 * diagonal part). */


	/* ----------------------------------------------------------------- */
	/* make sure that wflg is not too large. */
	/* ----------------------------------------------------------------- */

	/* With the current value of wflg, wflg+n must not cause integer
	 * overflow */

	wflg = clear_flag (wflg, wbig, W, n) ;

/* ========================================================================= */
/* COMPUTE (W [e] - wflg) = |Le\Lme| FOR ALL ELEMENTS */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * Scan 1:  compute the external degrees of previous elements with
	 * respect to the current element.  That is:
	 *       (W [e] - wflg) = |Le \ Lme|
	 * for each element e that appears in any supervariable in Lme.  The
	 * notation Le refers to the pattern (list of supervariables) of a
	 * previous element e, where e is not yet absorbed, stored in
	 * Iw [Pe [e] + 1 ... Pe [e] + Len [e]].  The notation Lme
	 * refers to the pattern of the current element (stored in
	 * Iw [pme1..pme2]).   If aggressive absorption is enabled, and
	 * (W [e] - wflg) becomes zero, then the element e will be absorbed
	 * in Scan 2.
	 * ----------------------------------------------------------------- */

	AMD_DEBUG2 (("me: ")) ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    eln = Elen [i] ;
	    AMD_DEBUG3 ((""ID" Elen "ID": \n", i, eln)) ;
	    if (eln > 0)
	    {
		/* note that Nv [i] has been negated to denote i in Lme: */
		nvi = -Nv [i] ;
		ASSERT (nvi > 0 && Pe [i] >= 0 && Pe [i] < iwlen) ;
		wnvi = wflg - nvi ;
		for (p = Pe [i] ; p <= Pe [i] + eln - 1 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    AMD_DEBUG4 (("    e "ID" we "ID" ", e, we)) ;
		    if (we >= wflg)
		    {
			/* unabsorbed element e has been seen in this loop */
			AMD_DEBUG4 (("    unabsorbed, first time seen")) ;
			we -= nvi ;
		    }
		    else if (we != 0)
		    {
			/* e is an unabsorbed element */
			/* this is the first we have seen e in all of Scan 1 */
			AMD_DEBUG4 (("    unabsorbed")) ;
			we = Degree [e] + wnvi ;
		    }
		    AMD_DEBUG4 (("\n")) ;
		    W [e] = we ;
		}
	    }
	}
	AMD_DEBUG2 (("\n")) ;

/* ========================================================================= */
/* DEGREE UPDATE AND ELEMENT ABSORPTION */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * Scan 2:  for each i in Lme, sum up the degree of Lme (which is
	 * degme), plus the sum of the external degrees of each Le for the
	 * elements e appearing within i, plus the supervariables in i.
	 * Place i in hash list.
	 * ----------------------------------------------------------------- */

	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n && Nv [i] < 0 && Elen [i] >= 0) ;
	    AMD_DEBUG2 (("Updating: i "ID" "ID" "ID"\n", i, Elen[i], Len [i]));
	    p1 = Pe [i] ;
	    p2 = p1 + Elen [i] - 1 ;
	    pn = p1 ;
	    hash = 0 ;
	    deg = 0 ;
	    ASSERT (p1 >= 0 && p1 < iwlen && p2 >= -1 && p2 < iwlen) ;

	    /* ------------------------------------------------------------- */
	    /* scan the element list associated with supervariable i */
	    /* ------------------------------------------------------------- */

	    /* UMFPACK/MA38-style approximate degree: */
	    if (aggressive)
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			/* dext = | Le \ Lme | */
			dext = we - wflg ;
			if (dext > 0)
			{
			    deg += dext ;
			    Iw [pn++] = e ;
			    hash += e ;
			    AMD_DEBUG4 ((" e: "ID" hash = "ID"\n",e,hash)) ;
			}
			else
			{
			    /* external degree of e is zero, absorb e into me*/
			    AMD_DEBUG1 ((" Element "ID" =>"ID" (aggressive)\n",
				e, me)) ;
			    ASSERT (dext == 0) ;
			    Pe [e] = FLIP (me) ;
			    W [e] = 0 ;
			}
		    }
		}
	    }
	    else
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			dext = we - wflg ;
			ASSERT (dext >= 0) ;
			deg += dext ;
			Iw [pn++] = e ;
			hash += e ;
			AMD_DEBUG4 (("	e: "ID" hash = "ID"\n",e,hash)) ;
		    }
		}
	    }

	    /* count the number of elements in i (including me): */
	    Elen [i] = pn - p1 + 1 ;

	    /* ------------------------------------------------------------- */
	    /* scan the supervariables in the list associated with i */
	    /* ------------------------------------------------------------- */

	    /* The bulk of the AMD run time is typically spent in this loop,
	     * particularly if the matrix has many dense rows that are not
	     * removed prior to ordering. */
	    p3 = pn ;
	    p4 = p1 + Len [i] ;
	    for (p = p2 + 1 ; p < p4 ; p++)
	    {
		j = Iw [p] ;
		ASSERT (j >= 0 && j < n) ;
		nvj = Nv [j] ;
		if (nvj > 0)
		{
		    /* j is unabsorbed, and not in Lme. */
		    /* add to degree and add to new list */
		    deg += nvj ;
		    Iw [pn++] = j ;
		    hash += j ;
		    AMD_DEBUG4 (("  s: "ID" hash "ID" Nv[j]= "ID"\n",
				j, hash, nvj)) ;
		}
	    }

	    /* ------------------------------------------------------------- */
	    /* update the degree and check for mass elimination */
	    /* ------------------------------------------------------------- */

	    /* with aggressive absorption, deg==0 is identical to the
	     * Elen [i] == 1 && p3 == pn test, below. */
	    ASSERT (IMPLIES (aggressive, (deg==0) == (Elen[i]==1 && p3==pn))) ;

	    if (Elen [i] == 1 && p3 == pn)
	    {

		/* --------------------------------------------------------- */
		/* mass elimination */
		/* --------------------------------------------------------- */

		/* There is nothing left of this node except for an edge to
		 * the current pivot element.  Elen [i] is 1, and there are
		 * no variables adjacent to node i.  Absorb i into the
		 * current pivot element, me.  Note that if there are two or
		 * more mass eliminations, fillin due to mass elimination is
		 * possible within the nvpiv-by-nvpiv pivot block.  It is this
		 * step that causes AMD's analysis to be an upper bound.
		 *
		 * The reason is that the selected pivot has a lower
		 * approximate degree than the true degree of the two mass
		 * eliminated nodes.  There is no edge between the two mass
		 * eliminated nodes.  They are merged with the current pivot
		 * anyway.
		 *
		 * No fillin occurs in the Schur complement, in any case,
		 * and this effect does not decrease the quality of the
		 * ordering itself, just the quality of the nonzero and
		 * flop count analysis.  It also means that the post-ordering
		 * is not an exact elimination tree post-ordering. */

		AMD_DEBUG1 (("  MASS i "ID" => parent e "ID"\n", i, me)) ;
		Pe [i] = FLIP (me) ;
		nvi = -Nv [i] ;
		degme -= nvi ;
		nvpiv += nvi ;
		nel += nvi ;
		Nv [i] = 0 ;
		Elen [i] = EMPTY ;

	    }
	    else
	    {

		/* --------------------------------------------------------- */
		/* update the upper-bound degree of i */
		/* --------------------------------------------------------- */

		/* the following degree does not yet include the size
		 * of the current element, which is added later: */

		Degree [i] = MIN (Degree [i], deg) ;

		/* --------------------------------------------------------- */
		/* add me to the list for i */
		/* --------------------------------------------------------- */

		/* move first supervariable to end of list */
		Iw [pn] = Iw [p3] ;
		/* move first element to end of element part of list */
		Iw [p3] = Iw [p1] ;
		/* add new element, me, to front of list. */
		Iw [p1] = me ;
		/* store the new length of the list in Len [i] */
		Len [i] = pn - p1 + 1 ;

		/* --------------------------------------------------------- */
		/* place in hash bucket.  Save hash key of i in Last [i]. */
		/* --------------------------------------------------------- */

		/* NOTE: this can fail if hash is negative, because the ANSI C
		 * standard does not define a % b when a and/or b are negative.
		 * That's why hash is defined as an unsigned Int, to avoid this
		 * problem. */
		hash = hash % n ;
		ASSERT (((Int) hash) >= 0 && ((Int) hash) < n) ;

		/* if the Hhead array is not used: */
		j = Head [hash] ;
		if (j <= EMPTY)
		{
		    /* degree list is empty, hash head is FLIP (j) */
		    Next [i] = FLIP (j) ;
		    Head [hash] = FLIP (i) ;
		}
		else
		{
		    /* degree list is not empty, use Last [Head [hash]] as
		     * hash head. */
		    Next [i] = Last [j] ;
		    Last [j] = i ;
		}

		/* if a separate Hhead array is used: *
		Next [i] = Hhead [hash] ;
		Hhead [hash] = i ;
		*/

		Last [i] = hash ;
	    }
	}

	Degree [me] = degme ;

	/* ----------------------------------------------------------------- */
	/* Clear the counter array, W [...], by incrementing wflg. */
	/* ----------------------------------------------------------------- */

	/* make sure that wflg+n does not cause integer overflow */
	lemax =  MAX (lemax, degme) ;
	wflg += lemax ;
	wflg = clear_flag (wflg, wbig, W, n) ;
	/*  at this point, W [0..n-1] < wflg holds */

/* ========================================================================= */
/* SUPERVARIABLE DETECTION */
/* ========================================================================= */

	AMD_DEBUG1 (("Detecting supervariables:\n")) ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    AMD_DEBUG2 (("Consider i "ID" nv "ID"\n", i, Nv [i])) ;
	    if (Nv [i] < 0)
	    {
		/* i is a principal variable in Lme */

		/* ---------------------------------------------------------
		 * examine all hash buckets with 2 or more variables.  We do
		 * this by examing all unique hash keys for supervariables in
		 * the pattern Lme of the current element, me
		 * --------------------------------------------------------- */

		/* let i = head of hash bucket, and empty the hash bucket */
		ASSERT (Last [i] >= 0 && Last [i] < n) ;
		hash = Last [i] ;

		/* if Hhead array is not used: */
		j = Head [hash] ;
		if (j == EMPTY)
		{
		    /* hash bucket and degree list are both empty */
		    i = EMPTY ;
		}
		else if (j < EMPTY)
		{
		    /* degree list is empty */
		    i = FLIP (j) ;
		    Head [hash] = EMPTY ;
		}
		else
		{
		    /* degree list is not empty, restore Last [j] of head j */
		    i = Last [j] ;
		    Last [j] = EMPTY ;
		}

		/* if separate Hhead array is used: *
		i = Hhead [hash] ;
		Hhead [hash] = EMPTY ;
		*/

		ASSERT (i >= EMPTY && i < n) ;
		AMD_DEBUG2 (("----i "ID" hash "ID"\n", i, hash)) ;

		while (i != EMPTY && Next [i] != EMPTY)
		{

		    /* -----------------------------------------------------
		     * this bucket has one or more variables following i.
		     * scan all of them to see if i can absorb any entries
		     * that follow i in hash bucket.  Scatter i into w.
		     * ----------------------------------------------------- */

		    ln = Len [i] ;
		    eln = Elen [i] ;
		    ASSERT (ln >= 0 && eln >= 0) ;
		    ASSERT (Pe [i] >= 0 && Pe [i] < iwlen) ;
		    /* do not flag the first element in the list (me) */
		    for (p = Pe [i] + 1 ; p <= Pe [i] + ln - 1 ; p++)
		    {
			ASSERT (Iw [p] >= 0 && Iw [p] < n) ;
			W [Iw [p]] = wflg ;
		    }

		    /* ----------------------------------------------------- */
		    /* scan every other entry j following i in bucket */
		    /* ----------------------------------------------------- */

		    jlast = i ;
		    j = Next [i] ;
		    ASSERT (j >= EMPTY && j < n) ;

		    while (j != EMPTY)
		    {
			/* ------------------------------------------------- */
			/* check if j and i have identical nonzero pattern */
			/* ------------------------------------------------- */

			AMD_DEBUG3 (("compare i "ID" and j "ID"\n", i,j)) ;

			/* check if i and j have the same Len and Elen */
			ASSERT (Len [j] >= 0 && Elen [j] >= 0) ;
			ASSERT (Pe [j] >= 0 && Pe [j] < iwlen) ;
			ok = (Len [j] == ln) && (Elen [j] == eln) ;
			/* skip the first element in the list (me) */
			for (p = Pe [j] + 1 ; ok && p <= Pe [j] + ln - 1 ; p++)
			{
			    ASSERT (Iw [p] >= 0 && Iw [p] < n) ;
			    if (W [Iw [p]] != wflg) ok = 0 ;
			}
			if (ok)
			{
			    /* --------------------------------------------- */
			    /* found it!  j can be absorbed into i */
			    /* --------------------------------------------- */

			    AMD_DEBUG1 (("found it! j "ID" => i "ID"\n", j,i));
			    Pe [j] = FLIP (i) ;
			    /* both Nv [i] and Nv [j] are negated since they */
			    /* are in Lme, and the absolute values of each */
			    /* are the number of variables in i and j: */
			    Nv [i] += Nv [j] ;
			    Nv [j] = 0 ;
			    Elen [j] = EMPTY ;
			    /* delete j from hash bucket */
			    ASSERT (j != Next [j]) ;
			    j = Next [j] ;
			    Next [jlast] = j ;

			}
			else
			{
			    /* j cannot be absorbed into i */
			    jlast = j ;
			    ASSERT (j != Next [j]) ;
			    j = Next [j] ;
			}
			ASSERT (j >= EMPTY && j < n) ;
		    }

		    /* -----------------------------------------------------
		     * no more variables can be absorbed into i
		     * go to next i in bucket and clear flag array
		     * ----------------------------------------------------- */

		    wflg++ ;
		    i = Next [i] ;
		    ASSERT (i >= EMPTY && i < n) ;

		}
	    }
	}
	AMD_DEBUG2 (("detect done\n")) ;

/* ========================================================================= */
/* RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVARIABLES FROM ELEMENT */
/* ========================================================================= */

	p = pme1 ;
	nleft = n - nel ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    nvi = -Nv [i] ;
	    AMD_DEBUG3 (("Restore i "ID" "ID"\n", i, nvi)) ;
	    if (nvi > 0)
	    {
		/* i is a principal variable in Lme */
		/* restore Nv [i] to signify that i is principal */
		Nv [i] = nvi ;

		/* --------------------------------------------------------- */
		/* compute the external degree (add size of current element) */
		/* --------------------------------------------------------- */

		deg = Degree [i] + degme - nvi ;
		deg = MIN (deg, nleft - nvi) ;
		ASSERT (IMPLIES (aggressive, deg > 0) && deg >= 0 && deg < n) ;

		/* --------------------------------------------------------- */
		/* place the supervariable at the head of the degree list */
		/* --------------------------------------------------------- */

		inext = Head [deg] ;
		ASSERT (inext >= EMPTY && inext < n) ;
		if (inext != EMPTY) Last [inext] = i ;
		Next [i] = inext ;
		Last [i] = EMPTY ;
		Head [deg] = i ;

		/* --------------------------------------------------------- */
		/* save the new degree, and find the minimum degree */
		/* --------------------------------------------------------- */

		mindeg = MIN (mindeg, deg) ;
		Degree [i] = deg ;

		/* --------------------------------------------------------- */
		/* place the supervariable in the element pattern */
		/* --------------------------------------------------------- */

		Iw [p++] = i ;

	    }
	}
	AMD_DEBUG2 (("restore done\n")) ;

/* ========================================================================= */
/* FINALIZE THE NEW ELEMENT */
/* ========================================================================= */

	AMD_DEBUG2 (("ME = "ID" DONE\n", me)) ;
	Nv [me] = nvpiv ;
	/* save the length of the list for the new element me */
	Len [me] = p - pme1 ;
	if (Len [me] == 0)
	{
	    /* there is nothing left of the current pivot element */
	    /* it is a root of the assembly tree */
	    Pe [me] = EMPTY ;
	    W [me] = 0 ;
	}
	if (elenme != 0)
	{
	    /* element was not constructed in place: deallocate part of */
	    /* it since newly nonprincipal variables may have been removed */
	    pfree = p ;
	}

	/* The new element has nvpiv pivots and the size of the contribution
	 * block for a multifrontal method is degme-by-degme, not including
	 * the "dense" rows/columns.  If the "dense" rows/columns are included,
	 * the frontal matrix is no larger than
	 * (degme+ndense)-by-(degme+ndense).
	 */

	if (Info != (double *) NULL)
	{
	    f = nvpiv ;
	    r = degme + ndense ;
	    dmax = MAX (dmax, f + r) ;

	    /* number of nonzeros in L (excluding the diagonal) */
	    lnzme = f*r + (f-1)*f/2 ;
	    lnz += lnzme ;

	    /* number of divide operations for LDL' and for LU */
	    ndiv += lnzme ;

	    /* number of multiply-subtract pairs for LU */
	    s = f*r*r + r*(f-1)*f + (f-1)*f*(2*f-1)/6 ;
	    nms_lu += s ;

	    /* number of multiply-subtract pairs for LDL' */
	    nms_ldl += (s + lnzme)/2 ;
	}

#ifndef NDEBUG
	AMD_DEBUG2 (("finalize done nel "ID" n "ID"\n   ::::\n", nel, n)) ;
	for (pme = Pe [me] ; pme <= Pe [me] + Len [me] - 1 ; pme++)
	{
	      AMD_DEBUG3 ((" "ID"", Iw [pme])) ;
	}
	AMD_DEBUG3 (("\n")) ;
#endif

    }

/* ========================================================================= */
/* DONE SELECTING PIVOTS */
/* ========================================================================= */

    if (Info != (double *) NULL)
    {

	/* count the work to factorize the ndense-by-ndense submatrix */
	f = ndense ;
	dmax = MAX (dmax, (double) ndense) ;

	/* number of nonzeros in L (excluding the diagonal) */
	lnzme = (f-1)*f/2 ;
	lnz += lnzme ;

	/* number of divide operations for LDL' and for LU */
	ndiv += lnzme ;

	/* number of multiply-subtract pairs for LU */
	s = (f-1)*f*(2*f-1)/6 ;
	nms_lu += s ;

	/* number of multiply-subtract pairs for LDL' */
	nms_ldl += (s + lnzme)/2 ;

	/* number of nz's in L (excl. diagonal) */
	Info [AMD_LNZ] = lnz ;

	/* number of divide ops for LU and LDL' */
	Info [AMD_NDIV] = ndiv ;

	/* number of multiply-subtract pairs for LDL' */
	Info [AMD_NMULTSUBS_LDL] = nms_ldl ;

	/* number of multiply-subtract pairs for LU */
	Info [AMD_NMULTSUBS_LU] = nms_lu ;

	/* number of "dense" rows/columns */
	Info [AMD_NDENSE] = ndense ;

	/* largest front is dmax-by-dmax */
	Info [AMD_DMAX] = dmax ;

	/* number of garbage collections in AMD */
	Info [AMD_NCMPA] = ncmpa ;

	/* successful ordering */
	Info [AMD_STATUS] = AMD_OK ;
    }

/* ========================================================================= */
/* POST-ORDERING */
/* ========================================================================= */

/* -------------------------------------------------------------------------
 * Variables at this point:
 *
 * Pe: holds the elimination tree.  The parent of j is FLIP (Pe [j]),
 *	or EMPTY if j is a root.  The tree holds both elements and
 *	non-principal (unordered) variables absorbed into them.
 *	Dense variables are non-principal and unordered.
 *
 * Elen: holds the size of each element, including the diagonal part.
 *	FLIP (Elen [e]) > 0 if e is an element.  For unordered
 *	variables i, Elen [i] is EMPTY.
 *
 * Nv: Nv [e] > 0 is the number of pivots represented by the element e.
 *	For unordered variables i, Nv [i] is zero.
 *
 * Contents no longer needed:
 *	W, Iw, Len, Degree, Head, Next, Last.
 *
 * The matrix itself has been destroyed.
 *
 * n: the size of the matrix.
 * No other scalars needed (pfree, iwlen, etc.)
 * ------------------------------------------------------------------------- */

    /* restore Pe */
    for (i = 0 ; i < n ; i++)
    {
	Pe [i] = FLIP (Pe [i]) ;
    }

    /* restore Elen, for output information, and for postordering */
    for (i = 0 ; i < n ; i++)
    {
	Elen [i] = FLIP (Elen [i]) ;
    }

/* Now the parent of j is Pe [j], or EMPTY if j is a root.  Elen [e] > 0
 * is the size of element e.  Elen [i] is EMPTY for unordered variable i. */

#ifndef NDEBUG
    AMD_DEBUG2 (("\nTree:\n")) ;
    for (i = 0 ; i < n ; i++)
    {
	AMD_DEBUG2 ((" "ID" parent: "ID"   ", i, Pe [i])) ;
	ASSERT (Pe [i] >= EMPTY && Pe [i] < n) ;
	if (Nv [i] > 0)
	{
	    /* this is an element */
	    e = i ;
	    AMD_DEBUG2 ((" element, size is "ID"\n", Elen [i])) ;
	    ASSERT (Elen [e] > 0) ;
	}
	AMD_DEBUG2 (("\n")) ;
    }
    AMD_DEBUG2 (("\nelements:\n")) ;
    for (e = 0 ; e < n ; e++)
    {
	if (Nv [e] > 0)
	{
	    AMD_DEBUG3 (("Element e= "ID" size "ID" nv "ID" \n", e,
		Elen [e], Nv [e])) ;
	}
    }
    AMD_DEBUG2 (("\nvariables:\n")) ;
    for (i = 0 ; i < n ; i++)
    {
	Int cnt ;
	if (Nv [i] == 0)
	{
	    AMD_DEBUG3 (("i unordered: "ID"\n", i)) ;
	    j = Pe [i] ;
	    cnt = 0 ;
	    AMD_DEBUG3 (("  j: "ID"\n", j)) ;
	    if (j == EMPTY)
	    {
		AMD_DEBUG3 (("	i is a dense variable\n")) ;
	    }
	    else
	    {
		ASSERT (j >= 0 && j < n) ;
		while (Nv [j] == 0)
		{
		    AMD_DEBUG3 (("	j : "ID"\n", j)) ;
		    j = Pe [j] ;
		    AMD_DEBUG3 (("	j:: "ID"\n", j)) ;
		    cnt++ ;
		    if (cnt > n) break ;
		}
		e = j ;
		AMD_DEBUG3 (("	got to e: "ID"\n", e)) ;
	    }
	}
    }
#endif

/* ========================================================================= */
/* compress the paths of the variables */
/* ========================================================================= */

    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{

	    /* -------------------------------------------------------------
	     * i is an un-ordered row.  Traverse the tree from i until
	     * reaching an element, e.  The element, e, was the principal
	     * supervariable of i and all nodes in the path from i to when e
	     * was selected as pivot.
	     * ------------------------------------------------------------- */

	    AMD_DEBUG1 (("Path compression, i unordered: "ID"\n", i)) ;
	    j = Pe [i] ;
	    ASSERT (j >= EMPTY && j < n) ;
	    AMD_DEBUG3 (("	j: "ID"\n", j)) ;
	    if (j == EMPTY)
	    {
		/* Skip a dense variable.  It has no parent. */
		AMD_DEBUG3 (("      i is a dense variable\n")) ;
		continue ;
	    }

	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		AMD_DEBUG3 (("		j : "ID"\n", j)) ;
		j = Pe [j] ;
		AMD_DEBUG3 (("		j:: "ID"\n", j)) ;
		ASSERT (j >= 0 && j < n) ;
	    }
	    /* got to an element e */
	    e = j ;
	    AMD_DEBUG3 (("got to e: "ID"\n", e)) ;

	    /* -------------------------------------------------------------
	     * traverse the path again from i to e, and compress the path
	     * (all nodes point to e).  Path compression allows this code to
	     * compute in O(n) time.
	     * ------------------------------------------------------------- */

	    j = i ;
	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		jnext = Pe [j] ;
		AMD_DEBUG3 (("j "ID" jnext "ID"\n", j, jnext)) ;
		Pe [j] = e ;
		j = jnext ;
		ASSERT (j >= 0 && j < n) ;
	    }
	}
    }

/* ========================================================================= */
/* postorder the assembly tree */
/* ========================================================================= */

    AMD_postorder (n, Pe, Nv, Elen,
	W,			/* output order */
	Head, Next, Last) ;	/* workspace */

/* ========================================================================= */
/* compute output permutation and inverse permutation */
/* ========================================================================= */

    /* W [e] = k means that element e is the kth element in the new
     * order.  e is in the range 0 to n-1, and k is in the range 0 to
     * the number of elements.  Use Head for inverse order. */

    for (k = 0 ; k < n ; k++)
    {
	Head [k] = EMPTY ;
	Next [k] = EMPTY ;
    }
    for (e = 0 ; e < n ; e++)
    {
	k = W [e] ;
	ASSERT ((k == EMPTY) == (Nv [e] == 0)) ;
	if (k != EMPTY)
	{
	    ASSERT (k >= 0 && k < n) ;
	    Head [k] = e ;
	}
    }

    /* construct output inverse permutation in Next,
     * and permutation in Last */
    nel = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	e = Head [k] ;
	if (e == EMPTY) break ;
	ASSERT (e >= 0 && e < n && Nv [e] > 0) ;
	Next [e] = nel ;
	nel += Nv [e] ;
    }
    ASSERT (nel == n - ndense) ;

    /* order non-principal variables (dense, & those merged into supervar's) */
    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{
	    e = Pe [i] ;
	    ASSERT (e >= EMPTY && e < n) ;
	    if (e != EMPTY)
	    {
		/* This is an unordered variable that was merged
		 * into element e via supernode detection or mass
		 * elimination of i when e became the pivot element.
		 * Place i in order just before e. */
		ASSERT (Next [i] == EMPTY && Nv [e] > 0) ;
		Next [i] = Next [e] ;
		Next [e]++ ;
	    }
	    else
	    {
		/* This is a dense unordered variable, with no parent.
		 * Place it last in the output order. */
		Next [i] = nel++ ;
	    }
	}
    }
    ASSERT (nel == n) ;

    AMD_DEBUG2 (("\n\nPerm:\n")) ;
    for (i = 0 ; i < n ; i++)
    {
	k = Next [i] ;
	ASSERT (k >= 0 && k < n) ;
	Last [k] = i ;
	AMD_DEBUG2 (("   perm ["ID"] = "ID"\n", k, i)) ;
    }
}



void AMD_preprocess
(
    Int n,		/* input matrix: A is n-by-n */
    const Int Ap [ ],	/* size n+1 */
    const Int Ai [ ],	/* size nz = Ap [n] */

    /* output matrix R: */
    Int Rp [ ],		/* size n+1 */
    Int Ri [ ],		/* size nz (or less, if duplicates present) */

    Int W [ ],		/* workspace of size n */
    Int Flag [ ]	/* workspace of size n */
)
{

    /* --------------------------------------------------------------------- */
    /* local variables */
    /* --------------------------------------------------------------------- */

    Int i, j, p, p2 ;

    ASSERT (AMD_valid (n, n, Ap, Ai) != AMD_INVALID) ;

    /* --------------------------------------------------------------------- */
    /* count the entries in each row of A (excluding duplicates) */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	W [i] = 0 ;		/* # of nonzeros in row i (excl duplicates) */
	Flag [i] = EMPTY ;	/* Flag [i] = j if i appears in column j */
    }
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		W [i]++ ;	    /* one more entry in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the row pointers for R */
    /* --------------------------------------------------------------------- */

    Rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Rp [i+1] = Rp [i] + W [i] ;
    }
    for (i = 0 ; i < n ; i++)
    {
	W [i] = Rp [i] ;
	Flag [i] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* construct the row form matrix R */
    /* --------------------------------------------------------------------- */

    /* R = row form of pattern of A */
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		Ri [W [i]++] = j ;  /* put col j in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }
}


void AMD_postorder
(
    /* inputs, not modified on output: */
    Int nn,		/* nodes are in the range 0..nn-1 */
    Int Parent [ ],	/* Parent [j] is the parent of j, or EMPTY if root */
    Int Nv [ ],		/* Nv [j] > 0 number of pivots represented by node j,
			 * or zero if j is not a node. */
    Int Fsize [ ],	/* Fsize [j]: size of node j */

    /* output, not defined on input: */
    Int Order [ ],	/* output post-order */

    /* workspaces of size nn: */
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
)
{
    Int i, j, k, parent, frsize, f, fprev, maxfrsize, bigfprev, bigf, fnext ;

    for (j = 0 ; j < nn ; j++)
    {
	Child [j] = EMPTY ;
	Sibling [j] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* place the children in link lists - bigger elements tend to be last */
    /* --------------------------------------------------------------------- */

    for (j = nn-1 ; j >= 0 ; j--)
    {
	if (Nv [j] > 0)
	{
	    /* this is an element */
	    parent = Parent [j] ;
	    if (parent != EMPTY)
	    {
		/* place the element in link list of the children its parent */
		/* bigger elements will tend to be at the end of the list */
		Sibling [j] = Child [parent] ;
		Child [parent] = j ;
	    }
	}
    }

#ifndef NDEBUG
    {
	Int nels, ff, nchild ;
	AMD_DEBUG1 (("\n\n================================ AMD_postorder:\n"));
	nels = 0 ;
	for (j = 0 ; j < nn ; j++)
	{
	    if (Nv [j] > 0)
	    {
		AMD_DEBUG1 (( ""ID" :  nels "ID" npiv "ID" size "ID
		    " parent "ID" maxfr "ID"\n", j, nels,
		    Nv [j], Fsize [j], Parent [j], Fsize [j])) ;
		/* this is an element */
		/* dump the link list of children */
		nchild = 0 ;
		AMD_DEBUG1 (("    Children: ")) ;
		for (ff = Child [j] ; ff != EMPTY ; ff = Sibling [ff])
		{
		    AMD_DEBUG1 ((ID" ", ff)) ;
		    ASSERT (Parent [ff] == j) ;
		    nchild++ ;
		    ASSERT (nchild < nn) ;
		}
		AMD_DEBUG1 (("\n")) ;
		parent = Parent [j] ;
		if (parent != EMPTY)
		{
		    ASSERT (Nv [parent] > 0) ;
		}
		nels++ ;
	    }
	}
    }
    AMD_DEBUG1 (("\n\nGo through the children of each node, and put\n"
		 "the biggest child last in each list:\n")) ;
#endif

    /* --------------------------------------------------------------------- */
    /* place the largest child last in the list of children for each node */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	if (Nv [i] > 0 && Child [i] != EMPTY)
	{

#ifndef NDEBUG
	    Int nchild ;
	    AMD_DEBUG1 (("Before partial sort, element "ID"\n", i)) ;
	    nchild = 0 ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (f >= 0 && f < nn) ;
		AMD_DEBUG1 (("      f: "ID"  size: "ID"\n", f, Fsize [f])) ;
		nchild++ ;
		ASSERT (nchild <= nn) ;
	    }
#endif

	    /* find the biggest element in the child list */
	    fprev = EMPTY ;
	    maxfrsize = EMPTY ;
	    bigfprev = EMPTY ;
	    bigf = EMPTY ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (f >= 0 && f < nn) ;
		frsize = Fsize [f] ;
		if (frsize >= maxfrsize)
		{
		    /* this is the biggest seen so far */
		    maxfrsize = frsize ;
		    bigfprev = fprev ;
		    bigf = f ;
		}
		fprev = f ;
	    }
	    ASSERT (bigf != EMPTY) ;

	    fnext = Sibling [bigf] ;

	    AMD_DEBUG1 (("bigf "ID" maxfrsize "ID" bigfprev "ID" fnext "ID
		" fprev " ID"\n", bigf, maxfrsize, bigfprev, fnext, fprev)) ;

	    if (fnext != EMPTY)
	    {
		/* if fnext is EMPTY then bigf is already at the end of list */

		if (bigfprev == EMPTY)
		{
		    /* delete bigf from the element of the list */
		    Child [i] = fnext ;
		}
		else
		{
		    /* delete bigf from the middle of the list */
		    Sibling [bigfprev] = fnext ;
		}

		/* put bigf at the end of the list */
		Sibling [bigf] = EMPTY ;
		ASSERT (Child [i] != EMPTY) ;
		ASSERT (fprev != bigf) ;
		ASSERT (fprev != EMPTY) ;
		Sibling [fprev] = bigf ;
	    }

#ifndef NDEBUG
	    AMD_DEBUG1 (("After partial sort, element "ID"\n", i)) ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (f >= 0 && f < nn) ;
		AMD_DEBUG1 (("        "ID"  "ID"\n", f, Fsize [f])) ;
		ASSERT (Nv [f] > 0) ;
		nchild-- ;
	    }
	    ASSERT (nchild == 0) ;
#endif

	}
    }

    /* --------------------------------------------------------------------- */
    /* postorder the assembly tree */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	Order [i] = EMPTY ;
    }

    k = 0 ;

    for (i = 0 ; i < nn ; i++)
    {
	if (Parent [i] == EMPTY && Nv [i] > 0)
	{
	    AMD_DEBUG1 (("Root of assembly tree "ID"\n", i)) ;
	    k = AMD_post_tree (i, k, Child, Sibling, Order, Stack
#ifndef NDEBUG
		, nn
#endif
		) ;
	}
    }
}


Int AMD_post_tree
(
    Int root,			/* root of the tree */
    Int k,			/* start numbering at k */
    Int Child [ ],		/* input argument of size nn, undefined on
				 * output.  Child [i] is the head of a link
				 * list of all nodes that are children of node
				 * i in the tree. */
    const Int Sibling [ ],	/* input argument of size nn, not modified.
				 * If f is a node in the link list of the
				 * children of node i, then Sibling [f] is the
				 * next child of node i.
				 */
    Int Order [ ],		/* output order, of size nn.  Order [i] = k
				 * if node i is the kth node of the reordered
				 * tree. */
    Int Stack [ ]		/* workspace of size nn */
#ifndef NDEBUG
    , Int nn			/* nodes are in the range 0..nn-1. */
#endif
)
{
    Int f, head, h, i ;

#if 0
    /* --------------------------------------------------------------------- */
    /* recursive version (Stack [ ] is not used): */
    /* --------------------------------------------------------------------- */

    /* this is simple, but can cause stack overflow if nn is large */
    i = root ;
    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
    {
	k = AMD_post_tree (f, k, Child, Sibling, Order, Stack, nn) ;
    }
    Order [i] = k++ ;
    return (k) ;
#endif

    /* --------------------------------------------------------------------- */
    /* non-recursive version, using an explicit stack */
    /* --------------------------------------------------------------------- */

    /* push root on the stack */
    head = 0 ;
    Stack [0] = root ;

    while (head >= 0)
    {
	/* get head of stack */
	ASSERT (head < nn) ;
	i = Stack [head] ;
	AMD_DEBUG1 (("head of stack "ID" \n", i)) ;
	ASSERT (i >= 0 && i < nn) ;

	if (Child [i] != EMPTY)
	{
	    /* the children of i are not yet ordered */
	    /* push each child onto the stack in reverse order */
	    /* so that small ones at the head of the list get popped first */
	    /* and the biggest one at the end of the list gets popped last */
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		head++ ;
		ASSERT (head < nn) ;
		ASSERT (f >= 0 && f < nn) ;
	    }
	    h = head ;
	    ASSERT (head < nn) ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (h > 0) ;
		Stack [h--] = f ;
		AMD_DEBUG1 (("push "ID" on stack\n", f)) ;
		ASSERT (f >= 0 && f < nn) ;
	    }
	    ASSERT (Stack [h] == i) ;

	    /* delete child list so that i gets ordered next time we see it */
	    Child [i] = EMPTY ;
	}
	else
	{
	    /* the children of i (if there were any) are already ordered */
	    /* remove i from the stack and order it.  Front i is kth front */
	    head-- ;
	    AMD_DEBUG1 (("pop "ID" order "ID"\n", i, k)) ;
	    Order [i] = k++ ;
	    ASSERT (k <= nn) ;
	}

#ifndef NDEBUG
	AMD_DEBUG1 (("\nStack:")) ;
	for (h = head ; h >= 0 ; h--)
	{
	    Int j = Stack [h] ;
	    AMD_DEBUG1 ((" "ID, j)) ;
	    ASSERT (j >= 0 && j < nn) ;
	}
	AMD_DEBUG1 (("\n\n")) ;
	ASSERT (head < nn) ;
#endif

    }
    return (k) ;
}


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

/* ========================================================================== */
/* === KLU_lsolve =========================================================== */
/* ========================================================================== */

/* Solve Lx=b.  Assumes L is unit lower triangular and where the unit diagonal
 * entry is NOT stored.  Overwrites B  with the solution X.  B is n-by-nrhs
 * and is stored in ROW form with row dimension nrhs.  nrhs must be in the
 * range 1 to 4. */
void KLU_lsolve
(
    /* inputs, not modified: */
    Int n,
    Int Lip [ ],
    Int Llen [ ],
    Unit LU [ ],
    Int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    Entry X [ ]
)
{
    Entry x [4], lik ;
    Int *Li ;
    Entry *Lx ;
    Int k, p, len, i ;

    switch (nrhs)
    {

        case 1:
            for (k = 0 ; k < n ; k++)
            {
                x [0] = X [k] ;
                GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
                /* unit diagonal of L is not stored*/
                for (p = 0 ; p < len ; p++)
                {
                    /* X [Li [p]] -= Lx [p] * x [0] ; */
                    MULT_SUB (X [Li [p]], Lx [p], x [0]) ;
                }
            }
            break ;

        case 2:

            for (k = 0 ; k < n ; k++)
            {
                x [0] = X [2*k    ] ;
                x [1] = X [2*k + 1] ;
                GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
                for (p = 0 ; p < len ; p++)
                {
                    i = Li [p] ;
                    lik = Lx [p] ;
                    MULT_SUB (X [2*i], lik, x [0]) ;
                    MULT_SUB (X [2*i + 1], lik, x [1]) ;
                }
            }
            break ;

        case 3:

            for (k = 0 ; k < n ; k++)
            {
                x [0] = X [3*k    ] ;
                x [1] = X [3*k + 1] ;
                x [2] = X [3*k + 2] ;
                GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
                for (p = 0 ; p < len ; p++)
                {
                    i = Li [p] ;
                    lik = Lx [p] ;
                    MULT_SUB (X [3*i], lik, x [0]) ;
                    MULT_SUB (X [3*i + 1], lik, x [1]) ;
                    MULT_SUB (X [3*i + 2], lik, x [2]) ;
                }
            }
            break ;

        case 4:

            for (k = 0 ; k < n ; k++)
            {
                x [0] = X [4*k    ] ;
                x [1] = X [4*k + 1] ;
                x [2] = X [4*k + 2] ;
                x [3] = X [4*k + 3] ;
                GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
                for (p = 0 ; p < len ; p++)
                {
                    i = Li [p] ;
                    lik = Lx [p] ;
                    MULT_SUB (X [4*i], lik, x [0]) ;
                    MULT_SUB (X [4*i + 1], lik, x [1]) ;
                    MULT_SUB (X [4*i + 2], lik, x [2]) ;
                    MULT_SUB (X [4*i + 3], lik, x [3]) ;
                }
            }
            break ;

    }
}

/* ========================================================================== */
/* === KLU_usolve =========================================================== */
/* ========================================================================== */

/* Solve Ux=b.  Assumes U is non-unit upper triangular and where the diagonal
 * entry is NOT stored.  Overwrites B with the solution X.  B is n-by-nrhs
 * and is stored in ROW form with row dimension nrhs.  nrhs must be in the
 * range 1 to 4. */

void KLU_usolve
(
    /* inputs, not modified: */
    Int n,
    Int Uip [ ],
    Int Ulen [ ],
    Unit LU [ ],
    Entry Udiag [ ],
    Int nrhs,
    /* right-hand-side on input, solution to Ux=b on output */
    Entry X [ ]
)
{
    Entry x [4], uik, ukk ;
    Int *Ui ;
    Entry *Ux ;
    Int k, p, len, i ;

    switch (nrhs)
    {

        case 1:

            for (k = n-1 ; k >= 0 ; k--)
            {
                GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, len) ;
                /* x [0] = X [k] / Udiag [k] ; */
                DIV (x [0], X [k], Udiag [k]) ;
                X [k] = x [0] ;
                for (p = 0 ; p < len ; p++)
                {
                    /* X [Ui [p]] -= Ux [p] * x [0] ; */
                    MULT_SUB (X [Ui [p]], Ux [p], x [0]) ;

                }
            }

            break ;

        case 2:

            for (k = n-1 ; k >= 0 ; k--)
            {
                GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, len) ;
                ukk = Udiag [k] ;
                /* x [0] = X [2*k    ] / ukk ;
                x [1] = X [2*k + 1] / ukk ; */
                DIV (x [0], X [2*k], ukk) ;
                DIV (x [1], X [2*k + 1], ukk) ;

                X [2*k    ] = x [0] ;
                X [2*k + 1] = x [1] ;
                for (p = 0 ; p < len ; p++)
                {
                    i = Ui [p] ;
                    uik = Ux [p] ;
                    /* X [2*i    ] -= uik * x [0] ;
                    X [2*i + 1] -= uik * x [1] ; */
                    MULT_SUB (X [2*i], uik, x [0]) ;
                    MULT_SUB (X [2*i + 1], uik, x [1]) ;
                }
            }

            break ;

        case 3:

            for (k = n-1 ; k >= 0 ; k--)
            {
                GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, len) ;
                ukk = Udiag [k] ;

                DIV (x [0], X [3*k], ukk) ;
                DIV (x [1], X [3*k + 1], ukk) ;
                DIV (x [2], X [3*k + 2], ukk) ;

                X [3*k    ] = x [0] ;
                X [3*k + 1] = x [1] ;
                X [3*k + 2] = x [2] ;
                for (p = 0 ; p < len ; p++)
                {
                    i = Ui [p] ;
                    uik = Ux [p] ;
                    MULT_SUB (X [3*i], uik, x [0]) ;
                    MULT_SUB (X [3*i + 1], uik, x [1]) ;
                    MULT_SUB (X [3*i + 2], uik, x [2]) ;
                }
            }

            break ;

        case 4:

            for (k = n-1 ; k >= 0 ; k--)
            {
                GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, len) ;
                ukk = Udiag [k] ;

                DIV (x [0], X [4*k], ukk) ;
                DIV (x [1], X [4*k + 1], ukk) ;
                DIV (x [2], X [4*k + 2], ukk) ;
                DIV (x [3], X [4*k + 3], ukk) ;

                X [4*k    ] = x [0] ;
                X [4*k + 1] = x [1] ;
                X [4*k + 2] = x [2] ;
                X [4*k + 3] = x [3] ;
                for (p = 0 ; p < len ; p++)
                {
                    i = Ui [p] ;
                    uik = Ux [p] ;

                    MULT_SUB (X [4*i], uik, x [0]) ;
                    MULT_SUB (X [4*i + 1], uik, x [1]) ;
                    MULT_SUB (X [4*i + 2], uik, x [2]) ;
                    MULT_SUB (X [4*i + 3], uik, x [3]) ;
                }
            }

            break ;

    }
}


int KLU_solve
(
    /* inputs, not modified */
    KLU_symbolic *Symbolic,
    KLU_numeric *Numeric,
    Int d,                  /* leading dimension of B */
    Int nrhs,               /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ],           /* size n*nrhs, in column-oriented form, with
                             * leading dimension d. */
    /* --------------- */
    KLU_common *Common
)
{
    Entry x [4], offik, s ;
    double rs, *Rs ;
    Entry *Offx, *X, *Bz, *Udiag ;
    Int *Q, *R, *Pnum, *Offp, *Offi, *Lip, *Uip, *Llen, *Ulen ;
    Unit **LUbx ;
    Int k1, k2, nk, k, block, pend, n, p, nblocks, chunk, nr, i ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Bz = (Entry *) B ;
    n = Symbolic->n ;
    nblocks = Symbolic->nblocks ;  /* number of blocks */
    Q = Symbolic->Q ;              /* size n */
    R = Symbolic->R ;               /* size n+1, but only R [0..nblocks] is used */

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    ASSERT (nblocks == Numeric->nblocks) ;
    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = (Entry *) Numeric->Offx ;

    Lip  = Numeric->Lip ;
    Llen = Numeric->Llen ;
    Uip  = Numeric->Uip ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;
    Udiag = Numeric->Udiag ;

    Rs = Numeric->Rs ;
    X = (Entry *) Numeric->Xwork ;

    /* ---------------------------------------------------------------------- */
    /* solve in chunks of 4 columns at a time */
    /* ---------------------------------------------------------------------- */

    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
    {

        /* ------------------------------------------------------------------ */
        /* get the size of the current chunk */
        /* ------------------------------------------------------------------ */

        nr = MIN (nrhs - chunk, 4) ;

        /* ------------------------------------------------------------------ */
        /* scale and permute the right hand side, X = P*(R\B) */
        /* ------------------------------------------------------------------ */

        if (Rs == NULL)
        {

            /* no scaling */
            switch (nr)
            {

                case 1:

                    for (k = 0 ; k < n ; k++)  
                    {
                        X [k] = Bz [Pnum [k]] ;
                    }
                    break ;

                case 2: 

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        X [2*k    ] = Bz [i      ] ;
                        X [2*k + 1] = Bz  [i + d  ] ;
                    }
                    break ;

                case 3:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        X [3*k    ] = Bz [i      ] ;
                        X [3*k + 1] = Bz [i + d  ] ;
                        X [3*k + 2] = Bz [i + d*2] ;
                    }
                    break ;

                case 4:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        X [4*k    ] = Bz [i      ] ;
                        X [4*k + 1] = Bz [i + d  ] ;
                        X [4*k + 2] = Bz [i + d*2] ;
                        X [4*k + 3] = Bz [i + d*3] ;
                    }
                    break ;
            }

        }
        else
        {

            switch (nr)
            {

                case 1:

                    for (k = 0 ; k < n ; k++)
                    {
                        SCALE_DIV_ASSIGN (X [k], Bz  [Pnum [k]], Rs [k]) ;
                    }
                    break ;

                case 2:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (X [2*k], Bz [i], rs) ;
                        SCALE_DIV_ASSIGN (X [2*k + 1], Bz [i + d], rs) ;
                    }
                    break ;

                case 3:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (X [3*k], Bz [i], rs) ;
                        SCALE_DIV_ASSIGN (X [3*k + 1], Bz [i + d], rs) ;
                        SCALE_DIV_ASSIGN (X [3*k + 2], Bz [i + d*2], rs) ;
                    }
                    break ;

                case 4:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (X [4*k], Bz [i], rs) ;
                        SCALE_DIV_ASSIGN (X [4*k + 1], Bz [i + d], rs) ;
                        SCALE_DIV_ASSIGN (X [4*k + 2], Bz [i + d*2], rs) ;
                        SCALE_DIV_ASSIGN (X [4*k + 3], Bz [i + d*3], rs) ;
                    }
                    break ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* solve X = (L*U + Off)\X */
        /* ------------------------------------------------------------------ */

        for (block = nblocks-1 ; block >= 0 ; block--)
        {

            /* -------------------------------------------------------------- */
            /* the block of size nk is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;

            /* solve the block system */
            if (nk == 1)  //快速求解
            {
                s = Udiag [k1] ;//U的对角线元素的指针
                switch (nr)
                {

                    case 1:
                        DIV (X [k1], X [k1], s) ; //DIV(c,a,b)  { (c) = (a) / (b) ; }
                        break ;

                    case 2:
                        DIV (X [2*k1], X [2*k1], s) ;
                        DIV (X [2*k1 + 1], X [2*k1 + 1], s) ;
                        break ;

                    case 3:
                        DIV (X [3*k1], X [3*k1], s) ;
                        DIV (X [3*k1 + 1], X [3*k1 + 1], s) ;
                        DIV (X [3*k1 + 2], X [3*k1 + 2], s) ;
                        break ;

                    case 4:
                        DIV (X [4*k1], X [4*k1], s) ;
                        DIV (X [4*k1 + 1], X [4*k1 + 1], s) ;
                        DIV (X [4*k1 + 2], X [4*k1 + 2], s) ;
                        DIV (X [4*k1 + 3], X [4*k1 + 3], s) ;
                        break ;

                }
            }
            else
            {
                KLU_lsolve (nk, Lip + k1, Llen + k1, LUbx [block], nr,
                        X + nr*k1) ;
                KLU_usolve (nk, Uip + k1, Ulen + k1, LUbx [block],
                        Udiag + k1, nr, X + nr*k1) ;
            }

            /* -------------------------------------------------------------- */
            /* block back-substitution for the off-diagonal-block entries */
            /* -------------------------------------------------------------- */

            if (block > 0)
            {
                switch (nr)
                {

                    case 1:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [k] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                MULT_SUB (X [Offi [p]], Offx [p], x [0]) ;
                            }
                        }
                        break ;

                    case 2:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [2*k    ] ;
                            x [1] = X [2*k + 1] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
                                offik = Offx [p] ;
                                MULT_SUB (X [2*i], offik, x [0]) ;
                                MULT_SUB (X [2*i + 1], offik, x [1]) ;
                            }
                        }
                        break ;

                    case 3:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [3*k    ] ;
                            x [1] = X [3*k + 1] ;
                            x [2] = X [3*k + 2] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
                                offik = Offx [p] ;
                                MULT_SUB (X [3*i], offik, x [0]) ;
                                MULT_SUB (X [3*i + 1], offik, x [1]) ;
                                MULT_SUB (X [3*i + 2], offik, x [2]) ;
                            }
                        }
                        break ;

                    case 4:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [4*k    ] ;
                            x [1] = X [4*k + 1] ;
                            x [2] = X [4*k + 2] ;
                            x [3] = X [4*k + 3] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
                                offik = Offx [p] ;
                                MULT_SUB (X [4*i], offik, x [0]) ;
                                MULT_SUB (X [4*i + 1], offik, x [1]) ;
                                MULT_SUB (X [4*i + 2], offik, x [2]) ;
                                MULT_SUB (X [4*i + 3], offik, x [3]) ;
                            }
                        }
                        break ;
                }
            }
        }

        /* ------------------------------------------------------------------ */
        /* permute the result, Bz  = Q*X */
        /* ------------------------------------------------------------------ */

        switch (nr)
        {

            case 1:

                for (k = 0 ; k < n ; k++)
                {
                    Bz  [Q [k]] = X [k] ;
                }
                break ;

            case 2:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    Bz  [i      ] = X [2*k    ] ;
                    Bz  [i + d  ] = X [2*k + 1] ;
                }
                break ;

            case 3:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    Bz  [i      ] = X [3*k    ] ;
                    Bz  [i + d  ] = X [3*k + 1] ;
                    Bz  [i + d*2] = X [3*k + 2] ;
                }
                break ;

            case 4:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    Bz  [i      ] = X [4*k    ] ;
                    Bz  [i + d  ] = X [4*k + 1] ;
                    Bz  [i + d*2] = X [4*k + 2] ;
                    Bz  [i + d*3] = X [4*k + 3] ;
                }
                break ;
        }

        /* ------------------------------------------------------------------ */
        /* go to the next chunk of B */
        /* ------------------------------------------------------------------ */

        Bz  += d*4 ;
    }
    return (TRUE) ;
}

/* ========================================================================== */
/* === KLU_refactor ========================================================= */
/* ========================================================================== */

int KLU_refactor        /* returns TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    double Ax [ ],
    KLU_symbolic *Symbolic,

    /* input/output */
    KLU_numeric *Numeric,
    KLU_common  *Common
)
{
    Entry ukk, ujk, s ;
    Entry *Offx, *Lx, *Ux, *X, *Az, *Udiag ;
    double *Rs ;
    Int *Q, *R, *Pnum, *Ui, *Li, *Pinv, *Lip, *Uip, *Llen, *Ulen ;
    Unit **LUbx ;
    Unit *LU ;
    Int k1, k2, nk, k, block, oldcol, pend, oldrow, n, p, newrow, scale,
        nblocks, poff, i, j, up, ulen, llen, maxblock, nzoff ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */
    //common和factor的一致
    Common->numerical_rank = EMPTY ;
    Common->singular_col = EMPTY ;

    Az = (Entry *) Ax ; //

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    n = Symbolic->n ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    Pnum = Numeric->Pnum ;  //最终行置换矩阵
    Offx = (Entry *) Numeric->Offx ; //大小，数值

    LUbx = (Unit **) Numeric->LUbx ; //指向L和U的索引和元素的指针数组

    scale = Common->scale ; //使用得是klu_common 默认得scale = 2

    Rs = Numeric->Rs ; //缩放因子

    Pinv = Numeric->Pinv ; //最终行置换得逆矩阵
    X = (Entry *) Numeric->Xwork ;
    Common->nrealloc = 0 ;
    Udiag = Numeric->Udiag ; //指向U得对角线元素得指针
    nzoff = Symbolic->nzoff ; //非对角块中的非零元素个数

    /* ---------------------------------------------------------------------- */
    /*s 根据一个名为scale的参数来决定是否检查输入矩阵并计算行比例因子Rs。如果scale小于 0，则不进行缩放操作，也不检查输入矩阵；
    如果scale大于等于 0，则检查输入矩阵的索引是否在范围内，并调用KLU_scale函数计算行比例因子，若计算过程中出现问题则返回FALSE表示失败。*/
    /* ---------------------------------------------------------------------- */

    /* do no scale, or check the input matrix, if scale < 0 小于0 就不用缩放 */
    if (scale >= 0)
    {
        /* 如果Ax的值发生变化，则计算Rs缩放因子 */
        if (!KLU_scale (scale, n, Ap, Ai, Ax, Rs, NULL, Common))
        {
            return (FALSE) ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace X */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < maxblock ; k++)
    {
        /* X [k] = 0 */
        CLEAR (X [k]) ;
    }

    poff = 0 ;

    /* ---------------------------------------------------------------------- */
    /* factor each block */
    /* ---------------------------------------------------------------------- */

    if (scale > 0)
    {

        /* ------------------------------------------------------------------ */
        /* scaling */
        /* ------------------------------------------------------------------ */

        for (block = 0 ; block < nblocks ; block++)
        {

            /* -------------------------------------------------------------- */
            /* the block is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;

            if (nk == 1) //单元素的情况
            {

                /* ---------------------------------------------------------- */
                /* singleton case */
                /* ---------------------------------------------------------- */

                oldcol = Q [k1] ;
                pend = Ap [oldcol+1] ;
                CLEAR (s) ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    oldrow = Ai [p] ;
                    newrow = Pinv [oldrow] - k1 ;
                    if (newrow < 0 && poff < nzoff)
                    {
                        /* entry in off-diagonal block */
                        /* Offx [poff] = Az [p] / Rs [oldrow] */
                        //和没缩放处理不同的 是  此时需要除以对应的缩放因子
                        SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]) ;
                        poff++ ;
                    }
                    else
                    {
                        /* singleton */
                        /* s = Az [p] / Rs [oldrow] */
                        SCALE_DIV_ASSIGN (s, Az [p], Rs [oldrow]) ;
                    }
                }
                Udiag [k1] = s ;

            }
            else //缩放因子多块的情况
            {

                /* ---------------------------------------------------------- */
                /* construct and factor the kth block */
                /* ---------------------------------------------------------- */

                Lip  = Numeric->Lip  + k1 ;
                Llen = Numeric->Llen + k1 ;
                Uip  = Numeric->Uip  + k1 ;
                Ulen = Numeric->Ulen + k1 ;
                LU = LUbx [block] ;

                for (k = 0 ; k < nk ; k++)
                {

                    /* ------------------------------------------------------ */
                    /* scatter kth column of the block into workspace X */
                    /* ------------------------------------------------------ */

                    oldcol = Q [k+k1] ;
                    pend = Ap [oldcol+1] ;
                    for (p = Ap [oldcol] ; p < pend ; p++)
                    {
                        oldrow = Ai [p] ;
                        newrow = Pinv [oldrow] - k1 ;
                        if (newrow < 0 && poff < nzoff)
                        {
                            /* entry in off-diagonal part */
                            /* Offx [poff] = Az [p] / Rs [oldrow] */
                            SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]);
                            poff++ ;
                        }
                        else
                        {
                            /* (newrow,k) is an entry in the block */
                            /* X [newrow] = Az [p] / Rs [oldrow] */
                            SCALE_DIV_ASSIGN (X [newrow], Az [p], Rs [oldrow]) ;
                        }
                    }

                    /* ------------------------------------------------------ */
                    /* 计算上三角矩阵U的第k列，并更新原始矩阵（这里记为A）的第k列。同时，它检查是否出现奇异情况，并计算下三角矩阵L的第k列*/
                    /* ------------------------------------------------------ */

                    GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
                    for (up = 0 ; up < ulen ; up++)
                    {
                        j = Ui [up] ;
                        ujk = X [j] ;
                        /* X [j] = 0 */
                        CLEAR (X [j]) ;
                        Ux [up] = ujk ;
                        GET_POINTER (LU, Lip, Llen, Li, Lx, j, llen) ;
                        for (p = 0 ; p < llen ; p++)
                        {
                            /* X [Li [p]] -= Lx [p] * ujk */
                            MULT_SUB (X [Li [p]], Lx [p], ujk) ;
                        }
                    }

                    //获取上三角矩阵的对角元素并检查奇异情况
                    /* get the diagonal entry of U */
                    ukk = X [k] ;
                    /* X [k] = 0 */
                    CLEAR (X [k]) ;
                    /* singular case */
                    if (IS_ZERO (ukk))
                    {
                        /* matrix is numerically singular */
                        Common->status = KLU_SINGULAR ;
                        if (Common->numerical_rank == EMPTY)
                        {
                            Common->numerical_rank = k+k1 ;
                            Common->singular_col = Q [k+k1] ;
                        }
                        if (Common->halt_if_singular)
                        {
                            /* do not continue the factorization */
                            return (FALSE) ;
                        }
                    }
                    Udiag [k+k1] = ukk ;

                    //再次获取下三角矩阵第k列的指针和长度信息。

                    /* gather and divide by pivot to get kth column of L */
                    GET_POINTER (LU, Lip, Llen, Li, Lx, k, llen) ;
                    for (p = 0 ; p < llen ; p++)
                    {
                        i = Li [p] ;
                        DIV (Lx [p], X [i], ukk) ;
                        CLEAR (X [i]) ;
                    }
                }
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /*根据关键行顺序排列比例因子Rs */
    /* ---------------------------------------------------------------------- */

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

    return (TRUE) ;
}



static inline uint64_t rte_rdtsc(void)
{
    union
    {
        uint64_t tsc_64;
        struct
        {
            uint32_t lo_32;
            uint32_t hi_32;
        };
    } tsc;

    asm volatile("rdtsc" : "=a"(tsc.lo_32),
                           "=d"(tsc.hi_32));
    return tsc.tsc_64;
}


/* ========================================================================== */
/* === main()函数  ========================================================== */
/* ========================================================================== */

typedef struct
{
    bool verbose;
    int32_t core;
    int32_t matrix_size;
    int64_t repeat;
    int32_t warm_steps;
} App_ctx;

App_ctx g_app_ctx = {
    .verbose = false,
    .core = 6,
    .matrix_size = 57,

    .warm_steps = 10};


typedef struct {
	double* A;	
	double* original_B;
	int matrix_size;
	int *Ap;
	int *Ai;
	double *Ax;
} tFuncBlockDemoContext;

static tFuncBlockDemoContext* context;

void init_context()
{		
	context = malloc(sizeof(tFuncBlockDemoContext));
	if (!context)
	{
		printf("malloc tFuncBlockDemoContext error");
		exit(-1);
	}
	memset(context, 0x00, sizeof(tFuncBlockDemoContext));
	
    context->matrix_size =  g_app_ctx.matrix_size; //	通过命令行的赋值传递给context->matrix_size 
	
	context->A =  (double *)malloc(context->matrix_size * context->matrix_size * sizeof(double));
	if (!context->A)
	{
		printf("emt_local_malloc A error");
		exit(-1);
	}
	
	context->original_B =  (double *)malloc(context->matrix_size * sizeof(double));
	if (!context->original_B)
	{
		printf("emt_local_malloc original_B error");
		exit(-1);
	}
}

int denseToCscSparse(const double A[], int n, int **Ap, int **Ai, double **Ax)
{
    int nz = 0;
    for (int i = 0; i < n * n; i++)
    {   
        if(fabs(A[i]) >= EPS)
        {   
            nz += 1;
        }
    }

    int *_Ap = (int *)malloc((n + 1) * sizeof(int));     // A的列指针
    int *_Ai = (int *)malloc(nz * sizeof(int));          // A的行索引
    double *_Ax = (double *)malloc(nz * sizeof(double)); // A的非零元素值
	
	if(!_Ap || !_Ai || !_Ax)
	{
		return 0;
	}

    memset(_Ap, 0, (n + 1) * sizeof(int));
    memset(_Ai, 0, nz * sizeof(int));
    memset(_Ax, 0, nz * sizeof(double));

    int k = 0;
    for (int j = 0; j < n; j++) // j为列号
    {
        _Ap[j + 1] += _Ap[j];
        for (int i = 0; i < n; i++) // i为行号
        {
            double value = A[i * n + j]; // 第i行第j列的值
            //记录非0值的三元组信息
            if (fabs(value) >= EPS)
            {
                _Ap[j + 1]++;
                _Ai[k] = i;  //记录行号
                _Ax[k] = value; //记录对应的值
                k++;
            }
        }
    }

    *Ap = _Ap;
    *Ai = _Ai;
    *Ax = _Ax;
	
	return 1;
}

//确保CSV文件德编码格式
int hasBOM(FILE *stream) {
    unsigned char bom[3];
    if (fread(bom, 1, 3, stream) == 3 && bom[0] == 0xEF && bom[1] == 0xBB && bom[2] == 0xBF) {
        return 1;
    } else {
        // Rewind the file pointer to the beginning
        fseek(stream, 0, SEEK_SET);
        return 0;
    }
}

void init_matrix(tFuncBlockDemoContext *context)
{
  
    FILE *stream;
    char *line = NULL ;
    size_t len = 0 ;
    ssize_t nread;
    double val;

   // *************************************************************
    //读取矩阵A 读取向量b
    char filename[64] = {0};
    snprintf(filename, 64, "matrix_%d.csv",context->matrix_size);
    stream = fopen(filename, "r"); //stream 表示文件指针 代表了打开的文件资源 
    if (stream == NULL) {
        printf("fopen matrix.csv error");
        exit(-1);
    }
    if (hasBOM(stream)) 
	{
        // Skip BOM if present
        fseek(stream, 3, SEEK_SET);
    }

    int idx = 0;
   
//通过strtok来分割 指向第一个以逗号分隔的子字符串
     while ((nread = getline(&line, &len, stream))!= -1) 
	{
        char *token = strtok(line, ",");
        while (token!= NULL) 
        {
	    val = atof(token);
            context->A[idx++]= val;
            token = strtok(NULL, ","); //从上次分割后的位置开始查找下一个分割符
        }
    }
   
        free (line);
        fclose(stream);
/*
     for (int i = 0 ; i < idx ; i++)
        {
          printf ("A [%d] = %g\n", i, context->A[i]) ;  
        } 
*/

	if(!denseToCscSparse(context->A, context->matrix_size, &context->Ap, &context->Ai, &context->Ax))
	{
		printf("dense matrix to csc matrix err\n");
		exit(-1);
	}

    //**************************************************************************** 读取向量b //
    FILE *stream_b;
    char *line_b = NULL;
    size_t len_b = 0;
    ssize_t nread_b;
    double val_b;
    int row_b = 0;

    char filename_b[64] = {0};
    snprintf(filename_b,64,"matrix_b_%d.csv",context->matrix_size);
    stream_b = fopen(filename_b,"r");
    if (stream_b == NULL) {
        printf("fopen matrix_b.csv error");
        exit(-1);
    }
       if (hasBOM(stream_b)) 
	{
        // Skip BOM if present
        fseek(stream_b, 3, SEEK_SET);
    }

    while( (nread_b = getline(&line_b,&len_b,stream_b)) != -1 ) //stream_b这个流中有多少行就会循环执行多少次
    {
   
        if( nread_b >0 && line_b[nread_b-1] == '\n' ) //这个line就是getline读取的一行的内容 
        {
            line_b[nread_b-1] = '\0';
        }
         val_b = atof(line_b);
         context->original_B[row_b] = val_b;
         row_b++;
    }
/*
      for (int i = 0 ; i < row_b ; i++)
        {
          printf ("B [%d] = %g\n", i,  context->original_B[i]) ;  
        } 
*/
    free(line_b);
    fclose(stream_b);

    //**************************************************************************** 读取普通矩阵则不需要设值 //
	 int nnz = context->Ap[context->matrix_size];
	 for(int i = 0; i < nnz; i++)
	 {
	 	context->Ax[i] = DRAND(); //对CSC格式中的非0元素数组重新随机赋值
         }
   
     for(int i = 0; i < context->matrix_size; i++)
    {
       context->original_B[i] = DRAND(); //对向量B数据随机赋值
    }
   
   //**************************************************************************** 读取普通矩阵则不需要设值 //

}

int main( int argc,char **argv )
{
	rsn_sched_setscheduler(SCHED_FIFO, 99);
	//通过命令行先给 matrix_size的参数
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-m") == 0)
        {
        if (i + 1 < argc) 
            {
                g_app_ctx.matrix_size = atoi(argv[i + 1]);
                break;
            }
         else
            {
                fprintf(stderr, "Error: --matrix_size option requires an argument.\n");
                return 1;
            }
        }
    }

    init_context(); //再将预设的值传递
    init_matrix(context);

    uint64_t time_begin;

    uint64_t time_factorize;
    uint64_t time_used_factorize;
    uint64_t total_time_factorize = 0;
    uint64_t factorize_min_time = 10000000000000ULL;
    uint64_t factorize_max_time = 0;

    uint64_t time_end;
    uint64_t time_used_solve;
    uint64_t total_time_solve = 0;
    uint64_t min_time = 10000000000000ULL;
    uint64_t max_time = 0;

    int repeat = 100000;
    int warm = 1000;
    double CPU_GHZ = 2.7;

  /*
    context->matrix_size = 5;
    int p [ ] = {0, 2, 5, 9, 10, 12} ;
    context->Ap = p ;
    int    i [ ] = { 0,  1,  0,   2,  4,  1,  2,  3,   4,  2,  1,  4} ;
    context->Ai = i ;
    double x [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
    context->Ax = x;
    double b [ ] = {8., 45., -3., 3., 19.} ;
    context->original_B = b ;
*/

    KLU_symbolic *Symbolic ;
    klu_common Common ;
    klu_numeric *Numeric ;

    klu_defaults (&Common) ;

    Symbolic = klu_analyze (context->matrix_size, context->Ap, context->Ai, &Common) ;
	
    Numeric = klu_factor (context->Ap, context->Ai, context->Ax, Symbolic, &Common) ;

    double * tmp = malloc(context->matrix_size * sizeof(double)); //tmp是这段分配空间的起始地址 i则可以代表地址偏移
    memcpy(tmp,context->original_B,context->matrix_size * sizeof(double)); //把original_B地址的内容cp到tmp这个地址
	
	int idx;
    for (idx = 0 ; idx < repeat ; idx++)
    {
		for (int i = 0;i<context->matrix_size;i++)                                                                                          
		{
			//tmp[i]表示基于首地址tmp的偏移量                                                                                        
			context->original_B[i]= tmp[i];      //如果是  context->original_B[i]= tmp[i] 则是将b向量还原                                                                                      
		}

		time_begin = rte_rdtsc();

		klu_refactor (context->Ap, context->Ai, context->Ax, Symbolic, Numeric, &Common) ;
		time_factorize = rte_rdtsc();

		klu_solve (Symbolic, Numeric, context->matrix_size, 1, context->original_B , &Common) ;
		time_end = rte_rdtsc();

		//预热阶段不统计
		if (idx >= warm)
		{
			time_used_factorize = (time_factorize - time_begin) / CPU_GHZ;
			total_time_factorize += time_used_factorize;
			if (time_used_factorize > factorize_max_time)
				factorize_max_time = time_used_factorize;
			if (time_used_factorize < factorize_min_time)
				factorize_min_time = time_used_factorize;

			time_used_solve = (time_end - time_factorize) / CPU_GHZ;
			total_time_solve += time_used_solve;
			if (time_used_solve > max_time)
				max_time = time_used_solve;
			if (time_used_solve < min_time)
				min_time = time_used_solve;
		}
		if (idx == repeat -1)
		{
			for (int i = 0 ; i < context->matrix_size ; i++)
				printf ("x [%d] = %g\n", i,  context->original_B[i]) ;
		}
    }

    klu_free_symbolic (&Symbolic, &Common) ;
    klu_free_numeric (&Numeric, &Common) ;

    printf(" refactor-time:%luns, refactor_min: %luns, refactor_max: %luns, refactor_avg: %luns\n",time_used_factorize,factorize_min_time,factorize_max_time,total_time_factorize / repeat );

    printf(" solve-time:%luns, min: %luns, max: %luns, avg: %luns\n",time_used_solve,min_time,max_time,total_time_solve / repeat );


    return(0);
}





