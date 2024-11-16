#include"klu.h"

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
    double amd_Info [20], lnz, lnz1, flops, flops1 ;
    Int k1, k2, nk, k, block, oldcol, pend, newcol, result, pc, p, newrow,
        maxnz, nzoff, cstats [20], ok, err = KLU_INVALID ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* compute the inverse of Pbtf */

    for (k = 0 ; k < n ; k++)
    {
        assert (Pbtf [k] >= 0 && Pbtf [k] < n) ;
        Pinv [Pbtf [k]] = k ;
    }

    nzoff = 0 ;
    lnz = 0 ;
    maxnz = 0 ;
    flops = 0 ;
    Symbolic->symmetry = EMPTY ;        /* only computed by AMD */

    /* ---------------------------------------------------------------------- */
    /* order each block */
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
                    assert (newrow < k2) ;
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
        else if (ordering == 0) //使用的是AMD排序算法
        {

            /* -------------------------------------------------------------- */
            /* order the block with AMD (C+C') */
            /* -------------------------------------------------------------- */

            result = amd_order (nk, Cp, Ci, Pblk, NULL, amd_Info) ;

            ok = (result >= 0) ;
            if (result == -1)
            {
                err = KLU_OUT_OF_MEMORY ;
            }

            /* account for memory usage in AMD */
            Common->mempeak = MAX (Common->mempeak,
                Common->memusage + amd_Info [7]) ;

            /* get the ordering statistics from AMD */
            lnz1 = (Int) (amd_Info [9]) + nk ;
            flops1 = 2 * amd_Info [12] + amd_Info [10] ;
            if (pc == maxnz)
            {
                /* get the symmetry of the biggest block */
                Symbolic->symmetry = amd_Info [3] ;
            }

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

        // PRINTF (("Pblk, 1-based:\n")) ;
        for (k = 0 ; k < nk ; k++)
        {
            assert (k + k1 < n) ;
            assert (Pblk [k] + k1 < n) ;
            Q [k + k1] = Qbtf [Pblk [k] + k1] ;
        }
        for (k = 0 ; k < nk ; k++)
        {
            assert (k + k1 < n) ;
            assert (Pblk [k] + k1 < n) ;
            P [k + k1] = Pbtf [Pblk [k] + k1] ;
        }
    }

    // PRINTF (("nzoff %d  Ap[n] %d\n", nzoff, Ap [n])) ;
    assert (nzoff >= 0 && nzoff <= Ap [n]) ;

    /* return estimates of # of nonzeros in L including diagonal */
    Symbolic->lnz = lnz ;           /* EMPTY if COLAMD used */
    Symbolic->unz = lnz ;
    Symbolic->nzoff = nzoff ;
    Symbolic->est_flops = flops ;   /* EMPTY if COLAMD or user-ordering used */
    return (KLU_OK) ;
}


// nblocks 总是为1  Common->ordering = 0  Common->status =0
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
    /* allocate the Symbolic object, and check input matrix */
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
    // if (ordering == 1) //ordering == 0
    // {
    //     /* COLAMD */
    //     Cilen = COLAMD_recommended (nz, n, n) ;
    // }
    if (ordering == 0 || (ordering == 3 && Common->user_order != NULL))
    {
        /* AMD or user ordering function */
        Cilen = nz+1 ;
    }
    //else
    // {
    //     /* invalid ordering */
    //     Common->status = KLU_INVALID ;
    //     KLU_free_symbolic (&Symbolic, Common) ;
    //     return (NULL) ;
    // }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace for BTF permutation */
    /* ---------------------------------------------------------------------- */

    Pbtf = KLU_malloc (n, sizeof (Int), Common) ;
    Qbtf = KLU_malloc (n, sizeof (Int), Common) ;

    /* ---------------------------------------------------------------------- */
    /* get the common parameters for BTF and ordering method */
    /* ---------------------------------------------------------------------- */

    do_btf = Common->btf ; //Common->btf=1
    do_btf = (do_btf) ? TRUE : FALSE ;
    Symbolic->ordering = ordering ;
    Symbolic->do_btf = do_btf ;
    Symbolic->structural_rank = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* find the block triangular form (if requested) */
    /* ---------------------------------------------------------------------- */

    Common->work = 0 ;
//do_btf= 0 即没进行预排序
    // if (do_btf)
    // {
    //     Work = KLU_malloc (5*n, sizeof (Int), Common) ;
    //     // if (Common->status < KLU_OK)
    //     // {
    //     //     /* out of memory */
    //     //     KLU_free (Pbtf, n, sizeof (Int), Common) ;
    //     //     KLU_free (Qbtf, n, sizeof (Int), Common) ;
    //     //     KLU_free_symbolic (&Symbolic, Common) ;
    //     //     return (NULL) ;
    //     // }
        
    //    // nblocks = 1; //一直为1
    //     nblocks = BTF_order (n, Ap, Ai, Common->maxwork, &work, Pbtf, Qbtf, R,
    //             &(Symbolic->structural_rank), Work) ;
    //     Common->structural_rank = Symbolic->structural_rank ;
    //     Common->work += work ;

    //     KLU_free (Work, 5*n, sizeof (Int), Common) ;

    //     /* unflip Qbtf if the matrix does not have full structural rank */
    //     if (Symbolic->structural_rank < n)
    //     {
    //         for (k = 0 ; k < n ; k++)
    //         {
    //             Qbtf [k] = BTF_UNFLIP (Qbtf [k]) ;
    //         }
    //     }

    //     /* find the size of the largest block */
    //     maxblock = 1 ;
    //     for (block = 0 ; block < nblocks ; block++)
    //     {
    //         k1 = R [block] ;
    //         k2 = R [block+1] ;
    //         nk = k2 - k1 ;
    //         PRINTF (("block %d size %d\n", block, nk)) ;
    //         maxblock = MAX (maxblock, nk) ;
    //     }
    // }

    if(!do_btf)
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

    // printf (("maxblock size %d\n", maxblock)) ;
    Symbolic->maxblock = maxblock ;

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, for analyze_worker */
    /* ---------------------------------------------------------------------- */

    Pblk = KLU_malloc (maxblock, sizeof (Int), Common) ;
    Cp   = KLU_malloc (maxblock + 1, sizeof (Int), Common) ;
    Ci   = KLU_malloc (MAX (Cilen, nz+1), sizeof (Int), Common) ;
    Pinv = KLU_malloc (n, sizeof (Int), Common) ;

    /* ---------------------------------------------------------------------- */
    /* order each block of the BTF ordering, and a fill-reducing ordering */
    /* ---------------------------------------------------------------------- */

    if (Common->status == KLU_OK)
    {
     
        //Common->status = 0 说明analyze_worker成功 对给定矩阵进行分块处理和排序操作
      
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

    // if (Common->status < KLU_OK)
    // {
    //     KLU_free_symbolic (&Symbolic, Common) ;
    // }
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
       
     order_and_analyze (n, Ap, Ai, Common) ;
}



int main()
{
    int    n = 5 ;
    int    Ap [ ] = {0, 2, 5, 9, 10, 12} ;
    int    Ai [ ] = { 0,  1,  0,   2,  4,  1,  2,  3,   4,  2,  1,  4} ;
    double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
    double b [ ] = {8., 45., -3., 3., 19.} ;
    klu_symbolic *Symbolic ;
    klu_common Common ;
    klu_defaults (&Common) ;
    Symbolic = klu_analyze (n, Ap, Ai, &Common) ;

}