#include "klu.h"

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
    
    //这段判断得意义并不大
    if (scale > 0) //大于0则说明存在缩放
    {
        /* factorization was not scaled, but refactorization is scaled */
        if (Numeric->Rs == NULL)
        {
            Numeric->Rs = KLU_malloc (n, sizeof (double), Common) ; //为缩放因子数组rs分配内存
            if (Common->status < KLU_OK)
            {
                Common->status = KLU_OUT_OF_MEMORY ;
                return (FALSE) ;
            }
        }
    }
    else
    {
        /* no scaling for refactorization; ensure Numeric->Rs is freed.  This
         * does nothing if Numeric->Rs is already NULL. */
        Numeric->Rs = KLU_free (Numeric->Rs, n, sizeof (double), Common) ;
    }


    Rs = Numeric->Rs ;

    Pinv = Numeric->Pinv ;
    X = (Entry *) Numeric->Xwork ;
    Common->nrealloc = 0 ;
    Udiag = Numeric->Udiag ;
    nzoff = Symbolic->nzoff ;

    /* ---------------------------------------------------------------------- */
    /*s 根据一个名为scale的参数来决定是否检查输入矩阵并计算行比例因子Rs。如果scale小于 0，则不进行缩放操作，也不检查输入矩阵；
    如果scale大于等于 0，则检查输入矩阵的索引是否在范围内，并调用KLU_scale函数计算行比例因子，若计算过程中出现问题则返回FALSE表示失败。*/
    /* ---------------------------------------------------------------------- */

    /* do no scale, or check the input matrix, if scale < 0 小于0 就不用缩放 */
    if (scale >= 0)
    {
        /* 检查超出范围的索引，但不检查重复项 */
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

    if (scale <= 0)
    {

        /* ------------------------------------------------------------------ */
        /* no scaling */
        /* ------------------------------------------------------------------ */

        for (block = 0 ; block < nblocks ; block++)
        {

            /* -------------------------------------------------------------- */
            /* the block is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;

            if (nk == 1) //单元素块
            {

                /* ---------------------------------------------------------- */
                /* 单元素块的情况*/
                /* ---------------------------------------------------------- */

                oldcol = Q [k1] ;
                pend = Ap [oldcol+1] ;
                CLEAR (s) ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    newrow = Pinv [Ai [p]] - k1 ;
                    if (newrow < 0 && poff < nzoff) //表示元素在非对角快中
                    {
                        /* entry in off-diagonal block */
                        Offx [poff] = Az [p] ;
                        poff++ ;
                    }
                    else
                    {
                        /* singleton */
                        s = Az [p] ;
                    }
                }
                Udiag [k1] = s ;

            }
            else
            {

                /* ---------------------------------------------------------- */
                /* construct and factor the kth block */
                /* ---------------------------------------------------------- */

                Lip  = Numeric->Lip  + k1 ;
                Llen = Numeric->Llen + k1 ;
                Uip  = Numeric->Uip  + k1 ;
                Ulen = Numeric->Ulen + k1 ;
                LU = LUbx [block] ;

                for (k = 0 ; k < nk ; k++) //遍历块的每一列
                {

                    /* ------------------------------------------------------ */
                    /* scatter kth column of the block into workspace X */
                    /* ------------------------------------------------------ */

                    oldcol = Q [k+k1] ; //通过列置换数组Q和块起始索引k1确定当前列在原始矩阵中的对应列索引
                    pend = Ap [oldcol+1] ; //确定该列的非零元素结束位置
                    
                    //
                    for (p = Ap [oldcol] ; p < pend ; p++)
                    {
                        newrow = Pinv [Ai [p]] - k1 ;
                        if (newrow < 0 && poff < nzoff)
                        {
                            /* entry in off-diagonal block */
                            Offx [poff] = Az [p] ;
                            poff++ ;
                        }
                        else
                        {
                            /* (newrow,k) is an entry in the block */
                            X [newrow] = Az [p] ;
                        }
                    }

                    /* ------------------------------------------------------ */
                    /* 计算上三角矩阵U的第k列，并在这个过程中更新矩阵A的第k列 */
                    /* ------------------------------------------------------ */

                    GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
                    //遍历上三角的第K列
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

                    //如果IS_ZERO(ukk)，表示矩阵可能是奇异的：
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

                    //计算下三角的第K列
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
    else  //进行缩放的情况
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
