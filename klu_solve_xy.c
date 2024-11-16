
#include "klu.h"

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
