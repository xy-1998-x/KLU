
#include "klu.h"

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

int main()
{
  /*  KLU_common *Common;
    KLU_defaults (Common) ;
    printf("tol: %f\n", Common->tol);
    printf("memgrow: %f\n",Common->memgrow);
    printf("initmem_amd: %f\n", Common->initmem_amd);
    printf("initmem: %f\n",Common->initmem);
    printf("maxwork: %f\n", Common->maxwork);
    printf("btf: %d\n", Common->btf);
    printf("ordering: %d\n", Common->ordering);
    printf("scale: %d\n", Common->scale);
    printf("user_order: %p\n", Common->user_order);
    printf("user_data: %p\n", Common->user_data);
    printf("halt_if_singular: %d\n", Common->halt_if_singular);
    printf("status: %d\n", Common->status);
    printf("nrealloc: %d\n", Common->nrealloc);
    printf("structural_rank: %d\n", Common->structural_rank);
    printf("numerical_rank: %d\n", Common->numerical_rank);
    printf("singular_col: %d\n", Common->singular_col);
    printf("noffdiag: %d\n", Common->noffdiag);
    printf("flops: %f\n", Common->flops);
    printf("rcond: %f\n", Common->rcond);
    printf("condest: %f\n", Common->condest);
    printf("rgrowth: %f\n", Common->rgrowth);
    printf("work: %f\n", Common->work);
    printf("memusage: %zu\n", Common->memusage);
    printf("mempeak: %zu\n", Common->mempeak);

    这么写 会有段错误：因为初始化一个指针的时候 没给指针分配内存而使用 这样会导致段错误
    而只是初始化变量的时候 系统会自动分配内存 然后使用   KLU_defaults (&Common)应用变量的地址
     则不会出现段错误                             */

//这种写法不会有段错误
    KLU_common Common;
    KLU_defaults (&Common) ;

    printf("tol: %f\n", Common.tol);
    printf("memgrow: %f\n",Common.memgrow);
    printf("initmem_amd: %f\n", Common.initmem_amd);
    printf("initmem: %f\n",Common.initmem);
    printf("maxwork: %f\n", Common.maxwork);
    printf("btf: %d\n", Common.btf);
    printf("ordering: %d\n", Common.ordering);
    printf("scale: %d\n", Common.scale);
    printf("user_order: %p\n", Common.user_order);
    printf("user_data: %p\n", Common.user_data);
    printf("halt_if_singular: %d\n", Common.halt_if_singular);
    printf("status: %d\n", Common.status);
    printf("nrealloc: %d\n", Common.nrealloc);
    printf("structural_rank: %d\n", Common.structural_rank);
    printf("numerical_rank: %d\n", Common.numerical_rank);
    printf("singular_col: %d\n", Common.singular_col);
    printf("noffdiag: %d\n", Common.noffdiag);
    printf("flops: %f\n", Common.flops);
    printf("rcond: %f\n", Common.rcond);
    printf("condest: %f\n", Common.condest);
    printf("rgrowth: %f\n", Common.rgrowth);
    printf("work: %f\n", Common.work);
    printf("memusage: %zu\n", Common.memusage);
    printf("mempeak: %zu\n", Common.mempeak);

    return 0;
}