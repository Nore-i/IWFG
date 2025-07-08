/* test_iwfg.c – drive iwfg_bottom_k() without MATLAB
 *
 * Build:
 *   gcc -std=c11 -I. test_iwfg.c iwfg.c wfg.c pq.c -lm -o iwfg_test
 */

#include "mex_stub.h"             /* must precede any source that #includes <mex.h> */
#include <stdio.h>
#include <math.h>

#include "wfg.h"
#include "iwfg.h"

/* -------------------------------------------
 *  Problem constants (compile-time literals)
 * ------------------------------------------- */
#define M_POINTS   8              /* #points   */
#define N_OBJS     2              /* #objectives */
#define K_BOTTOM   2              /* bottom-k to extract */

/* Convenience macro: WFG = maximise → store –y */
#define NEG(a) (-(a))

int main(void)
{
    /* ----------- test data (fixed-size array) ------------------- */
    double y[M_POINTS][N_OBJS] = {
        { -1.0915, -1.0945 },
        { -0.8615, -1.2125 },
        { -1.1357, -0.8069 },
        { -0.2163, -2.0744 },
        { -0.8615, -1.8341 },
        { -0.7525, -2.0353 },
        { -0.2000, -2.2195 },
        { -0.8615, -1.6057 }
    };

    /* ---------- global scratch allocation required by WFG ------- */
    // wfg_alloc(M_POINTS, N_OBJS);

    /* ---------- build FRONT ------------------------------------- */
    FRONT front;
    wfg_front_init(&front, M_POINTS, N_OBJS);

    for (int i = 0; i < M_POINTS; ++i)
        for (int j = 0; j < N_OBJS; ++j)
            front.points[i].objectives[j] = NEG(y[i][j]);   /* maximise */

    /* ---------- call IWFG --------------------------------------- */
    int idx[K_BOTTOM];     /* 1-based indices */

    int rc = iwfg_bottom_k(&front, K_BOTTOM, idx);
    if (rc != 0) {
        fprintf(stderr, "iwfg_bottom_k() returned error code %d\n", rc);
        return EXIT_FAILURE;
    }

    /* ---------- show result ------------------------------------- */
    printf("Bottom-%d hyper-volume contributors (1-based indices):", K_BOTTOM);
    for (int i = 0; i < K_BOTTOM; ++i) printf(" %d", idx[i]);
    puts("");

    /* ---------- tidy up ----------------------------------------- */
    wfg_front_destroy(&front);
    // wfg_free(M_POINTS, N_OBJS);

    return 0;
}
