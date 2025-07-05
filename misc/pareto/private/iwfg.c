#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "mex.h"
#include "wfg.h"
#include "iwfg.h"
#include "pq.h"

#define INITIAL_SLICE_CAPACITY 16

static int cmp_on_k_desc(const void *a, const void *b, void *kp)
{
    int k = *(const int*)kp;
    const POINT *p = a, *q = b;
    double diff = q->objectives[k] - p->objectives[k];   /* q – p (desc) */
    return (diff > 0.0) - (diff < 0.0);                  /* signum        */
}

/* Insert p into ql[0 .. *sz-1] so that
 *   – ql stays sorted descending on objective k;
 *   – any element dominated by p (on objectives 0..k) is discarded.
 * On return *sz is incremented by exactly 1 plus the #kept elements.         */
static void insert_on_k_desc(POINT *ql, int *sz,
                             const POINT *p, int k)
{
    /* ---------- 1.  find insertion point (first item that DOESN’T beat p[k]) */
    int ins = 0;
    while (ins < *sz && BEATS(ql[ins].objectives[k], p->objectives[k]))
        ++ins;

    /* ---------- 2.  shift tail right to make room for p --------------------- */
    for (int j = *sz; j > ins; --j)
        ql[j] = ql[j - 1];

    /* ---------- 3.  insert p ------------------------------------------------- */
    ql[ins] = *p;
    int new_size = ins + 1;

    /* ---------- 4.  append remaining tail, skipping dominated ones ---------- */
    for (int j = ins + 1; j <= *sz; ++j)            /* original elements */
        if (!dominates1way(*p, ql[j], k))
            ql[new_size++] = ql[j];                 /* keep if not dominated */

    *sz = new_size;                                 /* updated logical size */
}

int exgreater (const void *v1, const void *v2)
/* this sorts points worsening in the last objective */
{
    int i;
    POINT p = *(((EXPOINT*)v1)->point);
    POINT q = *(((EXPOINT*)v2)->point);
    int n = p.n; /* assume that p and q have the same size */

    for (i = n - 1; i >= 0; i--)
        if BEATS(p.objectives[i],q.objectives[i]) return -1;
        else if BEATS(q.objectives[i],p.objectives[i]) return  1;

    return 0;
}

int slgreater (const void *v1, const void *v2)
/* this sorts slices worsening base volume */
{
    SLICE p = *((SLICE*)v1);
    SLICE q = *((SLICE*)v2);

    if (p.base_volume > q.base_volume)
        return -1;
    else 
        return  1;
}

int iwfg_bottom_k(const FRONT *in_front, int k, int *idx_out)
{
    if (!in_front || !idx_out || k <= 0)
        return -1;

    const int m = in_front->nPoints;
    const int n = in_front->n;
    if (m == 0) return -2;
    if (k > m)  k = m;

    // wfg_alloc(m, n);

    printf("IWFG\n-------------------\n# of points: %d\nObjective dimension: %d\nBottom k: %d\n-------------------\n", m, n, k);
    printf("\n");
    FRONT slice_buf;
    wfg_front_init(&slice_buf, m, n-1);

    EXPOINT *expoints;
    expoints = (EXPOINT *)mxMalloc(m * sizeof(EXPOINT));
    
    /* Initialize excluded points */
    for (int i = 0; i < m; ++i) {
        expoints[i].idx = i+1;
        expoints[i].exclusive_hv = 0.0;
        expoints[i].point = &(in_front->points[i]);
        expoints[i].nslices = 0;
        expoints[i].nslices_alloc = INITIAL_SLICE_CAPACITY;
        expoints[i].slices = (SLICE *)mxMalloc(INITIAL_SLICE_CAPACITY * sizeof(SLICE));
        expoints[i].done_slices = 0;
    }

    /* Sort in the last objective */
    qsort(expoints, m, sizeof(EXPOINT), exgreater);

    /* Create the initial slice */
    for (int i = 0; i < m; ++i)
    {
        /* allocate & fill pl with **all points except i** */
        POINT *pl = (POINT*)mxMalloc((m - 1) * sizeof(POINT));
        int     p  = 0;
        for (int j = 0; j < m; ++j)
            if (j != i)                         /* skip the excluded point */
                pl[p++] = *(expoints[j].point);  /* copy that row */

        expoints[i].slices[0].points     = pl;
        expoints[i].slices[0].nPoints    = m - 1;
        expoints[i].slices[0].n          = n;
        expoints[i].slices[0].base_volume = 1.0;

        ++expoints[i].nslices;
    }

    /* Create the initial slices along the last objective */
    const int last_obj = n - 1;           /* last objective */

    for (int i = 0; i < m; ++i)
    {
        SLICE *root  = &expoints[i].slices[0];      /* {(1 , pl)}        */
        POINT *pl    = root->points;                /* the list  pl      */
        int    np    = root->nPoints;

        /* pl is already sorted worsening in obj n-1 */
        /* prepare empty ql ---------------------------------------------------- */
        POINT *ql = (POINT*)mxMalloc(np * sizeof(POINT));
        int     ql_size = 0;

        double v = expoints[i].point->objectives[last_obj];   /* z[k]           */
        int     dominated = 0;

        SLICE *s_prime = (SLICE*)mxMalloc((m) * sizeof(SLICE));
        int     s_prime_size = 0;

        for (int j = 0; j < np && !dominated; ++j)
        {
            POINT p = pl[j];
            double pk = p.objectives[last_obj];

            /* ---- if v beats p[k] → close a slice -------------------------- */
            if (BEATS(v, pk)) {
                add_slice(&s_prime[s_prime_size], ql, ql_size, fabs(v - pk), n);
                s_prime_size++;
                v = pk;
            }

            /* ---- insert p in ql, keeping it ordered on k-1 ---------------- */
            if (last_obj > 0)
                insert_on_k_desc(ql, &ql_size, &p, last_obj - 1);
            else                 /* last_obj==0 ⇒ no earlier objective */
                ql[ql_size++] = p;

            /* ---- dominance test vs. z on objectives ≥ k-1 ----------------- */
            dominated = dominates1way(p, *(expoints[i].point),
                                    (last_obj > 0) ? last_obj - 1 : 0);
        }

        /* ---- tail slice down to reference point (0) if not dominated ---- */
        if (!dominated) {
            add_slice(&s_prime[s_prime_size], ql, ql_size, fabs(v - 0.0), n);
            s_prime_size++;
        }

        mxFree(ql);

        /* Update slices for objective n-1 */
        ensure_slice_capacity(&expoints[i], s_prime_size);

        mxFree(root->points);
        mxFree(expoints[i].slices);                   // Free previous slice container
        expoints[i].slices = s_prime;                 // Transfer ownership
        expoints[i].nslices_alloc = s_prime_size;     // Update allocation size

        expoints[i].nslices = s_prime_size;
        s_prime = NULL; // Avoid double free

        /* Sort slices by their base volume */
        qsort(expoints[i].slices, expoints[i].nslices, sizeof(SLICE), slgreater);

        /* Remove last objective of expoint */
        expoints[i].point->n = n - 1;  /* reduce the number of objectives */

        compute_slice(&slice_buf, &expoints[i]);
    }

    printf("Printing expoint array:\n");

    for (int u = 0; u < m; ++u) {
        printf("Expoint %d (", expoints[u].idx);
        printf("ex_hv = %lf, done_slices = %d): ",
            expoints[u].exclusive_hv, expoints[u].done_slices);
        printf("{");
        for (int v = 0; v < expoints[u].point->n; ++v) {
            printf("%lf ", expoints[u].point[0].objectives[v]);
        }
        printf("}\n");
    }
    printf("\n");


    PQ pq;
    pq_init(&pq, m);

    /* build heap */
    for (int i = 0; i < m; ++i)
        pq_push(&pq, &expoints[i], i);

    /* main loop: extract k smallest contributors --------------------------- */
    int n_found = 0;
    while (n_found < k && pq_size(&pq) > 0)
    {
        EXPOINT *ep = pq_pop(&pq);           /* least contributor so far      */

        if (ep->done_slices == ep->nslices) {           /* finished point ?   */
            idx_out[n_found++] = ep->idx;               /* store 1-based idx  */
        } else {
            compute_slice(&slice_buf, ep);
            /* its contribution grew, push it back with new key */
            pq_push(&pq, ep, (int) (ep - expoints));            /* uses current hv */
        }
    }

    pq_free(&pq);


    for (int i = 0; i < m; ++i) {
        /* slice 0 owns only the POINT array, not the OBJECTIVE buffers        */
        if (expoints[i].nslices > 0 && expoints[i].slices[0].points) {
            printf("Freeing: %p\n", expoints[i].slices[0].points);
            mxFree(expoints[i].slices[0].points);
        }

        /* all the *other* slices were deep-copied with add_slice()            */
        for (int j = 1; j < expoints[i].nslices; ++j)
            slice_free(&expoints[i].slices[j]);

        mxFree(expoints[i].slices);        /* free slice container             */
    }

    wfg_front_destroy(&slice_buf);
    // wfg_free(m, n); 
    mxFree(expoints);                             /* finally the expoint list  */
    

    return 0;
}

static void compute_slice(FRONT *buffer, EXPOINT *ep)
{    
    int np_slice = ep->slices[ep->done_slices].nPoints;
    int dim_slice = ep->slices[ep->done_slices].n;
    double ex_hv;


    if (np_slice + 1 > buffer->nPoints_alloc || dim_slice > buffer->n_alloc)
        mexErrMsgTxt("internal error: slice_buf capacity exceeded");


    printf("Processing top slice for expoint %d: ", ep->idx);
    printf("{");
    for (int v = 0; v < ep->point->n; ++v) {
        printf("%lf ", ep->point[0].objectives[v]);
    }
    printf("}\n");
    printf("Slice (base_volume, point list, #points, #objectives):\n");
    printf(" (%lf, ", ep->slices[ep->done_slices].base_volume);
    printf("{");
    for (int u = 0; u < np_slice; ++u) {
        printf("(");
        for (int v = 0; v < dim_slice; ++v) {
            printf("%lf ", ep->slices[ep->done_slices].points[u].objectives[v]);
        }
        printf("),");
    }
    printf("}, %d, %d)\n", np_slice, dim_slice);


    wfg_front_resize (buffer, np_slice + 1, dim_slice);

    for (int u = np_slice + 1; u < buffer->nPoints; ++u)
        buffer->points[u].n = dim_slice; /* set the number of objectives */

    for (int u = 0; u < np_slice; u++) {
        for (int v = 0; v < dim_slice; v++)
            buffer->points[u].objectives[v] = ep->slices[ep->done_slices].points[u].objectives[v];
    }



    for (int v = 0; v < dim_slice; v++)
        buffer->points[np_slice].objectives[v] = ep->point->objectives[v];

    if (np_slice == 0)
    {
        ex_hv = inclhv(*ep->point);  /* if no slice, use the point itself */
    } else 
    {
        ex_hv = exclhv (buffer, np_slice);
    }

    
    printf("Exclusive hypervolume for this slice: %lf\n", ex_hv * ep->slices[ep->done_slices].base_volume);
    printf("\n");

    ep->exclusive_hv += ex_hv * ep->slices[ep->done_slices].base_volume;  /* store exclusive hypervolume */
    // slice_free(&ep->slices[ep->done_slices]);
    ep->done_slices++;  /* increment done slices */
}


static void ensure_slice_capacity(EXPOINT *ep, int capacity)
{
    if (capacity > ep->nslices_alloc) {
        ep->nslices_alloc = capacity;
        ep->slices = (SLICE*)mxRealloc(ep->slices,
                        capacity * sizeof(SLICE));
    }
}

static void add_slice(SLICE   *sl,
                      POINT   *ql,      /* input points to copy      */
                      int      ql_size, /* #points in ql             */
                      double   base,    /* slice thickness           */
                      int      n)       /* original #objectives      */
{
    if (ql_size <= 0) {                 /* nothing to copy           */
        sl->points      = NULL;
        sl->nPoints     = 0;
        sl->n           = n - 1;
        sl->base_volume = base;
        return;
    }

    /* allocate array of POINTs inside the slice ------------------- */
    sl->points      = (POINT*)mxMalloc(ql_size * sizeof(POINT));
    sl->nPoints     = ql_size;
    sl->n           = n - 1;            /* slice lives in n-1 dims   */
    sl->base_volume = base;

    /* deep-copy each point, trimming the last objective ----------- */
    for (int j = 0; j < ql_size; ++j)
    {
        sl->points[j].n = n - 1;        /* new dimensionality        */

        /* allocate objectives array for this point */
        sl->points[j].objectives =
            (OBJECTIVE*)mxMalloc((n - 1) * sizeof(OBJECTIVE));

        /* copy the first n-1 coordinates */
        memcpy(sl->points[j].objectives,
               ql[j].objectives,
               (n - 1) * sizeof(OBJECTIVE));
    }
}

static void slice_free (SLICE *sl)
{
    if (!sl->points) return;                 /* empty slice */
    for (int j = 0; j < sl->nPoints; ++j)
        mxFree(sl->points[j].objectives);
    mxFree(sl->points);
    sl->points = NULL;
}