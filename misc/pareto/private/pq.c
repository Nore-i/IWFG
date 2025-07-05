/* ---------- pq.c ------------------------------------------------------- */
#include <stdlib.h>
#include <string.h>
#include "pq.h"
#include "iwfg.h"      /* for EXPOINT                                      */

#define PARENT(i)   ((i) >> 1)
#define LEFT(i)     ((i) << 1)
#define RIGHT(i)    (LEFT(i) + 1)

/* ------------ helpers -------------------------------------------------- */
static void swap(PQ *pq, int i, int j)
{
    EXPOINT *tmp = pq->data[i];
    pq->data[i]  = pq->data[j];
    pq->data[j]  = tmp;
    pq->pos[pq->data[i]->idx - 1] = i;
    pq->pos[pq->data[j]->idx - 1] = j;
}

static void sift_up(PQ *pq, int i)
{
    while (i > 1 && pq->data[i]->exclusive_hv < pq->data[PARENT(i)]->exclusive_hv)
    {
        swap(pq, i, PARENT(i));
        i = PARENT(i);
    }
}

static void sift_down(PQ *pq, int i)
{
    for (;;)
    {
        int l = LEFT(i), r = RIGHT(i), smallest = i;

        if (l <= pq->size &&
            pq->data[l]->exclusive_hv < pq->data[smallest]->exclusive_hv)
            smallest = l;

        if (r <= pq->size &&
            pq->data[r]->exclusive_hv < pq->data[smallest]->exclusive_hv)
            smallest = r;

        if (smallest == i) break;
        swap(pq, i, smallest);
        i = smallest;
    }
}

/* ------------ public interface ---------------------------------------- */
int pq_init(PQ *pq, int cap)
{
    pq->data = malloc((cap + 1) * sizeof(EXPOINT*)); /* +1 because 1-based */
    pq->pos  = malloc( cap        * sizeof(int));
    if (!pq->data || !pq->pos) return -1;
    pq->size = 0;
    pq->cap  = cap;
    memset(pq->pos, 0, cap * sizeof(int));
    return 0;
}

void pq_free(PQ *pq)
{
    free(pq->data);
    free(pq->pos);
    pq->data = NULL; pq->pos = NULL; pq->size = pq->cap = 0;
}

int pq_push(PQ *pq, EXPOINT *ep, int idx)
{
    if (pq->size == pq->cap) return -1;            /* full */
    int i = ++pq->size;
    pq->data[i] = ep;
    pq->pos[idx] = i;
    sift_up(pq, i);
    return 0;
}

EXPOINT *pq_pop(PQ *pq)
{
    if (pq->size == 0) return NULL;
    EXPOINT *min = pq->data[1];
    pq->pos[min - pq->data[1]] = 0;

    pq->data[1] = pq->data[pq->size];
    pq->pos[pq->data[1] - pq->data[1]] = 1;
    pq->size--;
    sift_down(pq, 1);
    return min;
}

void pq_increase_key(PQ *pq, int idx, double new_hv)
{
    int i = pq->pos[idx];
    if (i == 0) return;                 /* already popped */
    pq->data[i]->exclusive_hv = new_hv; /* key increased  */
    sift_down(pq, i);                  /* only need down heapify */
}
