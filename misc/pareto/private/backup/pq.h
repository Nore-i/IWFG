/* ---------- pq.h ------------------------------------------------------- */
#ifndef PQ_H
#define PQ_H

struct EXPOINT;

#include "iwfg.h"

typedef struct {
    EXPOINT **data;
    int      *pos;       /* pos[idx_in_expoints] = heap position or 0 */
    int       size, cap;
} PQ;

/* initialise a PQ able to hold up to <cap> elements                       */
int  pq_init(PQ *pq, int cap);

/* release memory held by pq                                              */
void pq_free(PQ *pq);

/* current number of elements                                             */
static inline int pq_size(const PQ *pq) { return pq->size; }

/* insert <ep> (whose exclusive_hv is already up-to-date)                  */
int  pq_push(PQ *pq, EXPOINT *ep, int idx_in_expoints);

/* extract element with smallest exclusive_hv; returns NULL if empty       */
EXPOINT *pq_pop(PQ *pq);

/* raise the key of element <ep> to <new_hv>  (new_hv  ≥ old value)        */
void pq_increase_key(PQ *pq, int idx_in_expoints, double new_hv);

#endif
