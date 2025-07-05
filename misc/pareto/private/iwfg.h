/*****************************************************************************
 *                                                                           *
 *                  Small (Matlab/Octave) Toolbox for Kriging                *
 *                                                                           *
 * Copyright Notice                                                          *
 *                                                                           *
 *    Copyright (C) 2015, 2016 CentraleSupelec                               *
 *                                                                           *
 *    Author:  Julien Bect  <julien.bect@centralesupelec.fr>                 *
 *                                                                           *
 *    Based on the file wfg.h from WFG 1.10 by Lyndon While, Lucas           *
 *    Bradstreet, Luigi Barone, released under the GPLv2+ licence. The       *
 *    original copyright notice is:                                          *
 *                                                                           *
 *       Copyright (C) 2010 Lyndon While, Lucas Bradstreet                   *
 *                                                                           *
 * Copying Permission Statement                                              *
 *                                                                           *
 *    This file is part of                                                   *
 *                                                                           *
 *            STK: a Small (Matlab/Octave) Toolbox for Kriging               *
 *               (https://github.com/stk-kriging/stk/)                   *
 *                                                                           *
 *    STK is free software: you can redistribute it and/or modify it under   *
 *    the terms of the GNU General Public License as published by the Free   *
 *    Software Foundation,  either version 3  of the License, or  (at your   *
 *    option) any later version.                                             *
 *                                                                           *
 *    STK is distributed  in the hope that it will  be useful, but WITHOUT   *
 *    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY   *
 *    or FITNESS  FOR A  PARTICULAR PURPOSE.  See  the GNU  General Public   *
 *    License for more details.                                              *
 *                                                                           *
 *    You should  have received a copy  of the GNU  General Public License   *
 *    along with STK.  If not, see <http://www.gnu.org/licenses/>.           *
 *                                                                           *
 ****************************************************************************/

#ifndef ___IWFG_H___
#define ___IWFG_H___

#include "wfg.h"

typedef struct
{
    int nPoints;
    int n;
    double base_volume;
    POINT *points;
}
SLICE;

typedef struct 
{
    int idx;   /* The index of the point in the original front */
    POINT *point;   /* The excluded point in question */
    SLICE *slices; /* Currently associated slices */
    int nslices;
    int nslices_alloc;
    int done_slices; /* The number of slices already processed */
    double exclusive_hv; /* The exclusive hypervolume contribution of this point */
}
EXPOINT;

int exgreater (const void *v1, const void *v2);

int slgreater (const void *v1, const void *v2);

static void insert_on_k_desc(POINT *ql, int *sz,
                             const POINT *p, int k);
                          
static int cmp_on_k_desc(const void *a, const void *b, void *kp);

int iwfg_bottom_k(const FRONT *in_front, int k, int *idx_out);

static void add_slice(SLICE *slice,
                      POINT   *ql,  int ql_size,
                      double   base, int n);

static void ensure_slice_capacity(EXPOINT *ep, int capacity);

static void slice_free (SLICE *sl);

static void compute_slice(FRONT *buffer, EXPOINT *ep);

#endif
