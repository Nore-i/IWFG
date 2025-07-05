/*****************************************************************************
 *                                                                           *
 *                  Small (Matlab/Octave) Toolbox for Kriging                *
 *                                                                           *
 * Copyright Notice                                                          *
 *                                                                           *
 *    Copyright (C) 2015-2017 CentraleSupelec                                *
 *                                                                           *
 *    Author:  Julien Bect  <julien.bect@centralesupelec.fr>                 *
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

#include "stk_mex.h"
#include <stdlib.h>
#include "wfg.h"
#include "iwfg.h"

typedef struct {
    double   val;   /* objective value */
    mwIndex  idx;   /* 1-based MATLAB row index */
} pair_t;

static int cmp_pair_asc(const void *a, const void *b)
/* return <0 if a<b, >0 if a>b (ascending order) */
{
    double d = ((const pair_t*)a)->val - ((const pair_t*)b)->val;
    return (d < 0.0) ? -1 : (d > 0.0);
}

/* We assume in this file that OBJECTIVE is the same as double */

static mxArray*
bottom_k_indices(mxArray* f,   /* M-by-N MATLAB matrix  */
                 FRONT*   buffer,  /* already initialised   */
                 int32_T        k)
{
  size_t i, j;           /* loop indices */
  size_t nb_points;      /* number of points */
  size_t nb_objectives;  /* number of objectives */
  double *data;          /* pointer to input data */
  int *indices;        /* indices of bottom-k contributors */
  int result = 0;     /* result of iwfg_bottom_k() */
  

  nb_points = mxGetM (f);
  if (nb_points == 0)
      mexErrMsgTxt("No points in the front.");

  nb_objectives = mxGetN (f);
  data = mxGetPr (f);

  indices = mxMalloc(nb_points, sizeof(int));
  if (indices != NULL)
    for (i = 0; i < nb_points; ++i)
        indices[i] = (int)i + 1;  /* MATLAB is 1-based */


  if (nb_objectives == 0)
    {
      mexErrMsgTxt("No objectives for points in the front.");
    }
  else if (nb_objectives == 1)
    {
      /* one objective: sort and return k minimum */
      pair_t *pairs = (pair_t*)mxMalloc(nb_points * sizeof(pair_t));
      for (i = 0; i < nb_points; ++i) {
          pairs[i].val = data[i];    /* column-major; only one column */
          pairs[i].idx = (mwIndex)(i + 1);   /* MATLAB is 1-based */
      }

      /* ---------- sort by objective value (ascending) ----------- */
      qsort(pairs, nb_points, sizeof(pair_t), cmp_pair_asc);

      /* ---------- saturate k and create output ------------------ */
      if (k > (int32_T)nb_points) k = (int32_T)nb_points;

      mxArray *idx_arr = mxCreateNumericMatrix(k, 1, mxINT32_CLASS, mxREAL);
      int32_T *idx_data = (int32_T*)mxGetData(idx_arr);

      for (i = 0; i < (size_t)k; ++i)
          idx_data[i] = (int32_T)pairs[i].idx;

      mxFree(pairs);
      return idx_arr;           /* <-- done for 1-D case */
    }
  else /* two ore more objectives */
    {

      wfg_front_resize (buffer, nb_points, nb_objectives);

      for (i = 0; i < nb_points; i++)
        for (j = 0; j < nb_objectives; j++)
          buffer->points[i].objectives[j] = data[j * nb_points + i];

      result = iwfg_bottom_k(buffer, k, indices);
      if (result < 0)
      {
          mxFree(indices);
          mexErrMsgIdAndTxt("stk:iwfg:bottom_k_indices",
                            "Error in iwfg_bottom_k: %d", result);
      }

      mxArray *idx_arr = mxCreateNumericMatrix(k, 1, mxINT32_CLASS, mxREAL);
      int32_T *idx_data = (int32_T*)mxGetData(idx_arr);
      for (i = 0; i < (size_t)k; ++i)
          idx_data[i] = (int32_T)indices[i];

      mxFree(indices);

      return idx_arr;
    }

  // FRONT buf = *buf_proto;   /* shallow copy of sizes / pointers     */
  // mwSize m  = mxGetM(front_mx);
  // if (k > (int32_T)m) k = (int32_T)m;       /* saturate            */
  // mxArray* idx = mxCreateNumericMatrix(k, 1, mxINT32_CLASS, mxREAL);
  // int32_T* p   = (int32_T*)mxGetData(idx);
  // for (int32_T i = 0; i < k; ++i) p[i] = i+1;   /* <- placeholder  */
  // return idx;
}

#define Y_IN     prhs[0]
#define K_IN     prhs[1]
#define IDX_OUT  plhs[0]

/* ---------------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  /* ----------- variables shared with the old code --------------- */
  size_t i;                    /* loop index over fronts             */
  mxArray **fronts;            /* array of N mxArray* (each a front)  */
  FRONT buffer;                /* reusable WFG FRONT                  */
  size_t nb_fronts;            /* #fronts (1 or numel(cell))          */
  size_t nb_points, nb_objectives;
  size_t maxm = 0, maxn = 0;
  bool   must_free_fronts;

  /* ---------------------- new:  K -------------------------------- */
  int32_T k;    /* number of indices to return per front             */

  /* ######################## argument checks ##################### */

  /* ------ #lhs / #rhs ------------------------------------------ */
  if (nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");
  if (nrhs != 2)
      mexErrMsgTxt("Exactly two input arguments expected: Y and K.");

  /* ------ prhs[0]  → fronts[] ---------------------------------- */
  if (stk_is_realmatrix(Y_IN))                 /* single matrix front */
  {
      nb_fronts      = 1;
      fronts         = (mxArray**)prhs;        /* alias */
      must_free_fronts = false;
  }
  else if (mxIsCell(Y_IN))                     /* cell array of fronts */
  {
      nb_fronts      = mxGetNumberOfElements(Y_IN);
      fronts         = (mxArray**)mxMalloc(nb_fronts * sizeof(mxArray*));
      for (i = 0; i < nb_fronts; ++i)
          fronts[i]  = mxGetCell(Y_IN, i);
      must_free_fronts = true;
  }
  else
      mexErrMsgTxt("Y must be a real matrix or a cell array of matrices.");

  /* ------ prhs[1]  →  k (scalar int32) -------------------------- */
  if (mxIsInt32(K_IN) && mxGetNumberOfElements(K_IN) == 1)
  {
      k = *((int32_T*)mxGetData(K_IN));
  }
  else if (mxIsDouble(K_IN) && mxGetNumberOfElements(K_IN) == 1)
  {
      double tmp = mxGetScalar(K_IN);
      if (tmp < 1.0  || floor(tmp) != tmp)
          mexErrMsgTxt("K must be a positive integer.");
      k = (int32_T)tmp;
  }
  else
      mexErrMsgTxt("K must be a scalar int32.");

  if (k <= 0)
      mexErrMsgTxt("K must be strictly positive.");

  /* ########### figure out max #points / #objectives ############# */

  for (i = 0; i < nb_fronts; ++i)
  {
      nb_points      = mxGetM(fronts[i]);
      nb_objectives  = mxGetN(fronts[i]);
      if (nb_points      > maxm) maxm = nb_points;
      if (nb_objectives  > maxn) maxn = nb_objectives;

      if (k > nb_points)
          mexErrMsgTxt("K must not be larger than the number of points in any front.");
  }

  /* Reusable buffer compatible with the largest front we will see */
  wfg_front_init(&buffer, (int)maxm, (int)maxn);
  wfg_alloc((int)maxm, (int)maxn);

  /* ########################   main loop   ####################### */

  if (nb_fronts == 1)                  /* single front → single vector */
  {
      IDX_OUT = bottom_k_indices(fronts[0], &buffer, k);
  }
  else                                 /* cell array in = cell array out */
  {
      IDX_OUT = mxCreateCellArray(mxGetNumberOfDimensions(Y_IN),
                                  mxGetDimensions(Y_IN));
      for (i = 0; i < nb_fronts; ++i)
          mxSetCell(IDX_OUT, i, bottom_k_indices(fronts[i], &buffer, k));
  }

  /* ########################  clean-up  ########################### */

  wfg_free((int)maxm, (int)maxn);
  wfg_front_destroy(&buffer);
  if (must_free_fronts) mxFree(fronts);
}

