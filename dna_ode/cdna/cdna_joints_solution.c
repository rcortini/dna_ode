/* dna_ode: simulating DNA using the ODE physics engine
 *
 * Copyright (C) 2014, 2015  Ruggero Cortini

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "dna_ode_core.h"
#include "dna.h"
#include "cdna.h"



/* alloc all the memory necessary for the joint solution operations */
ccdna_joint_solution_workspace * ccdna_joint_solution_workspace_alloc (const unsigned int N) {
  ccdna_joint_solution_workspace *ws = (ccdna_joint_solution_workspace *) malloc (sizeof (*ws));

  ws->t = (t_vector *) malloc (N*sizeof (t_vector));
  ws->rhs = (t_vector *) malloc (N*sizeof (t_vector));
  ws->Hii = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->Hij = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->lambda = (t_vector *) malloc (N*sizeof (t_vector));
  ws->buffer = (t_vector *) malloc (N*sizeof (t_vector));
  ws->Dii = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->Dii_inv = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->Lij = (t_matrix *) malloc ((N-1)*sizeof (t_matrix));
  ws->LNj = (t_matrix *) malloc ((N-1)*sizeof (t_matrix));

  return ws;
}



/* init the joint solution workspace */
void ccdna_joint_solution_workspace_init (cdna_simcontext *cdna_simcon) {
  (void) cdna_simcon;
}



/* free all the memory associated with the joint solution workspace */
void ccdna_joint_solution_workspace_free (ccdna_joint_solution_workspace *ws) {
  free (ws->t);
  free (ws->rhs);
  free (ws->Hii);
  free (ws->Hij);
  free (ws->lambda);
  free (ws->buffer);
  free (ws->Dii);
  free (ws->Dii_inv);
  free (ws->Lij);
  free (ws->LNj);
  free (ws);
}



/* macros to allow for compact writing */
#define WS_T(x) (cdna_simcon->ws->t [x])
#define WS_RHS(x) (cdna_simcon->ws->rhs [x])
#define WS_HII(x) (cdna_simcon->ws->Hii [x])
#define WS_HIJ(x) (cdna_simcon->ws->Hij [x])
#define WS_LAMBDA(x) (cdna_simcon->ws->lambda [x])
#define WS_BUFFER(x) (cdna_simcon->ws->buffer [x])
#define WS_DII(x) (cdna_simcon->ws->Dii [x])
#define WS_Dinv(x) (cdna_simcon->ws->Dii_inv [x])
#define WS_LIJ(x) (cdna_simcon->ws->Lij [x])
#define WS_LNJ(x) (cdna_simcon->ws->LNj [x])
#define DNA_SEG(x) (cdna_simcon->dna->seg [x])



/* this function calculates the exact solution to the problem of the DNA constraints. It is invoked only
 * if the dSpaceCollide function returned that there are no collisions (simcon->numc==0). */
void ccdna_joints_solution (simcontext *simcon, cdna_simcontext *cdna_simcon) {
  int c, c_next;
  const int N = cdna_simcon->dna->nsegments;
  t_matrix prev_T, next_T;
  t_matrix LDLt, matrix_buff;
  t_vector buff, buff2;
  dna_segment *prev_seg, *next_seg;
  const t_real ks = PAR(erp_glob)/(PAR(ode_step)*PAR(ode_step));
  const t_real kd = 1./PAR(ode_step);
  const t_real cfm = PAR(cfm_glob);

  /* cycle on all DNA joints */
  next_seg = DNA_SEG (0);
  segment_next_t (WS_T (0), DNA_SEG (0));
  matrix_vector_cross (next_T, WS_T (0));
  for (c=0; c<N; c++) {
    /* "circularity" indices: if c+1==N, then c_next = 0, otherwise it is equal to c+1 */
    c_next = (c+1)%N;

    /* initialize the bodies */
    prev_seg = next_seg;
    next_seg = DNA_SEG (c_next);

    /* fetch t_(c+1) */
    if (c!=N-1) segment_next_t (WS_T (c_next), next_seg);

    /* right-hand side of equation */
    segment_rhs (WS_RHS (c), prev_seg, next_seg, WS_T (c), WS_T (c_next), ks, kd);

    /* calculate T matrix */
    matrix_copy (prev_T, next_T);
    matrix_vector_cross (next_T, WS_T (c_next));

    /* calculating H = J * M^-1 * J^T */
    segment_H (WS_HII (c), WS_HIJ (c), prev_seg, next_seg, prev_T, next_T, cfm, kd);
  }

  /* A = L D L^T decomposition */
  matrix_copy (WS_DII (0), WS_HII (0));
  dCopyVector3 (WS_BUFFER (0), WS_RHS (0));
  matrix_inverse (WS_Dinv (0), WS_DII (0));
  matrix_product (WS_LNJ (0), WS_HIJ (N-1), WS_Dinv (0));
  matrix_vector_product (buff, WS_LNJ (0), WS_BUFFER (0));
  vector_diff (WS_BUFFER (N-1), WS_RHS (N-1), buff);
  matrix_ABAt (LDLt, WS_LNJ (0), WS_DII (0));
  matrix_diff (WS_DII (N-1), WS_HII (N-1), LDLt);

  /* forward substitution */
  for (c=1; c<N-1; c++) {
    matrix_product (WS_LIJ (c-1), WS_HIJ (c-1), WS_Dinv (c-1));
    matrix_ABAt (LDLt, WS_LIJ (c-1), WS_DII (c-1));
    matrix_diff (WS_DII (c), WS_HII (c), LDLt);
    matrix_inverse (WS_Dinv (c), WS_DII (c));
    if (c<N-2) {
      t_matrix HDinv, minus_LNJ;
      matrix_scale (minus_LNJ, WS_LNJ (c-1),-1.);
      matrix_product (HDinv, WS_HIJ (c-1), WS_Dinv (c));
      matrix_product (WS_LNJ (c), minus_LNJ, HDinv);
    }
    else {
      t_matrix LH, H_minus_LH;
      matrix_product (LH, WS_LNJ (c-1), WS_HIJ (c-1));
      matrix_diff (H_minus_LH, WS_HIJ (N-2), LH);
      matrix_product (WS_LNJ (c), H_minus_LH, WS_Dinv (c));
    }
    matrix_ABAt (LDLt, WS_LNJ (c), WS_DII (c));
    matrix_diff (WS_DII (N-1), WS_DII (N-1), LDLt);

    /* L Y' = X */
    matrix_vector_product (buff, WS_LIJ (c-1), WS_BUFFER (c-1));
    vector_diff (WS_BUFFER (c), WS_RHS (c), buff);
    matrix_vector_product (buff, WS_LNJ (c), WS_BUFFER (c));
    vector_inplace_diff (WS_BUFFER (N-1), buff);
  }

  /* back substitution */
  matrix_inverse (WS_Dinv (N-1), WS_DII (N-1));
  matrix_vector_product (WS_LAMBDA (N-1), WS_Dinv (N-1), WS_BUFFER (N-1));
  matrix_inverse (WS_Dinv (N-2), WS_DII (N-2));
  matrix_vector_product (buff, WS_Dinv (N-2), WS_BUFFER (N-2));

  matrix_transpose (matrix_buff, WS_LNJ (N-2));
  matrix_vector_product (buff2, matrix_buff, WS_LAMBDA (N-1));
  vector_diff (WS_LAMBDA (N-2), buff, buff2);

  for (c=N-3; c>=0; c--) {
    matrix_vector_product (buff, WS_Dinv (c), WS_BUFFER (c));
    matrix_transpose (matrix_buff, WS_LIJ (c));
    matrix_vector_product (buff2, matrix_buff, WS_LAMBDA (c+1));
    vector_diff (WS_LAMBDA (c), buff, buff2);
    matrix_transpose (matrix_buff, WS_LNJ (c));
    matrix_vector_product (buff, matrix_buff, WS_LAMBDA (N-1));
    vector_inplace_diff (WS_LAMBDA (c), buff);
  }

  /* generalised forces */
  for (c=0; c<N; c++) {
    t_vector fi;
    c_next = (c+1)%N;

    /* add forces and torques */
    body_add_force (DNA_SEG (c)->b_id, WS_LAMBDA (c));
    vector_cross (buff, WS_T (c), WS_LAMBDA (c));
    body_add_torque (DNA_SEG (c)->b_id, buff);

    dCopyVector3 (fi, WS_LAMBDA (c));
    dNegateVector3 (fi);
    body_add_force (DNA_SEG (c_next)->b_id, fi);
    vector_cross (buff, WS_T (c_next), WS_LAMBDA (c));
    body_add_torque (DNA_SEG (c_next)->b_id, buff);
  }
}
