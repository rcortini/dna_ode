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
#include "fdna.h"


/* alloc all the memory necessary for the joint solution operations */
fdna_joint_solution_workspace * fdna_joint_solution_workspace_alloc (const unsigned int N) {
  fdna_joint_solution_workspace *ws = (fdna_joint_solution_workspace *) malloc (sizeof (*ws));

  ws->t = (t_vector *) malloc (N*sizeof (t_vector));
  ws->rhs = (t_vector *) malloc (N*sizeof (t_vector));
  ws->Hii = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->Hij = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->lambda = (t_vector *) malloc (N*sizeof (t_vector));
  ws->buffer = (t_vector *) malloc (N*sizeof (t_vector));
  ws->Dii = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->Dii_inv = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->Lij = (t_matrix *) malloc ((N-1)*sizeof (t_matrix));

  return ws;
}

/* initialize the joint solution workspace */
void fdna_joint_solution_workspace_init (fdna_simcontext *fdna_simcon) {
  (void) fdna_simcon;
}

/* free all the memory associated with the joint solution workspace */
void fdna_joint_solution_workspace_free (fdna_joint_solution_workspace *ws) {
  free (ws->t);
  free (ws->rhs);
  free (ws->Hii);
  free (ws->Hij);
  free (ws->lambda);
  free (ws->buffer);
  free (ws->Dii);
  free (ws->Dii_inv);
  free (ws->Lij);
  free (ws);
}



/* macros to allow for compact writing */
#define WS_T(x) (fdna_simcon->ws->t [x])
#define WS_ERR(x) (fdna_simcon->ws->err [x])
#define WS_JV(x) (fdna_simcon->ws->JV [x])
#define WS_RHS(x) (fdna_simcon->ws->rhs [x])
#define WS_HII(x) (fdna_simcon->ws->Hii [x])
#define WS_HIJ(x) (fdna_simcon->ws->Hij [x])
#define WS_LAMBDA(x) (fdna_simcon->ws->lambda [x])
#define WS_BUFFER(x) (fdna_simcon->ws->buffer [x])
#define WS_DII(x) (fdna_simcon->ws->Dii [x])
#define WS_Dinv(x) (fdna_simcon->ws->Dii_inv [x])
#define WS_LIJ(x) (fdna_simcon->ws->Lij [x])
#define WS_LNJ(x) (fdna_simcon->ws->LNj [x])
#define WS_I(x) (fdna_simcon->ws->I [x])
#define DNA_SEG(x) (fdna_simcon->dna->seg [x])



/* this function calculates the exact solution to the problem of the DNA constraints. It is invoked only
 * if the dSpaceCollide function returned that there are no collisions (simcon->numc==0). */
void fdna_joints_solution (simcontext *simcon, fdna_simcontext *fdna_simcon) {
  int c;
  const int N = fdna_simcon->dna->nsegments;
  t_matrix prev_T, next_T;
  t_matrix LDLt;
  t_vector buff;
  dna_segment *prev_seg, *next_seg;
  const t_real ks = PAR(erp_glob)/(PAR(ode_step)*PAR(ode_step));
  const t_real kd = 1./PAR(ode_step);
  const t_real cfm = PAR(cfm_glob);

  /* cycle on all DNA joints */
  next_seg = DNA_SEG (0);
  segment_next_t (WS_T (0), DNA_SEG (0));
  matrix_vector_cross (next_T, WS_T (0));
  for (c=0; c<N-1; c++) {
    /* initialize the bodies */
    prev_seg = next_seg;
    next_seg = DNA_SEG (c+1);

    /* fetch t_(c+1) */
    segment_next_t (WS_T (c+1), next_seg);

    /* right-hand side of equation */
    segment_rhs (WS_RHS (c), prev_seg, next_seg, WS_T (c), WS_T (c+1), ks, kd);

    /* calculate T matrix */
    matrix_copy (prev_T, next_T);
    matrix_vector_cross (next_T, WS_T (c+1));

    /* calculating H = J * M^-1 * J^T */
    segment_H (WS_HII (c), WS_HIJ (c), prev_seg, next_seg, prev_T, next_T, cfm, kd);
  }

  /* A = L D L^T decomposition */
  matrix_copy (WS_DII (0), WS_HII (0));
  dCopyVector3 (WS_BUFFER (0), WS_RHS (0));

  /* forward substitution */
  for (c=1; c<N-1; c++) {
    matrix_inverse (WS_Dinv (c-1), WS_DII (c-1));
    matrix_product (WS_LIJ (c-1), WS_HIJ (c-1), WS_Dinv (c-1));
    matrix_ABAt (LDLt, WS_LIJ (c-1), WS_DII (c-1));
    matrix_diff (WS_DII (c), WS_HII (c), LDLt);
    matrix_vector_product (buff, WS_LIJ (c-1), WS_BUFFER (c-1));

    /* L Y' = X */
    vector_diff (WS_BUFFER (c), WS_RHS (c), buff);
  }

  /* check_LDL_decomposition (N-1,
      fdna_simcon->ws->Lij,
      fdna_simcon->ws->Dii, 
      fdna_simcon->ws->Hij,
      fdna_simcon->ws->Hii); */

  /* check_forward_substitution (N-1,
      fdna_simcon->ws->Lij,
      fdna_simcon->ws->rhs,
      fdna_simcon->ws->buffer); */

  /* back substitution */
  matrix_inverse (WS_Dinv (N-2), WS_DII (N-2));
  matrix_vector_product (WS_LAMBDA (N-2), WS_Dinv (N-2), WS_BUFFER (N-2));
  for (c=N-3; c>=0; c--) {
    t_vector buff2;
    t_matrix Lt;

    matrix_vector_product (buff, WS_Dinv (c), WS_BUFFER (c));
    matrix_transpose (Lt, WS_LIJ (c));
    matrix_vector_product (buff2, Lt, WS_LAMBDA (c+1));
    vector_diff (WS_LAMBDA (c), buff, buff2);
  }

  /* check_back_substitution (N-1,
      fdna_simcon->ws->Lij,
      fdna_simcon->ws->Dii,
      fdna_simcon->ws->lambda,
      fdna_simcon->ws->buffer); */

  /* check_LDL_solution (N-1,
      fdna_simcon->ws->Hii,
      fdna_simcon->ws->Hij,
      fdna_simcon->ws->lambda,
      fdna_simcon->ws->rhs); */

  /* generalised forces */
  for (c=0; c<N-1; c++) {
    t_vector fi;

    /* add forces and torques */
    body_add_force (DNA_SEG (c)->b_id, WS_LAMBDA (c));
    vector_cross (buff, WS_T (c), WS_LAMBDA (c));
    body_add_torque (DNA_SEG (c)->b_id, buff);

    vector_copy (fi, WS_LAMBDA (c));
    dNegateVector3 (fi);
    body_add_force (DNA_SEG (c+1)->b_id, fi);
    vector_cross (buff, WS_T (c+1), WS_LAMBDA (c));
    body_add_torque (DNA_SEG (c+1)->b_id, buff);
  }
}
