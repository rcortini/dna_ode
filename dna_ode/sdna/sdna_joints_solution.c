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
#include "sdna.h"


/* alloc all the memory necessary for the joint solution operations */
sdna_joint_solution_workspace * sdna_joint_solution_workspace_alloc (const unsigned int N) {
  sdna_joint_solution_workspace *ws = (sdna_joint_solution_workspace *) malloc (sizeof (*ws));

  ws->t = (t_vector *) malloc ((N+1)*sizeof (t_vector));
  ws->lambda  = (t_vector *) malloc ((N+1)*sizeof (t_vector));
  ws->buffer  = (t_vector *) malloc ((N+1)*sizeof (t_vector));
  ws->rhs = (t_vector *) malloc ((N+1)*sizeof (t_vector));
  ws->Hii = (t_matrix *) malloc ((N+1)*sizeof (t_matrix));
  ws->Hij = (t_matrix *) malloc (N*sizeof (t_matrix));
  ws->Dii = (t_matrix *) malloc ((N+1)*sizeof (t_matrix));
  ws->Dii_inv = (t_matrix *) malloc ((N+1)*sizeof (t_matrix));
  ws->Lij = (t_matrix *) malloc (N*sizeof (t_matrix));

  return ws;
}


/* initialize workspace */
void sdna_joint_solution_workspace_init (sdna_simcontext *sdna_simcon) {
  (void) sdna_simcon;
}


/* free all the memory associated with the joint solution workspace */
void sdna_joint_solution_workspace_free (sdna_joint_solution_workspace *ws) {
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
#define WS_T(x) (sdna_simcon->ws->t [x])
#define WS_RHS(x) (sdna_simcon->ws->rhs [x])
#define WS_HII(x) (sdna_simcon->ws->Hii [x])
#define WS_HIJ(x) (sdna_simcon->ws->Hij [x])
#define WS_LAMBDA(x) (sdna_simcon->ws->lambda [x])
#define WS_BUFFER(x) (sdna_simcon->ws->buffer [x])
#define WS_DII(x) (sdna_simcon->ws->Dii [x])
#define WS_Dinv(x) (sdna_simcon->ws->Dii_inv [x])
#define WS_LIJ(x) (sdna_simcon->ws->Lij [x])
#define DNA_SEG(x) (sdna_simcon->dna->seg [x])



/* this function calculates the exact solution to the problem of the DNA constraints. It is invoked only
 * if the dSpaceCollide function returned that there are no collisions (simcon->numc==0). */
void sdna_joints_solution (simcontext *simcon, sdna_simcontext *sdna_simcon) {
  int c;
  const int N = sdna_simcon->dna->njoints;
  t_matrix prev_T, next_T;
  t_matrix LDLt, tTt;
  t_vector buff, fi;
  dna_segment *prev_seg, *next_seg;
  const t_real ks = PAR(erp_glob) / (PAR(ode_step)*PAR(ode_step));
  const t_real kd = 1./PAR(ode_step);
  const t_real cfm = PAR(cfm_glob);
  const t_real diagonal_factor = 1./PAR (dna_segment_mass) + cfm*kd;

  /* first DNA segment has a different matrix element */
  next_seg = DNA_SEG (0);
  segment_next_t (WS_T (0), DNA_SEG (0));
  segment_rhs_ground (WS_RHS (0), next_seg, WS_T (0), ks, kd);

  /* calculate H00 term, which is different from the others */
  matrix_vector_cross (next_T, WS_T (0));
  matrix_negate (next_T);
  matrix_ABAt (tTt, next_T, next_seg->inertia);
  matrix_copy (WS_HII (0), tTt);
  WS_HII (0)[0] += diagonal_factor;
  WS_HII (0)[5] += diagonal_factor;
  WS_HII (0)[10] += diagonal_factor;
  matrix_copy (WS_HIJ (0), tTt);
  WS_HIJ (0)[0] -= 1./sdna_simcon->bead->mass;
  WS_HIJ (0)[5] -= 1./sdna_simcon->bead->mass;
  WS_HIJ (0)[10] -= 1./sdna_simcon->bead->mass;

  /* cycle on all DNA segments */
  for (c=1; c<N-1; c++) {
    /* initialize the bodies */
    prev_seg = next_seg;
    next_seg = DNA_SEG (c);

    /* fetch t_(c) */
    segment_next_t (WS_T (c), next_seg);

    /* right-hand side of equation */
    segment_rhs (WS_RHS (c), prev_seg, next_seg, WS_T (c-1), WS_T (c), -ks, -kd);

    /* calculate T matrix */
    matrix_copy (prev_T, next_T);
    matrix_vector_cross (next_T, WS_T (c));

    /* calculating H = J * M^-1 * J^T */
    segment_H (WS_HII (c), WS_HIJ (c), prev_seg, next_seg, prev_T, next_T, cfm, kd);
  }

  /* rhs, Hii, and Hij of the bead */
  bead_joint (WS_T (N-1), sdna_simcon->bead);
  segment_rhs_bead (WS_RHS (N-1), next_seg, sdna_simcon->bead, WS_T (N-2), WS_T (N-1), -ks, -kd);
  matrix_copy (prev_T, next_T);
  matrix_vector_cross (next_T, WS_T (N-1));
  bead_H (WS_HII (N-1), next_seg, sdna_simcon->bead, prev_T, next_T, cfm, kd);

  /* A = L D L^T decomposition */
  matrix_copy (WS_DII (0), WS_HII (0));
  dCopyVector3 (WS_BUFFER (0), WS_RHS (0));

  /* forward substitution */
  for (c=1; c<N; c++) {
    matrix_inverse (WS_Dinv (c-1), WS_DII (c-1));
    matrix_product (WS_LIJ (c-1), WS_HIJ (c-1), WS_Dinv (c-1));
    matrix_ABAt (LDLt, WS_LIJ (c-1), WS_DII (c-1));
    matrix_diff (WS_DII (c), WS_HII (c), LDLt);
    matrix_vector_product (buff, WS_LIJ (c-1), WS_BUFFER (c-1));

    /* L Y' = X */
    vector_diff (WS_BUFFER (c), WS_RHS (c), buff);
  }

  /* check_LDL_decomposition (N,
      sdna_simcon->ws->Lij,
      sdna_simcon->ws->Dii, 
      sdna_simcon->ws->Hij,
      sdna_simcon->ws->Hii); */

  /* check_forward_substitution (N,
      sdna_simcon->ws->Lij,
      sdna_simcon->ws->rhs,
      sdna_simcon->ws->buffer); */

  /* back substitution */
  matrix_inverse (WS_Dinv (N-1), WS_DII (N-1));
  matrix_vector_product (WS_LAMBDA (N-1), WS_Dinv (N-1), WS_BUFFER (N-1));
  for (c=N-2; c>=0; c--) {
    t_vector buff2;
    t_matrix Lt;

    matrix_vector_product (buff, WS_Dinv (c), WS_BUFFER (c));
    matrix_transpose (Lt, WS_LIJ (c));
    matrix_vector_product (buff2, Lt, WS_LAMBDA (c+1));
    vector_diff (WS_LAMBDA (c), buff, buff2);
  }

  /* check_back_substitution (N,
      sdna_simcon->ws->Lij,
      sdna_simcon->ws->Dii,
      sdna_simcon->ws->lambda,
      sdna_simcon->ws->buffer); */

  /* check_LDL_solution (N,
      sdna_simcon->ws->Hii,
      sdna_simcon->ws->Hij,
      sdna_simcon->ws->lambda,
      sdna_simcon->ws->rhs); */

  /* generalised forces */

  /* first segment */
  body_add_force (DNA_SEG (0)->b_id, WS_LAMBDA (0));
  vector_copy (fi, WS_LAMBDA (0));
  dNegateVector3 (fi);
  vector_cross (buff, WS_T (0), fi);
  body_add_torque (DNA_SEG (0)->b_id, buff);

  /* central DNA segments */
  for (c=0; c<(N-2); c++) {
    vector_copy (fi, WS_LAMBDA (c+1));

    body_add_force (DNA_SEG (c+1)->b_id, fi);

    dNegateVector3 (fi);

    vector_cross (buff, WS_T (c), fi);
    body_add_torque (DNA_SEG (c)->b_id, buff);

    vector_cross (buff, WS_T (c+1), fi);
    body_add_torque (DNA_SEG (c+1)->b_id, buff);

    body_add_force (DNA_SEG (c)->b_id, fi);
  }

  /* bead */
  vector_copy (fi, WS_LAMBDA (c+1));

  body_add_force (sdna_simcon->bead->b_id, fi);

  dNegateVector3 (fi);
  body_add_force (DNA_SEG (c)->b_id, fi);

  vector_cross (buff, WS_T (c), fi);
  body_add_torque (DNA_SEG (c)->b_id, buff);

  vector_cross (buff, WS_T (c+1), fi);
  body_add_torque (sdna_simcon->bead->b_id, buff);
}
