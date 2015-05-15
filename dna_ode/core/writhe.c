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

#include "writhe.h"

t_real writhe_subchain_Klenin (t_vector *r, const t_real L, const unsigned int N) {
  unsigned int i, j;
  t_vector r1, r2, r3, r4;
  t_vector r13, r14, r24, r23, r12, r34;
  t_vector r13_r14, r14_r24, r24_r23, r23_r13, r34_r12;
  t_vector n1, n2, n3, n4;
  t_real r13_r14_norm, r14_r24_norm, r24_r23_norm, r23_r13_norm;
  t_real a1, a2, a3, a4;
  t_real Wr = 0.;

  for (i=1; i<N; i++) {
    for (j=0; j<i; j++) {
      /* r1, r2, r3, r4 */
      dCopyVector3 (r1, r [i]);
      dCopyVector3 (r2, r [i+1]);
      dCopyVector3 (r3, r [j]);
      dCopyVector3 (r4, r [j+1]);

      /* r13, r14, r24, r23 */
      vector_diff (r13, r1, r3);
      vector_diff (r14, r1, r4);
      vector_diff (r24, r2, r4);
      vector_diff (r23, r2, r3);
      vector_diff (r12, r1, r2);
      vector_diff (r34, r3, r4);

      /* r13_r14, r14_r24, r24_r23, r23_r13 */
      vector_cross (r13_r14, r13, r14);
      vector_cross (r14_r24, r14, r24);
      vector_cross (r24_r23, r24, r23);
      vector_cross (r23_r13, r23, r13);
      vector_cross (r34_r12, r34, r12);

      /* norm of r13_r14, r14_r24, r24_r23, r23_r13 */
      r13_r14_norm = vector_norm (r13_r14);
      r14_r24_norm = vector_norm (r14_r24);
      r24_r23_norm = vector_norm (r24_r23);
      r23_r13_norm = vector_norm (r23_r13);

      /* n1, n2, n3, n4 */
      vector_scale (n1, r13_r14, r13_r14_norm);
      vector_scale (n2, r14_r24, r14_r24_norm);
      vector_scale (n3, r24_r23, r24_r23_norm);
      vector_scale (n4, r23_r13, r23_r13_norm);

      /* print_vector3d (n1);
	 printf ("\n");
	 print_vector3d (n2);
	 printf ("\n");
	 print_vector3d (n3);
	 printf ("\n");
	 print_vector3d (n4);
	 printf ("\n"); */

      a1 = safe_asin (vector_dot (n1, n2));
      a2 = safe_asin (vector_dot (n2, n3));
      a3 = safe_asin (vector_dot (n3, n4));
      a4 = safe_asin (vector_dot (n4, n1));

      /* printf ("a1 = %f\n", a1);
	 printf ("a2 = %f\n", a2);
	 printf ("a3 = %f\n", a3);
	 printf ("a4 = %f\n", a4); */

      Wr += (a1 + a2 + a3 + a4) * signum (vector_dot (r34_r12, r13));
    }
  }
  Wr *= L*L/(2*M_PI);

  return Wr;
}



/* calculates the writhe of a DNA subchain. This function may be used even to 
 * calculate the writhe of a whole DNA molecule, but care needs to be taken
 * to discriminate the cases in which the DNA molecule is circular or not
 * (see the function writhe_dna) */
t_real writhe_subchain_Braun (t_vector *r, const t_real L, const unsigned int N) {
  unsigned int i, j;
  t_real Wr = 0.;
  t_vector e_ij, r_ij, t_i, t_j, ti_eij;
  t_real norm_r_ij, norm2_r_ij;

  for (i=1; i<N-2; i++) {
    for (j=i+1; j<N-1; j++) {
      /* r_ij */
      vector_diff (r_ij, r [j], r [i]);
      norm2_r_ij = dCalcVectorLengthSquare3 (r_ij);
      norm_r_ij = sqrt (norm2_r_ij);

      /* e_ij */
      vector_scale (e_ij, r_ij, 1./norm_r_ij);

      /* t_i */
      vector_diff (t_i, r [i+1], r [i-1]);
      dNormalize3 (t_i);

      /* t_j */
      vector_diff (t_j, r [j+1], r [j-1]);
      dNormalize3 (t_j);

      /* writhe component */
      vector_cross (ti_eij, t_i, e_ij);
      Wr += vector_dot (ti_eij, t_j)/norm2_r_ij;
    }
  }
  Wr *= L*L/(2*M_PI);

  return Wr;
}



/* calculates the writhe of an entire DNA molecule */
t_real writhe_dna (dna_ode *dna, const unsigned int is_circular) {
  const unsigned int N = dna->nsegments;
  unsigned int i;
  t_real Wr, L = 1.;

  if (is_circular) {
    t_vector *r = malloc ((N+2)*sizeof (t_vector)); /* TODO: put this elsewhere so we do only a single alloc */

    /* Define the r_i vectors: circular DNA case.
     * Here we define the chain as having two more nodes. The first one is the last
     * joint position, the last one is the first joint position. This way we can invoke
     * the "writhe_subchain" function with N+2 segments, and obtain the correct result. */
    body_point_position (r [0], dna->seg [N-1]->b_id, 0., 0., -0.5*L);
    for (i=1; i<(N+1); i++) {
      body_point_position (r [i], dna->seg [i-1]->b_id, 0., 0., -0.5*L);
    }
    body_point_position (r [N+1], dna->seg [0]->b_id, 0., 0., -0.5*L);

    /* calculate writhe */
    Wr = writhe_subchain_Braun (r, L, N+2);
    free (r);
  }
  else {
    t_vector *r = malloc (N*sizeof (t_vector)); /* TODO: put this elsewhere so we do only a single alloc */

    /* Define the r_i vectors: linear DNA case. We need to add a special case
     * to the last segment, as its end follows a different rule. */
    for (i=0; i<N-1; i++) {
      body_point_position (r [i], dna->seg [i]->b_id, 0., 0., -0.5*L);
    }
    body_point_position (r [N-1], dna->seg [i]->b_id, 0., 0., 0.5*L);

    /* calculate writhe */
    Wr = writhe_subchain_Braun (r, L, N);
    free (r);
  }

  return Wr;
}
