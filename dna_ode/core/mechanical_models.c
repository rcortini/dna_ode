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

#include "mechanical_models.h"

/* calculates the torque and twist between two DNA segments, given their tangent vectors. Returns
 * the vector with the torque, and also assigns the Tw value with the current value of the twist.
 * NONLINEAR MODEL (P. Carrivain's thesis, equation 2.13) */
void body_twist_and_torque_kinkable (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p) {
  t_vector bn, t1_plus_t2;
  t_real one_plus_t1_dot_t2;
  t_real *par = (t_real *) p;

  /* calculate 1 + t1.t2 and twist*/
  one_plus_t1_dot_t2 = 1. + T1_T2 (R1, R2);
  twist_12 (R1, R2, one_plus_t1_dot_t2, Tw);

  /* calculate the binormal to t2 and t1 */
  T1_x_T2 (bn, R1, R2);

  t1_plus_t2_vec (t1_plus_t2, R1, R2);
  vector_lincomb (torque, bn, t1_plus_t2, par [0], par [1] * (*Tw) / one_plus_t1_dot_t2);
}



/* calculates the torque and twist between two DNA segments, given their tangent vectors. Returns
 * the vector with the torque, and also assigns the Tw value with the current value of the twist.
 * LINEAR MODEL. (Derived from a bending energy which is ~gb theta^2. */
void body_twist_and_torque_harmonic (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p) {
  t_vector bn, t1_plus_t2;
  t_real one_plus_t1_dot_t2, theta;
  t_real *par = (t_real *) p;

  /* calculate 1 + t1.t2 */
  one_plus_t1_dot_t2 = 1. + T1_T2 (R1, R2);

  /* calculate the binormal to t2 and t1 */
  T1_x_T2 (bn, R1, R2);
  vector_unit (bn);

  /* angle between the tangent vectors */
  theta = theta_angle (R1, R2);

  /* calculate twist */
  twist_12 (R1, R2, one_plus_t1_dot_t2, Tw);

  t1_plus_t2_vec (t1_plus_t2, R1, R2);
  vector_lincomb (torque, bn, t1_plus_t2, par [0] * theta, par [1] * (*Tw) / one_plus_t1_dot_t2);
}



/* calculates the torque and twist between two DNA segments, given their tangent vectors. Returns
 * the vector with the torque, and also assigns the Tw value with the current value of the twist */
void body_twist_and_torque_nicked_kinkable (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p) {
  t_real *gb = (t_real *) p;

  /* calculate the binormal to t2 and t1 */
  T1_x_T2 (torque, R1, R2);

  /* this is the torque: no component along the twist */
  dScaleVector3 (torque, *gb);

  /* assign zero twist */
  *Tw = 0.;
}



/* calculates the torque and twist between two DNA segments, given their tangent vectors. Returns
 * the vector with the torque, and also assigns the Tw value with the current value of the twist */
void body_twist_and_torque_nicked_harmonic (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p) {
  t_vector bn;
  t_real theta;
  t_real *gb = (t_real *) p;

  /* calculate the binormal to t2 and t1 */
  T1_x_T2 (bn, R1, R2);
  vector_unit (bn);

  /* angle between the tangent vectors */
  theta = theta_angle (R1, R2);

  /* this is the torque: no component along the twist */
  vector_scale (torque, bn, (*gb) * theta);

  /* assign zero twist */
  *Tw = 0.;
}



/* calculates the torque and twist between two DNA segments, given their tangent vectors. Returns
 * the vector with the torque, and also assigns the Tw value with the current value of the twist.
 * HARMONIC4 MODEL. (Derived from a bending energy which is ~gb theta^2 + k4 theta^4. */
void body_twist_and_torque_harmonic4 (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p) {
  t_vector bn, t1_plus_t2;
  t_real one_plus_t1_dot_t2, theta;
  t_real *par = (t_real *) p;

  /* calculate 1 + t1.t2 */
  one_plus_t1_dot_t2 = 1. + T1_T2 (R1, R2);

  /* calculate the binormal to t2 and t1 */
  T1_x_T2 (bn, R1, R2);
  vector_unit (bn);

  /* angle between tangent vectors */
  theta = theta_angle (R1, R2);

  /* calculate the twist */
  twist_12 (R1, R2, one_plus_t1_dot_t2, Tw);

  t1_plus_t2_vec (t1_plus_t2, R1, R2);
  vector_lincomb (torque, bn, t1_plus_t2, par [0] * theta + par [2] * theta*theta*theta, par [1] * (*Tw) / one_plus_t1_dot_t2);
}



/* calculates the torque and twist between two DNA segments, given their tangent vectors. Returns
 * the vector with the torque, and also assigns the Tw value with the current value of the twist */
void body_twist_and_torque_nicked_harmonic4 (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p) {
  t_real theta;
  t_real *par = (t_real *) p;

  /* calculate the binormal to t2 and t1 */
  T1_x_T2 (torque, R1, R2);
  vector_unit (torque);

  /* angle between the tangent vectors */
  theta = theta_angle (R1, R2);

  /* this is the torque: no component along the twist */
  dScaleVector3 (torque, par [0] * theta + par [1] * theta*theta*theta);

  /* assign zero twist */
  *Tw = 0.;
}



/* calculates the torque and twist between two DNA segments, given their tangent vectors. Returns
 * the vector with the torque, and also assigns the Tw value with the current value of the twist.
 * INHOMOGENEOUS LINEAR MODEL. (Derived from a bending energy which is ~gb theta^2. */
void body_twist_and_torque_harmonic_inhom (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p) {
  t_vector bn, t1_plus_t2;
  t_real one_plus_t1_dot_t2, theta;
  inhom_par *par = (inhom_par *) p;

  /* calculate 1 + t1.t2 */
  one_plus_t1_dot_t2 = 1. + T1_T2 (R1, R2);

  /* calculate the binormal to t2 and t1 */
  T1_x_T2 (bn, R1, R2);
  vector_unit (bn);

  /* angle between the tangent vectors */
  theta = theta_angle (R1, R2);

  /* calculate twist */
  twist_12 (R1, R2, one_plus_t1_dot_t2, Tw);

  t1_plus_t2_vec (t1_plus_t2, R1, R2);

  /* get the segment identity */
  unsigned int count = par->count;
  unsigned int njoints = par->njoints;
  t_real gb = par->gb [count];
  t_real gt = par->gt [count];

  /* calculate torque */
  vector_lincomb (torque, bn, t1_plus_t2, gb * theta, gt * (*Tw) / one_plus_t1_dot_t2);

  /* increment counter */
  par->count++;
  if (par->count==njoints)
    par->count = 0;
}
