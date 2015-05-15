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


#include "dna.h"


/* this is passed to the root finder to extract the value of gb */
t_real gb_function (t_real x, void *p) {
  t_real *mcos = (t_real *) p;
  return langevin_function (x) - *mcos;
}



/* calculates the mechanical constants gb and gt */
void dna_gb_gt (t_real kuhn_length, t_real segment_length, t_real ltp, t_real tol, t_real *gb, t_real *gt) {
  t_real k_s, mcos;

  /* calculate cos (m) */
  k_s = kuhn_length/segment_length;
  mcos = (k_s-1.)/(k_s+1.);

  /* now we calculate gb. If there is only one segment per Kuhn length, there is no bending rigidity, and gb and gt are zero */
  if (kuhn_length/segment_length < 2.) {
    *gb=0.;
    *gt=0.;
  }
  else {
    int f_root_exit_status;
    struct f_root_params f_root_p;
    f_root_p.f = &gb_function;
    f_root_p.p = &mcos;
    f_root_p.type = gsl_root_fsolver_brent;
    f_root_p.eps_rel = 0.;
    f_root_p.eps_abs = tol;
    f_root_p.max_iter = 1000;

    /* we calculate iteratively the value of gb, using a "brent" root finding method */
    f_root_exit_status = f_root (3./8.*kuhn_length/segment_length, kuhn_length/segment_length, gb, &f_root_p);

    if (f_root_exit_status != GSL_SUCCESS) {
      err_message ("Error: impossible to find gb\n");
      exit (EXIT_FAILURE);
    }

    /* this is gt */
    *gt = ltp/segment_length;
  }
}



/* calculates DNA total twist, and adds mechanical torques to every segment.
 * Note that if there is a bead, it will be taken into account, because
 * of the cycle on the joints. The last joint will point to the bead in next_body */
t_real dna_twist_and_torque (dna_ode *dna, mechanical_model_f *mechanical_model) {
  unsigned int i;
  const unsigned int J = dna->njoints;
  t_real twist, total_twist = 0.;
  t_vector torque;

  /* cycle on all DNA joints */
  for (i=0; i<J; i++) {
    t_matrix prev_R, next_R;
    joint_ode *joint = dna->joint [i];
    dBodyID prev_body = *joint->previous_body;
    dBodyID next_body = *joint->next_body;

    /* enable the joint */
    dJointEnable (joint->j_id);

    /* get rotation matrices */
    body_rotation (prev_R, prev_body);
    body_rotation (next_R, next_body);

    /* calculate the torque acting between the two bodies */
    mechanical_model->f (torque, prev_R, next_R, &twist, mechanical_model->p);

    /* printf ("torque [%d] = ", i);
       print_vector3d (torque);
       printf ("\n"); */

    /* add torque to the segment, and also to the previous segment for mechanical equilibrium */
    body_add_torque (prev_body, torque);
    dNegateVector3 (torque);
    body_add_torque (next_body, torque);

    /* print a warning if the twist changed a lot during this time step */
    if ((!dna->nicked) && fabs (dna->joint [i]->twist - twist) > M_PI)
      err_message ("Warning: joint %d previous twist = %f actual twist = %f\n", i, joint->twist, twist);

    /* set the value of the twist of the current segment in the twist table and increment the total twist variation */
    joint->twist = twist;
    total_twist += twist;
  }
  return total_twist;
}



/* calculates DNA total twist */
t_real dna_twist (dna_ode *dna) {
  unsigned int i;
  const unsigned int J = dna->njoints;
  t_real twist, total_twist = 0.;

  /* cycle on all DNA joints */
  for (i=0; i<J; i++) {
    joint_ode *joint = dna->joint [i];
    t_matrix prev_R, next_R; 

    body_rotation (prev_R, *joint->previous_body);
    body_rotation (next_R, *joint->next_body);

    /* calculate the torque acting between the two bodies */
    twist = body_twist (prev_R, next_R);

    /* set the value of the twist of the current segment in the twist table and increment the total twist variation */
    joint->twist = twist;
    total_twist += twist;
  }
  return total_twist;
}



/* takes a DNA and initializes its velocities and forces */
void dna_shake (dna_ode *dna, sim_body *body_list, unsigned int nstart, mechanical_model_f *m) {
  /* cycle on all the dna segments and perform a Langevin-Euler step */
  for (unsigned int i=0; i<dna->njoints; i++) {
    t_vector this_segment_torque;
    t_matrix prev_R, next_R;
    t_real this_segment_twist;
    joint_ode *joint = dna->joint [i];
    const dBodyID prev_body = *joint->previous_body;
    const dBodyID next_body = *joint->next_body;

    /* perform a single Langevin local thermostat step on the DNA segment */
    langevin_euler_step (&body_list [nstart+i]);

    /* now update the current body position and axes */
    body_rotation (prev_R, prev_body);
    body_rotation (next_R, next_body);

    /* calculate the torque acting between the two bodies */
    m->f (this_segment_torque, prev_R, next_R, &this_segment_twist, m->p);

    /* add torque to the segment */
    body_add_torque (prev_body, this_segment_torque);
    dNegateVector3 (this_segment_torque);
    body_add_torque (next_body, this_segment_torque);

    /* set the value of the twist of the current segment in the twist table */
    joint->twist = this_segment_twist;
  }
}



/* Function to place the DNA segments along a curve */
void place_dna_segments_on_curve (unsigned int nsegments, dna_segment **seg, curve_f *f, t_real total_twist) {
  unsigned int i;
  t_vector z_prime, rs0, r0, r1, r2, t0, u0, v0, u, v, t, rotation_axis, curve_tangent, r_DNA;
  t_real phi, theta, s1, s2, d = 1., rotation_angle; /* TODO generalize d */
  curve_f *rot_trans;
  struct rot_trans_curve_parameters rot_trans_curve_p;

  /* initialize position and orientation of the first segment */
  CURVE_EVAL (rs0, 0, f);
  body_point_position (r_DNA, seg [0]->b_id, 0., 0., d/2.);
  vector_diff (r0, r_DNA, rs0);
  body_axis (u, seg [0]->b_id, 0);
  body_axis (v, seg [0]->b_id, 1);
  body_axis (t, seg [0]->b_id, 2);

  /* define the translated curve */
  curve_t_eval (curve_tangent, 0, f);
  vector_cross (rotation_axis, curve_tangent, t);
  dNormalize3 (rotation_axis);
  rotation_angle = angle_between_vectors (curve_tangent, t);
  rot_trans_curve_p.f = f;
  dCopyVector3 (rot_trans_curve_p.R0, r0);
  dCopyVector3 (rot_trans_curve_p.rotation_axis, rotation_axis);
  dCopyVector3 (rot_trans_curve_p.rotation_origin, rs0);
  rot_trans_curve_p.rotation_angle = rotation_angle;
  rot_trans = curve_init (rot_trans_curve, &rot_trans_curve_p);

  /* initialize */
  s2 = 0.;
  phi = total_twist/(t_real) nsegments;
  for (i=1; i<nsegments; i++) {
    t_vector buff; 

    /* init */
    s1 = s2;
    dCopyVector3 (u0, u);
    dCopyVector3 (v0, v);
    dCopyVector3 (t0, t);

    /* calculate at which value of s2 we have that |r (s1) - r (s2)| = d */
    s2 = s2_s1_d (s1, d, rot_trans);

    /* r1 and r2 are the positions of the joints */
    CURVE_EVAL (r1, s1, rot_trans);
    CURVE_EVAL (r2, s2, rot_trans);

    /* t is the tangent vector, defined as the normalized
     * vector connecting r1 and r2 */
    vector_diff (t, r2, r1);
    dNormalize3 (t);

    /* now we calculate u: we take u0, rotate it a first time by phi along z, and a second time by
     * theta along the axis orthogonal to t and z */
    theta = angle_between_vectors (t, t0);
    vector_cross (z_prime, t0, t);
    dNormalize3 (z_prime);
    vector_rotate_around_axis (buff, t0, phi, u0);
    vector_rotate_around_axis (u, z_prime, theta, buff);

    /* finally, we set u */
    vector_cross (v, t, u);

    /* r0 is the position of the center of the DNA segment */
    vector_sum (r0, r1, r2);
    dScaleVector3 (r0, 0.5);

    /* place the segment */
    body_set_position (seg [i]->b_id, r0);
    body_set_axes (seg [i]->b_id, u, v, t);
  }
  curve_free (rot_trans);
}



/* Creates a loop in the DNA molecule. This function takes nsegments segments
 * and makes a loop out of it. Outside this function one has to take care of what 
 * happens to the molecule extremities. */
void dna_loop_create (unsigned int nsegments, t_real loop_radius, t_real loop_pitch, t_real total_twist, dna_segment **seg) {
  curve_f *loop;
  struct helix_parameters helix_par;

  /* define helix parameters */
  helix_par.radius = loop_radius; 
  helix_par.pitch = loop_pitch;
  loop = curve_init (helix, &helix_par);
  curve_set_t (loop, helix_t);

  /* go */
  place_dna_segments_on_curve (nsegments, seg, loop, total_twist);

  /* free curve */
  curve_free (loop);
}


/* creates a WLC section with given force and force axis */
void dna_wlc_create (unsigned int nsegments, t_real beta, t_real force, t_real gb, const t_real *wlc_axis, dna_segment **seg) {
  unsigned int i;
  t_vector previous_position, next_position, joint_position, buff;
  t_vector previous_u, previous_v, previous_t, next_u, next_v, next_t;
  t_real phi, theta, E_final, E_initial;

  /* this cycle starts from the segment with i=0. This means that the segment with i=0,
   * the "first" segment, will not be touched */
  for (i=1; i<nsegments; i++) {
    /* get the previous segment's position and orientation */
    body_position (previous_position, seg [i-1]->b_id);
    body_axis (previous_u, seg [i-1]->b_id, 0);
    body_axis (previous_v, seg [i-1]->b_id, 1);
    body_axis (previous_t, seg [i-1]->b_id, 2);

    /* generate two random angles: phi = [0, 2pi] and theta = [0, pi] */
    do {
      dCopyVector3 (next_t, previous_t);
      phi = ran_uniform (0., 2*M_PI);
      theta = ran_uniform (0., M_PI);
      vector_rotate_around_axis (buff, previous_u, theta, previous_t);
      vector_rotate_around_axis (next_t, previous_t, phi, buff);
      E_initial = -force * vector_dot (previous_t, wlc_axis);
      E_final = gb * (1.-vector_dot (next_t, previous_t)) - force * vector_dot (next_t, wlc_axis);
    } while (!metropolis (beta, E_initial, E_final));

    /* now rotate the axes of the previous segment by an angle given by the angles */
    vector_rotate_around_axis (next_u, previous_t, phi, previous_u);
    vector_rotate_around_axis (buff, previous_u, theta, previous_v);
    vector_rotate_around_axis (next_v, previous_t, phi, buff);

    /* now assign the new position of the segment */
    segment_next_joint (joint_position, seg [i-1]);
    vector_lincomb (next_position, joint_position, next_t, 1., .5*seg [i]->length);
    body_set_position (seg [i]->b_id, next_position);

    /* and finally its orientation */
    body_set_axes (seg [i]->b_id, next_u, next_v, next_t);
  }
}



/* gets the position of the next joint */
void segment_next_joint (t_real *res, dna_segment *seg) {
  const t_real *pos = dBodyGetPosition (seg->b_id);
  t_vector axis;

  /* get scaled t vector */
  segment_next_t (axis, seg);

  /* sum position and scaled t vector */
  vector_sum (res, pos, axis);
}



/* gets the scaled tangent vector of a DNA segment */
void segment_next_t (t_real *res, dna_segment *seg) {
  body_axis (res, seg->b_id, 2);
  dScaleVector3 (res, .5*seg->length);
}



/* gets the bead joint */
void bead_joint (t_real *res, bead_ode *bead) {
  body_axis (res, bead->b_id, 2);
  dScaleVector3 (res, bead->radius);
}

/* returns the vector res = previous_anchor - next_anchor */
void segment_joint_pos_error (t_real *res, dna_segment *previous, dna_segment *next) {
  dVector3 previous_anchor, next_anchor;
  dBodyGetRelPointPos (previous->b_id, 0., 0., .5*previous->length, previous_anchor);
  dBodyGetRelPointPos (next->b_id, 0., 0., -.5*next->length, next_anchor);
  dSubtractVectors3 (res, previous_anchor, next_anchor);
}

/* returns the vector res = previous_vel - next_vel,
 * where these are the velocities of the anchor points */
void segment_joint_vel_error (t_real *res, dna_segment *previous, dna_segment *next) {
  dVector3 previous_vel, next_vel;
  dBodyGetRelPointVel (previous->b_id, 0., 0., .5*previous->length, previous_vel);
  dBodyGetRelPointVel (next->b_id, 0., 0., -.5*next->length, next_vel);
  dSubtractVectors3 (res, previous_vel, next_vel);
}

/* returns the vector res = previous_anchor - next_anchor  for the bead joint*/
void bead_joint_pos_error (t_real *res, dna_segment *last, bead_ode *bead) {
  dVector3 previous_anchor, next_anchor;
  dBodyGetRelPointPos (last->b_id, 0., 0., .5*last->length, previous_anchor);
  dBodyGetRelPointPos (bead->b_id, 0., 0., -bead->radius, next_anchor);
  dSubtractVectors3 (res, previous_anchor, next_anchor);
}

/* returns the vector res = previous_vel - next_vel,
 * where these are the velocities of the anchor points  for the bead joint*/
void bead_joint_vel_error (t_real *res, dna_segment *last, bead_ode *bead) {
  dVector3 previous_vel, next_vel;
  dBodyGetRelPointVel (last->b_id, 0., 0., .5*last->length, previous_vel);
  dBodyGetRelPointVel (bead->b_id, 0., 0., -bead->radius, next_vel);
  dSubtractVectors3 (res, previous_vel, next_vel);
}


/* calculates the right-hand side term for an internal joint of a DNA, for
 * the exact solution of the joint forces problem */
void segment_rhs (t_real *rhs,
    dna_segment *prev_seg,
    dna_segment *next_seg,
    const t_real *prev_t,
    const t_real *next_t,
    const t_real ks,
    const t_real kd) {

  dBodyID prev_body = prev_seg->b_id;
  dBodyID next_body = next_seg->b_id;
  t_vector err, jv, prev_force, prev_torque, next_force, next_torque, buff;

  /* get scaled forces */
  body_force (prev_force, prev_body);
  dScaleVector3 (prev_force, 1./prev_seg->mass);
  body_force (next_force, next_body);
  dScaleVector3 (next_force, 1./next_seg->mass);

  /* get scaled torques */
  matrix_vector_product (buff, prev_seg->inertia, dBodyGetTorque (prev_body));
  vector_cross (prev_torque, prev_t, buff);
  body_torque (next_torque, next_body);
  matrix_vector_product (buff, next_seg->inertia, dBodyGetTorque (next_body));
  vector_cross (next_torque, next_t, buff);

  /* error and Jacobian */
  segment_joint_pos_error (err, prev_seg, next_seg);
  segment_joint_vel_error (jv, prev_seg, next_seg);
  dScaleVector3 (err, -ks);
  dScaleVector3 (jv, -kd);

  /* calculate rhs */
  vector_sum (rhs, err, jv);
  vector_inplace_sum (rhs, prev_force);
  vector_inplace_diff (rhs, prev_torque);
  vector_inplace_diff (rhs, next_force);
  vector_inplace_diff (rhs, next_torque);
}



/* calculates the right-hand side of the equation for a segment that 
 * is anchored to the ground */
void segment_rhs_ground (t_real *rhs,
    dna_segment *seg,
    const t_real *t,
    const t_real ks,
    const t_real kd) {
  dBodyID body = seg->b_id;
  t_vector err, jv, force, torque, buff;

  /* get scaled force */
  body_force (force, body);
  dScaleVector3 (force, -1./seg->mass);

  /* get scaled torques */
  matrix_vector_product (buff, seg->inertia, dBodyGetTorque (body));
  vector_cross (torque, t, buff);

  /* error and Jacobian */
  body_point_position (err, body, 0., 0., -.5*seg->length);
  body_point_velocity (jv, body, 0., 0., -.5*seg->length);
  dScaleVector3 (err, -ks);
  dScaleVector3 (jv, -kd);

  /* calculate rhs */
  vector_sum (rhs, err, jv);
  vector_inplace_sum (rhs, force);
  vector_inplace_diff (rhs, torque);
}


void segment_rhs_bead (t_real *rhs,
    dna_segment *last_seg,
    bead_ode *bead,
    const t_real *last_t,
    const t_real *bead_t,
    const t_real ks,
    const t_real kd) {

  dBodyID prev_body = last_seg->b_id;
  dBodyID next_body = bead->b_id;
  t_vector err, jv, prev_force, prev_torque, next_force, next_torque, buff;

  /* get scaled forces */
  body_force (prev_force, prev_body);
  dScaleVector3 (prev_force, 1./last_seg->mass);
  body_force (next_force, next_body);
  dScaleVector3 (next_force, 1./bead->mass);

  /* get scaled torques */
  matrix_vector_product (buff, last_seg->inertia, dBodyGetTorque (prev_body));
  vector_cross (prev_torque, last_t, buff);
  body_torque (next_torque, next_body);
  matrix_vector_product (buff, bead->inertia, dBodyGetTorque (next_body));
  vector_cross (next_torque, bead_t, buff);

  /* error and Jacobian */
  bead_joint_pos_error (err, last_seg, bead);
  bead_joint_vel_error (jv, last_seg, bead);
  dScaleVector3 (err, -ks);
  dScaleVector3 (jv, -kd);

  /* calculate rhs */
  vector_sum (rhs, err, jv);
  vector_inplace_sum (rhs, prev_force);
  vector_inplace_diff (rhs, prev_torque);
  vector_inplace_diff (rhs, next_force);
  vector_inplace_diff (rhs, next_torque);

}

/* calculates the H matrix for the exact solution of the joint forces problem */
void segment_H (t_real *Hii,
    t_real *Hij,
    dna_segment *prev_seg,
    dna_segment *next_seg,
    t_real *prev_T,
    t_real *next_T,
    const t_real cfm,
    const t_real kd) {

  const t_real diagonal_factor = 1./prev_seg->mass + 1./next_seg->mass + cfm*kd;
  t_matrix prev_tTt, next_tTt;

  /* calculate the t x I x t^T matrices */
  matrix_ABAt (prev_tTt, prev_T, prev_seg->inertia);
  matrix_ABAt (next_tTt, next_T, next_seg->inertia);

  /* Hii */
  matrix_sum (Hii, prev_tTt, next_tTt);
  Hii[0] += diagonal_factor;
  Hii[5] += diagonal_factor;
  Hii[10] += diagonal_factor;

  /* Hij = Hji term */
  matrix_copy (Hij, next_tTt);
  Hij[0] -= 1./next_seg->mass;
  Hij[5] -= 1./next_seg->mass;
  Hij[10] -= 1./next_seg->mass;
}



/* calculates the H matrix for the exact solution of the joint forces problem */
void bead_H (t_real *Hii,
    dna_segment *last_seg,
    bead_ode *bead,
    t_real *prev_T,
    t_real *next_T,
    const t_real cfm,
    const t_real kd) {

  const t_real diagonal_factor = 1./last_seg->mass + 1./bead->mass + cfm*kd;
  t_matrix prev_tTt, next_tTt;

  /* calculate the t x I x t^T matrices */
  matrix_ABAt (prev_tTt, prev_T, last_seg->inertia);
  matrix_ABAt (next_tTt, next_T, bead->inertia);

  /* Hii */
  matrix_sum (Hii, prev_tTt, next_tTt);
  Hii[0] += diagonal_factor;
  Hii[5] += diagonal_factor;
  Hii[10] += diagonal_factor;
}
