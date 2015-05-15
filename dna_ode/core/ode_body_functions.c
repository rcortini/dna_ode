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

#include "ode_body_functions.h"

/* 
 * ODE-RELATED FUNCTIONS THAT ARE USED FOR CONVENIENCE
 */

/* returns the vector containing the position of the body */
void body_position (t_real *res, dBodyID b) {
  dCopyVector3 (res, dBodyGetPosition (b));
}

/* returns the position of a point according to the local frame: r = x u + y v + z t */
void body_point_position (t_real *res, dBodyID b, t_real x, t_real y, t_real z) {
  dBodyGetRelPointPos (b, x, y, z, res);
}

/* sets the position of the body to vec */
void body_set_position (dBodyID b, const t_real *vec) {
  dBodySetPosition (b, vec[0], vec[1], vec[2]);
}

/* returns the a-th principal axis of a body */
void body_axis (t_real *axis, dBodyID b, int a) {
  if (b == SURFACE_BODY_ID) {
    if (a==0) {
      axis[0] = 0.;
      axis[1] = 0.;
      axis[2] = -1.;
    }
    else if (a==1) {
      axis[0] = 0.;
      axis[1] = 1.;
      axis[2] = 0.;
    }
    else if (a==2) {
      axis[0] = 1.;
      axis[1] = 0.;
      axis[2] = 0.;
    }
  }
  else {
    const t_real *o = dBodyGetRotation (b);
    axis[0] = o[0+a];
    axis[1] = o[4+a];
    axis[2] = o[8+a];
  }
}

/* sets the orientation axes of a body */
void body_set_axes (dBodyID b, const t_real *u, const t_real *v, const t_real *t) {
  dMatrix3 r;
  r [0] = u[0]; r [4] = u[1]; r [8]  = u[2];
  r [1] = v[0]; r [5] = v[1]; r [9]  = v[2];
  r [2] = t[0]; r [6] = t[1]; r [10] = t[2];
  r [3] = 0.;   r [7] = 0.;   r [11] = 0.;
  dBodySetRotation (b, r);
}

/* sets position and rotation of a body based on position and quaternion values */
void body_set_position_and_quaternion (dBodyID b, t_real x, t_real y, t_real z, t_real q1, t_real q2, t_real q3, t_real q4) {
  dQuaternion q;
  dBodySetPosition (b, x, y, z);
  q [0] = q1;
  q [1] = q2;
  q [2] = q3;
  q [3] = q4;
  dBodySetQuaternion (b, q);
}

/* rotates a body by a certain angle around one of its axes */
void body_rotate_around_axis (dBodyID b, t_real angle, unsigned int axis_id) {
  t_vector u, v, t, rotation_axis;
  t_vector new_u, new_v, new_t;

  /* get old axes */
  body_axis (u, b, 0);
  body_axis (v, b, 1);
  body_axis (t, b, 2);

  if (axis_id==0)
    dCopyVector3 (rotation_axis, u);
  else if (axis_id==1)
    dCopyVector3 (rotation_axis, v);
  else if (axis_id==2)
    dCopyVector3 (rotation_axis, t);
  else {
    err_message ("Invalid axis ID\n");
    exit (EXIT_FAILURE);
  }

  /* rotate axes */
  vector_rotate_around_axis (new_u, u, angle, rotation_axis);
  vector_rotate_around_axis (new_v, v, angle, rotation_axis);
  vector_rotate_around_axis (new_t, t, angle, rotation_axis);

  /* set the new axes to the body */
  body_set_axes (b, new_u, new_v, new_t);
}

/* returns the mass of the body */
t_real body_mass (dBodyID b) {
  dMass m_id;
  dBodyGetMass (b, &m_id);
  return m_id.mass;
}


/* returns the principal components of the inertia tensor */
void body_inertia (t_real *res, dBodyID b) {
  dMass m;
  dBodyGetMass (b,&m);
  res [0] = m.I [0];
  res [1] = m.I [5];
  res [2] = m.I [10];
}


/* return the principal inertia tensor of a body (inverse = 0) or its inverse (inverse = 1) */
void body_inertia_tensor (t_real *RIRt, dBodyID b, int inverse) {
  const t_real *R;
  t_matrix mi;
  t_vector vi;

  /* get principal inertia components */
  body_inertia (vi, b);

  if (inverse==0) {
    mi[0] = vi[0]; mi[1] = 0.   ; mi[2]  = 0.   ;
    mi[4] = 0.   ; mi[5] = vi[1]; mi[6]  = 0.   ;
    mi[8] = 0.   ; mi[9] = 0.   ; mi[10] = vi[2];
  }
  else {
    mi[0] = 1./vi[0]; mi[1] = 0.      ; mi[2]  = 0.      ;
    mi[4] = 0.      ; mi[5] = 1./vi[1]; mi[6]  = 0.      ;
    mi[8] = 0.      ; mi[9] = 0.      ; mi[10] = 1./vi[2];
  }

  /* fix extra elements */
  mi[3] = mi[7] = mi[11] = 0.;

  /* perform the matrix multiplication to bring it to the
   * world reference frame */
  R = dBodyGetRotation (b);
  matrix_ABAt (RIRt, R, mi);
}



/* returns the rotation matrix of a body */
void body_rotation (t_real *res, dBodyID b) {
  if (b == SURFACE_BODY_ID) {
    res [0] = 0.;
    res [4] = 0.;
    res [8] = -1.;

    res [1] = 0.;
    res [5] = 1.;
    res [9] = 0.;

    res [2] = 1.;
    res [6] = 0.;
    res [10] = 0.;

    res [3] = res [7] = res [11] = 0.;
  }
  else
    matrix_copy (res, dBodyGetRotation(b));
}



/* returns the velocity of the body b in the vector vec */
void body_velocity (t_real *res, dBodyID b) {
  dCopyVector3 (res, dBodyGetLinearVel (b));
}


/* sets the position of the body to vec */
void body_set_velocity (dBodyID b, const t_real *vec, const t_real *omega) {
  dBodySetLinearVel (b, vec[0], vec[1], vec[2]);
  dBodySetAngularVel (b, omega[0], omega[1], omega[2]);
}


/* returns the angular velocity of body b in the lab frame reference */
void body_angular_velocity (t_real *res, dBodyID b) {
  dCopyVector3 (res, dBodyGetAngularVel (b));
}


/* returns the angular velocity of body b in the body reference frame */
void body_angular_velocity_principal_axes (t_real *res, dBodyID b) {
  const t_real * omega = dBodyGetAngularVel (b);
  dBodyVectorFromWorld (b, omega[0], omega[1], omega[2], res);
}


/* returns a vector that expresses the angular velocity of a point that has coordinates
 * (x, y, z) in the material reference frame (u, v, t) */
void body_point_velocity (t_real *res, dBodyID b, t_real x, t_real y, t_real z) {
  dBodyGetRelPointVel (b, x, y, z, res);
}


/* assigns both velocity and angular velocity to a body */
void body_set_velocity_and_angular_velocity (dBodyID b, t_real vx, t_real vy, t_real vz, t_real omega_x, t_real omega_y, t_real omega_z) {
  dBodySetLinearVel (b, vx, vy, vz);
  dBodySetAngularVel (b, omega_x, omega_y, omega_z);
}


/* returns the value of the torque applied on a body */
void body_torque (t_real *res, dBodyID b) {
  dCopyVector3 (res, dBodyGetTorque (b));
}


/* sets a torque given by the vector vec to the body b */
void body_set_torque (dBodyID b, const t_real *vec) {
  /* let's check if we have the surface, in which case we exit doing nothing */
  if (b == SURFACE_BODY_ID) return;

  /* in all other cases, we add the torque to the body */
  dBodySetTorque (b, vec[0], vec[1], vec[2]);
}


/* adds a torque given by the vector vec to the body b */
void body_add_torque (dBodyID b, const t_real *vec) {
  /* let's check if we have the surface, in which case we exit doing nothing */
  if (b == SURFACE_BODY_ID) return;

  /* in all other cases, we add the torque to the body */
  dBodyAddTorque (b, vec[0], vec[1], vec[2]);
}


/* adds a torque given by the vector vec to the body b,
 * relative to the body reference frame */
void body_add_rel_torque (dBodyID b, const t_real *vec) {
  /* let's check if we have the surface, in which case we exit doing nothing */
  if (b == SURFACE_BODY_ID) return;

  /* in all other cases, we add the torque to the body */
  dBodyAddRelTorque (b, vec[0], vec[1], vec[2]);
}


/* returns the force acting on the body */
void body_force (t_real *res, dBodyID b) {
  dCopyVector3 (res, dBodyGetForce (b));
}


/* adds a force given by the vector vec to the body b */
void body_add_force (dBodyID b, const t_real *vec) {
  dBodyAddForce (b, vec[0], vec[1], vec[2]);
}


/* returns the total kinetic energy of a body */
t_real body_total_kinetic_energy (dBodyID b) {
  const t_real *v, *omega;
  dVector3 wb;
  t_real trans_energy, rot_energy, energy;
  dMass m;

  /* get the linear and angular velocities */
  dBodyGetMass (b, &m);
  v = dBodyGetLinearVel (b);
  omega = dBodyGetAngularVel (b);
  dBodyVectorFromWorld (b, omega [0], omega [1], omega [2], wb);

  /* and assign the energy */
  /* TODO: restore */
  /* energy = 0.5 * m.mass * (v [0] * v [0] + v [1] * v [1] + v [2] * v [2]);
  energy += 0.5 * (m.I [0] * wb[0] * wb[0] + m.I[5] * wb[1] * wb[1] + m.I[10] * wb[2] * wb[2]); */
  trans_energy = 0.5 * m.mass * (v [0] * v [0] + v [1] * v [1] + v [2] * v [2]);
  rot_energy = 0.5 * (m.I [0] * wb[0] * wb[0] + m.I[5] * wb[1] * wb[1] + m.I[10] * wb[2] * wb[2]);
  energy = trans_energy + rot_energy;

  /* printf ("trans: %.3f    rot: %.3f\n", trans_energy, rot_energy); */

  /* and return */
  return energy;
}


/* returns the rotational z kinetic energy of a body */
t_real body_zrot_kinetic_energy (dBodyID b) {
  dMass m;
  const t_real *omega;
  dVector3 wb;

  dBodyGetMass (b, &m);
  omega = dBodyGetAngularVel (b);
  dBodyVectorFromWorld (b, omega [0], omega [1], omega [2], wb);

  return 0.5*m.I[10]*wb[2]*wb[2];
}

/* returns the joint anchoring point */
void joint_anchor (t_real *anchor, dJointID j) {
  dVector3 joint_position;
  int joint_type;
  joint_type = dJointGetType (j);
  if (joint_type == dJointTypeBall) {
    dJointGetBallAnchor (j, joint_position);
  }
  else {
    err_message ("Unrecognized joint type\n");
    exit (EXIT_FAILURE);
  }
  vector_set (anchor, joint_position [0], joint_position [1], joint_position [2]);
}

/* returns the joint second anchoring point */
void joint_anchor2 (t_real *anchor, dJointID j) {
  dVector3 joint_position;
  int joint_type;
  joint_type = dJointGetType (j);
  if (joint_type == dJointTypeBall) {
    dJointGetBallAnchor2 (j, joint_position);
  }
  else {
    err_message ("Unrecognized joint type\n");
    exit (EXIT_FAILURE);
  }
  vector_set (anchor, joint_position [0], joint_position [1], joint_position [2]);
}
