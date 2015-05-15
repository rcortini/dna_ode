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

#ifndef __ODE_BODY_FUNCTIONS_H__
#define __ODE_BODY_FUNCTIONS_H__

#include "dna_ode_core.h"
#include "dna_ode_math.h"

void body_position (t_real *res, dBodyID b);

void body_point_position (t_real *res, dBodyID b, t_real x, t_real y, t_real z);

void body_set_position (dBodyID b, const t_real *v);

void body_axis (t_real *res, dBodyID b, int a);

void body_set_axes (dBodyID b, const t_real *u, const t_real *v, const t_real *t);

void body_set_position_and_quaternion (dBodyID b, t_real x, t_real y, t_real z, t_real q1, t_real q2, t_real q3, t_real q4);

void body_rotate_around_axis (dBodyID b, t_real angle, unsigned int axis_id);

void body_set_velocity_and_angular_velocity (dBodyID b, t_real vx, t_real vy, t_real vz, t_real omega_x, t_real omega_y, t_real omega_z);

t_real body_mass (dBodyID b);

void body_inertia (t_real *res, dBodyID b);

void body_inertia_tensor (t_real *res, dBodyID b, int inverse);

void body_rotation (t_real *res, dBodyID b);

void body_velocity (t_real *res, dBodyID b);

void body_set_velocity (dBodyID b, const t_real *vec, const t_real *omega);

void body_angular_velocity (t_real *res, dBodyID b);

void body_angular_velocity_principal_axes (t_real *res, dBodyID b);

void body_point_velocity (t_real *res, dBodyID b, t_real x, t_real y, t_real z);

void body_torque (t_real *res, dBodyID b);

void body_set_torque (dBodyID b, const t_real *vec);

void body_add_torque (dBodyID b, const t_real *vec);

void body_add_rel_torque (dBodyID b, const t_real *vec);

void body_force (t_real *res, dBodyID b);

void body_add_force (dBodyID b, const t_real *vec);

t_real body_total_kinetic_energy (dBodyID b);

t_real body_zrot_kinetic_energy (dBodyID b);

void joint_anchor (t_real *res, dJointID j);

void joint_anchor2 (t_real *res, dJointID j);

#endif
