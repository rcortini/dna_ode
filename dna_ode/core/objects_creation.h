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

#ifndef __OBJECTS_CREATION_H__
#define __OBJECTS_CREATION_H__

#include "dna_ode_core.h"
#include "bead.h"
#include "ode_body_functions.h"

surface_ode * surface_create (dSpaceID space_id, t_real a, t_real b, t_real c, t_real d);

bead_ode * bead_create (unsigned int odepdb_id,
    FILE *odepdb_file,
    dWorldID world_id,
    dSpaceID space_id,
    int finite_rotation,
    int gyroscopic_mode,
    t_real mass,
    t_real radius,
    t_real eff_radius,
    t_real *position);

dna_segment * dna_segment_create (unsigned int odepdb_id,
    FILE *odepdb_file,
    dWorldID world_id,
    dSpaceID space_id,
    t_real mass,
    t_real radius,
    t_real length,
    unsigned int axis,
    int finite_rotation,
    int gyroscopic_mode,
    t_real *center,
    t_real *u,
    t_real *v,
    t_real *t);

void joint_create (unsigned int id,
    dWorldID world_id,
    dBodyID *prev,
    dBodyID *next,
    joint_ode *joint,
    const t_real *position,
    unsigned int assign_feedback);

#endif
