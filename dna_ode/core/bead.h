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

#ifndef __BEAD_H__
#define __BEAD_H__

#include "dna_ode_core.h"
#include "config_read.h"
#include "ode_body_functions.h"
#include <libconfig.h>

/* possible bead modes */
typedef enum {
  FIXED_N,
  FIXED_TORQUE,
  FIXED_EXT}
bead_mode;

/* a structure to contain the information on the bead */
struct bead_ode {
  bead_mode mode;
  unsigned int odepdb_id;
  dBodyID b_id;
  dGeomID g_id;
  t_real radius;
  t_real radius_effective;
  t_real mass;
  t_real *Lk;
  t_matrix inertia;

  /* mode-specific parameters */
  t_real bead_turns;
  t_real bead_angular_k;
  t_real bead_pos_k;
  t_real bead_align_k;
  t_real rho;
  t_real trap_x0;
  t_real trap_y0;
  t_real trap_z0;
  t_real bead_force;
  t_real bead_torque;
};
typedef struct bead_ode bead_ode;

void bead_init (bead_ode *bead,
    unsigned int odepdb_id,
    FILE *odepdb_file,
    dWorldID world_id,
    dSpaceID space_id,
    int finite_rotation,
    int gyroscopic_mode,
    t_real scale_l,
    t_real scale_e,
    t_real scale_f,
    config_t *simconfig);

void bead_force_and_torque (bead_ode *bead);

#endif
