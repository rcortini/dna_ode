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

#ifndef __SIMCONTEXT_H__
#define __SIMCONTEXT_H__



#include "config_read.h"
#include "dna_ode_core.h"
#include "langevin.h"
#include "mechanical_models.h"
#include "potentials.h"
#include "external_potentials.h"
#include "objects_creation.h"
#include "simparameters.h"



/* this structure contains all the information that is needed to perform the simulation. It contains
 * pointers both to the "simparameters" and to the "simstatistics" data structures, along with
 * other important data */
struct simcontext {
  /* pointers to simulation parameters and statistics */
  simparameters * simparams;

  /* ODE general objects */
  dWorldID world_id;
  dSpaceID space_id;
  dJointGroupID contactjoints;

  /* number of collisions */
  unsigned int numc;

  /* number of bodies in simulation */
  unsigned int nbodies;

  /* thermostat */
  thermostat my_thermostat;

  /* mechanical DNA model */
  mechanical_model_f mechanical_model;

  /* list of bodies to thermostat */
  sim_body *body_list;

  /* a structure to contain the inter-DNA forces parameters */
  potential_f U;

  /* a structure to contain the external forces */
  external_potential_f external_U;
  surface_ode *surface;

  /* files for input and output */
  FILE * odepdb_file;
  FILE * restart_file;
  FILE * initial_configuration_file;
  FILE * trajectory_file;
  FILE * data_file;

  /* additional data required */
  t_real sim_clock_0;
  t_real sim_clock;
  t_real elapsed_time;
  unsigned long sim_tick;
};
typedef struct simcontext simcontext;

void init_simulation_context (simcontext *simcon, simparameters *simparams);

#endif
