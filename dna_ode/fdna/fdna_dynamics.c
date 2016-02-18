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



/* this function performs a single simulation step */
void fdna_simstep (simcontext * simcon, fdna_simcontext * fdna_simcon) {
  unsigned int i;
  const unsigned int N = fdna_simcon->dna->nsegments;

  /* MECHANICAL TORQUES */

  fdna_simcon->dTw = dna_twist_and_torque (fdna_simcon->dna, &simcon->mechanical_model);
  fdna_simcon->dTw /= 2.*M_PI;

  /* ENERGY */

  /* kinetic energy */
  fdna_simcon->Kt = 0.;
  for (i=0; i<N; i++) fdna_simcon->Kt += body_total_kinetic_energy (fdna_simcon->dna->seg [i]->b_id);

  /* potential energy */
  if (simcon->simparams->forces_on) fdna_simcon->Ut = force_loop (&simcon->U, fdna_simcon->dna);

  /* external energy */
  if (simcon->simparams->external_forces_on) fdna_simcon->ext_Ut = external_forces_add (&simcon->external_U, fdna_simcon->dna);

  /* total energy */
  fdna_simcon->Et = fdna_simcon->Kt + fdna_simcon->Ut;

  /* THERMOSTAT */

  thermostat_all_bodies (&simcon->my_thermostat, simcon->body_list, fdna_simcon->dna->nsegments);

  /* WORLD STEP */

  simcon->numc = 0; /* reset the number of collisions */
  dSpaceCollide (simcon->space_id, simcon, &collision_callback);
  if (simcon->numc == 0) {
    for (i=0; i<fdna_simcon->dna->njoints; i++) dJointDisable (fdna_simcon->dna->joint [i]->j_id);
    fdna_joints_solution (simcon, fdna_simcon);
  }
  dWorldQuickStep (simcon->world_id, simcon->simparams->ode_step);
  if (simcon->contactjoints)
    dJointGroupEmpty (simcon->contactjoints);

  /* INCREMENT SIMULATION CLOCKS */

  simcon->sim_tick++;
  simcon->elapsed_time += simcon->simparams->ode_step;
  simcon->sim_clock = simcon->sim_clock_0 + simcon->elapsed_time;
}
