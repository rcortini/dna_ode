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
#include "cdna.h"



/* this function performs a single simulation step */
void cdna_simstep (simcontext * simcon, cdna_simcontext *cdna_simcon) {
  unsigned int i;
  const unsigned int N = cdna_simcon->dna->nsegments;

  /* MECHANICAL TORQUES */

  cdna_simcon->dTw = dna_twist_and_torque (cdna_simcon->dna, &simcon->mechanical_model);
  cdna_simcon->dTw /= 2.*M_PI;

  /* ENERGY */

  /* kinetic energy */
  cdna_simcon->Kt = 0.;
  for (i=0; i<N; i++) cdna_simcon->Kt += body_total_kinetic_energy (cdna_simcon->dna->seg [i]->b_id);

  /* potential energy */
  if (simcon->simparams->forces_on) cdna_simcon->Ut = force_loop (&simcon->U, cdna_simcon->dna);

  /* external energy */
  if (simcon->simparams->external_forces_on) cdna_simcon->ext_Ut = external_forces_add (&simcon->external_U, cdna_simcon->dna);

  /* total energy */
  cdna_simcon->Et = cdna_simcon->Kt + cdna_simcon->Ut;

  /* THERMOSTAT */

  thermostat_all_bodies (&simcon->my_thermostat, simcon->body_list, cdna_simcon->dna->nsegments);

  /* WORLD STEP */

  simcon->numc = 0; /* reset the number of collisions */
  dSpaceCollide (simcon->space_id, simcon, &collision_callback);
  if (simcon->numc == 0) {
    for (i=0; i<cdna_simcon->dna->nsegments; i++) dJointDisable (cdna_simcon->dna->joint [i]->j_id);
    ccdna_joints_solution (simcon, cdna_simcon);
  }
  dWorldQuickStep (simcon->world_id, simcon->simparams->ode_step);
  if (simcon->contactjoints)
    dJointGroupEmpty (simcon->contactjoints);

  /* INCREMENT SIMULATION CLOCKS */

  simcon->sim_tick++;
  simcon->elapsed_time += simcon->simparams->ode_step;
  simcon->sim_clock = simcon->sim_clock_0 + simcon->elapsed_time;
}
