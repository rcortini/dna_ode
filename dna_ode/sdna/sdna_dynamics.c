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
#include "sdna.h"



/* ****************************** */
/* SIMULATION STEP                */
/* ****************************** */



/* this function performs a single simulation step */
void sdna_simstep (simcontext * simcon, sdna_simcontext *sdna_simcon) {
  unsigned int i;
  const unsigned int N = sdna_simcon->dna->nsegments;
  t_vector ang_vel;

  /* MECHANICAL TORQUES */

  sdna_simcon->Tw = dna_twist_and_torque (sdna_simcon->dna, &simcon->mechanical_model);
  sdna_simcon->Tw /= 2.*M_PI;

  /* ENERGY */

  /* kinetic energy */
  sdna_simcon->Kt = 0.;
  for (i=0; i<N; i++) sdna_simcon->Kt += body_total_kinetic_energy (sdna_simcon->dna->seg [i]->b_id);
  sdna_simcon->Kt += body_total_kinetic_energy (sdna_simcon->bead->b_id);

  /* potential energy */
  if (simcon->simparams->forces_on) sdna_simcon->Ut = force_loop (&simcon->U, sdna_simcon->dna);

  /* total energy */
  sdna_simcon->Et = sdna_simcon->Kt + sdna_simcon->Ut;

  /* BEAD FORCE AND TORQUE */
  bead_force_and_torque (sdna_simcon->bead);

  /* THERMOSTAT */

  /* and invoke the thermostat */
  thermostat_all_bodies (&simcon->my_thermostat, simcon->body_list, sdna_simcon->dna->nsegments+1);

  /* DETECT COLLISIONS */
  simcon->numc = 0; /* reset the number of collisions */
  dSpaceCollide (simcon->space_id, simcon, &collision_callback);

  /* if we have no collisions, use the exact solution for the joints */
   if (simcon->numc==0) {
    for (i=0; i<sdna_simcon->dna->njoints; i++) dJointDisable (sdna_simcon->dna->joint [i]->j_id);
    sdna_joints_solution (simcon, sdna_simcon);
  }

  /* STEP THE WORLD */

  dWorldQuickStep (simcon->world_id, simcon->simparams->ode_step);
  if (simcon->contactjoints)
    dJointGroupEmpty (simcon->contactjoints);

  /* LINKING NUMBER CONSTRAINT */
  body_angular_velocity_principal_axes (ang_vel, sdna_simcon->bead->b_id);
  sdna_simcon->Lk += ang_vel [2]*(simcon->simparams->ode_step / (2.*M_PI));

  /* INCREMENT SIMULATION CLOCKS */

  simcon->sim_tick++;
  simcon->elapsed_time += simcon->simparams->ode_step;
  simcon->sim_clock = simcon->sim_clock_0 + simcon->elapsed_time;
}
