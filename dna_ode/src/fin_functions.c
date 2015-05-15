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
#include "simcontext.h"



/* closes the output files that were opened during the simulation interval */
void close_output_files (simcontext * simcon) {
  fclose (simcon->trajectory_file);
  fclose (simcon->data_file);
  fclose (simcon->odepdb_file);
}



/* frees memory allocated during system-independent initialization */
void free_memory (simcontext * simcon) {
  unsigned int i;

  /* free random number generator */
  my_free_RNG ();

  /* destroy ODE-related stuff */
  dSpaceDestroy (simcon->space_id);
  dJointGroupDestroy (simcon->contactjoints);
  dWorldDestroy (simcon->world_id);

  /* destroy thermostat */
  if (simcon->my_thermostat.type == GLOBAL) {
    free (simcon->my_thermostat.params);
  }

  /* destroy external forces memory */
  if (simcon->simparams->external_forces_on) {
    if (strcmp (simcon->simparams->external_potential_type, "harmonic_z")==0) {
      free (simcon->external_U.p);
    }
    if (strcmp (simcon->simparams->external_potential_type, "morse_z")==0) {
      free (simcon->external_U.p);
    }
    if (strcmp (simcon->simparams->external_potential_type, "gravity")==0) {
      free (simcon->surface);
    }
  }

  /* destroy mechanical model memory */
  if (simcon->mechanical_model.f == body_twist_and_torque_harmonic ||
      simcon->mechanical_model.f == body_twist_and_torque_kinkable) {
    free (simcon->mechanical_model.p);
  }

  /* destroy body_list */
  for (i=0; i<simcon->nbodies; i++) free (simcon->body_list [i].name);
  free (simcon->body_list);

  /* destroy simcon and simparams */
  free (simcon->simparams);
  free (simcon);
}
