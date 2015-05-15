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
#include "sdna.h"



/* 
 * FINALIZATION FUNCTIONS
 */



/* free all the memory associated with the simulation */
void sdna_free_memory (sdna_simcontext * sdna_simcon) {
  unsigned int i;

  /* free memory */
  sdna_joint_solution_workspace_free (sdna_simcon->ws);
  for (i=0; i<sdna_simcon->dna->nsegments; i++) free (sdna_simcon->dna->seg [i]);
  for (i=0; i<sdna_simcon->dna->njoints; i++) free (sdna_simcon->dna->joint [i]);
  free (sdna_simcon->dna->seg);
  free (sdna_simcon->dna->joint);
  free (sdna_simcon->dna);
  free (sdna_simcon->bead);
  free (sdna_simcon->surface);
  free (sdna_simcon);
}
