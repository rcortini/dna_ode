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
#include "fdna.h"



/* 
 * FINALIZATION FUNCTIONS
 */



void fdna_free_memory (fdna_simcontext *fdna_simcon) {
  unsigned int i;
  for (i=0; i<fdna_simcon->dna->nsegments; i++) free (fdna_simcon->dna->seg [i]);
  for (i=0; i<fdna_simcon->dna->njoints; i++) free (fdna_simcon->dna->joint [i]);
  free (fdna_simcon->dna->seg);
  free (fdna_simcon->dna->joint);
  free (fdna_simcon->dna);
  fdna_joint_solution_workspace_free (fdna_simcon->ws);
  free (fdna_simcon);
}
