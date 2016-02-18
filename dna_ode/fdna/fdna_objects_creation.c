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



/* 
 * SIMULATED OBJECT CREATION FUNCTIONS
 */



/* creates the DNA in the simulation context */
dna_ode * fdna_create (simcontext *simcon) {
  int i, nsegments, njoints;
  t_vector center_i, t0, u0, v0;
  dna_ode *dna;

  /* memory allocation */
  dna = (dna_ode *) malloc (sizeof (dna_ode));

  /* get the number of segments and of joints from the simulation context. For convenience */
  nsegments = simcon->simparams->dna_nsegments;
  njoints = simcon->simparams->dna_nsegments-1;

  /* assign the number of segments and of joints to the DNA data structure */
  dna->nsegments = nsegments;
  dna->njoints = njoints;

  /* allocate the memory for all the simulated DNA segments */
  dna->seg = (dna_segment **) malloc (nsegments * sizeof (dna_segment *));
  dna->joint = (joint_ode **) malloc (njoints * sizeof (joint_ode *));

  /* assigns the direction vectors */
  vector_set (u0, 1., 0., 0.);
  vector_set (v0, 0., 1., 0.);
  vector_set (t0, 0., 0., 1.);

  /* creates the DNA as an array of cylinders connected by "Ball-in-Socket" type of joints */
  for (i=0; i<nsegments; i++) {
    /* calculates the position of the center of mass of segment i */
    vector_set (center_i, 0., 0., (t_real) i);

    /* assigns position and orientation of the i-th DNA segment */
    dna->seg [i] = dna_segment_create (i,
	simcon->odepdb_file,
	simcon->world_id,
	simcon->space_id,
	simcon->simparams->dna_segment_mass,
	simcon->simparams->dna_effective_rad,
	simcon->simparams->dna_segment_length,
	3,
	simcon->simparams->finite_rotation,
	simcon->simparams->gyroscopic_mode,
	center_i, u0, v0, t0);

    dna->seg [i]->id = i;
  }

  return dna;
}



/* creates and assigns all the DNA joints */
void fdna_joints_create (simcontext * simcon, dna_ode * dna) {
  unsigned int i;

  /* this is a cycle on all the joints */
  for (i=0; i<dna->njoints; i++) {
    t_vector joint_position;

    /* memory allocation for both the joint and the joint feedback */
    dna->joint [i] = (joint_ode *) malloc (sizeof (joint_ode));

    /* now calculate the position of the joint from the point of view of the "previous" and "next" segments */
    segment_next_joint (joint_position, dna->seg [i]);

    /* creates the joint. The function will know what joint type to create based on its identity */
    joint_create (i, simcon->world_id, &dna->seg [i]->b_id, &dna->seg [i+1]->b_id, dna->joint [i], joint_position, 1);
    dna->joint [i]->next_body_size = simcon->simparams->dna_segment_length;
  }
}
