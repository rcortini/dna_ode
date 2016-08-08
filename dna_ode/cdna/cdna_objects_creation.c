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



/* 
 * SIMULATED OBJECT CREATION FUNCTIONS
 */



/* creates the DNA in the simulation context */
dna_ode * cdna_create (simcontext *simcon) {
  int i, nsegments, njoints;
  t_vector center_i, t0, u0, v0;
  dna_ode *dna;
  t_real alpha = 2*M_PI/simcon->simparams->dna_nsegments;
  t_real R = simcon->simparams->dna_segment_length / sin (alpha);
  t_real h = R * cos (alpha);

  /* memory allocation */
  dna = (dna_ode *) malloc (sizeof (dna_ode));

  /* get the number of segments and of joints from the simulation context. For convenience */
  nsegments = simcon->simparams->dna_nsegments;
  njoints = simcon->simparams->dna_nsegments;

  /* assign the number of segments and of joints to the DNA data structure */
  dna->nsegments = nsegments;
  dna->njoints = njoints;

  /* allocate the memory for all the simulated DNA segments */
  dna->seg = (dna_segment **) malloc (nsegments * sizeof (dna_segment *));
  dna->joint = (joint_ode **) malloc (njoints * sizeof (joint_ode *));

  /* assign the zero linking number of the DNA */
  dna->Lk0 = (t_real ) PAR(dna_nbasepairs) / simcon->simparams->dna_basepairs_per_turn;

  /* assigns the direction vectors */
  vector_set (u0, 0., 0, -1.);
  vector_set (v0, 0., 1., 0.);
  vector_set (t0, 1., 0., 0.);

  /* creates the DNA as an array of cylinders connected by "Ball-in-Socket" type of joints */
  for (i=0; i<nsegments; i++) {
    /* calculates the position of the center of mass of segment i */
    vector_set (center_i, (t_real) (i + 0.5), -h, 0.);

    /* assigns position and orientation of the i-th DNA segment */
    dna->seg [i] = dna_segment_create (i,
	simcon->odepdb_file,
	simcon->world_id,
	simcon->space_id,
	PAR (dna_segment_mass),
	PAR (dna_effective_rad),
	PAR (dna_segment_length),
	3,
	PAR (finite_rotation),
	PAR (gyroscopic_mode),
	center_i,
	u0,
	v0,
	t0);
    dna->seg [i]->id = i;
  }

  return dna;
}



/* creates and assigns all the DNA joints */
void cdna_joints_create (simcontext * simcon, dna_ode * dna) {
  unsigned int i;
  const unsigned int N = dna->nsegments;

  /* this is a cycle on all the joints */
  for (i=0; i<dna->njoints; i++) {
    dBodyID *previous, *next;
    t_vector position;

    /* memory allocation for both the joint and the joint feedback */
    dna->joint [i] = (joint_ode *) malloc (sizeof (joint_ode));

    /* position of the joint */
    body_point_position (position, dna->seg [i]->b_id, 0., 0., PAR (dna_segment_length)/2.);

    /* previous and next bodies */
    previous = &dna->seg [i]->b_id;
    next = &dna->seg [(i+1)%N]->b_id; /* circularity indices */

    /* creates the joint. The function will know what joint type to create based on its identity */
    joint_create (i, simcon->world_id, previous, next, dna->joint [i], position, 1);
  }
}
