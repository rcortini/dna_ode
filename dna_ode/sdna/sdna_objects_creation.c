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
#include "sdna.h"



/* creates the DNA in the simulation context */
dna_ode * sdna_create (simcontext *simcon) {
  int i, nsegments, njoints;
  t_vector center_i, t0, u0, v0;
  dna_ode *dna;

  /* memory allocation */
  dna = (dna_ode *) malloc (sizeof (dna_ode));

  /* get the number of segments and of joints from the simulation context. For convenience */
  nsegments = simcon->simparams->dna_nsegments;
  njoints = simcon->simparams->dna_nsegments + 1;

  /* assign the number of segments and of joints to the DNA data structure */
  dna->nsegments = nsegments;
  dna->njoints = njoints;

  /* assign the zero linking number of the DNA */
  dna->Lk0 = dna->nsegments / simcon->simparams->dna_basepairs_per_turn;

  /* allocate the memory for all the simulated DNA segments */
  dna->seg = (dna_segment **) malloc (nsegments * sizeof (dna_segment *));
  dna->joint = (joint_ode **) malloc (njoints * sizeof (joint_ode *));

  /* assigns the direction vectors */
  vector_set (u0, 0., 0, -1.);
  vector_set (v0, 0., 1., 0.);
  vector_set (t0, 1., 0., 0.);

  /* creates the DNA as an array of cylinders connected by "Ball-in-Socket" type of joints */
  for (i=0; i<nsegments; i++) {
    /* calculates the position of the center of mass of segment i */
    vector_set (center_i, (t_real) (i + 0.5), 0., 0.);

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

  /* NOTE: we do not create the joints here because we do not know yet the position of the anchoring point on the bead */
  return dna;
}



/* creates and assigns all the DNA joints */
void sdna_joints_create (simcontext * simcon, sdna_simcontext *sdna_simcon, dna_ode * dna) {
  const unsigned int nsegments = sdna_simcon->dna->nsegments;
  unsigned int i;

  /* this is a cycle on all the joints */
  for (i=0; i<dna->njoints; i++) {
    dBodyID *previous, *next;
    t_vector joint_position;

    /* memory allocation for both the joint and the joint feedback */
    dna->joint [i] = (joint_ode *) malloc (sizeof (joint_ode));

    /* case of the DNA-surface joint */
    if (i==0) {
      /* first we evaluate what is the position of the joint */
      body_point_position (joint_position, dna->seg [0]->b_id, 0., 0., -.5*sdna_simcon->dna->seg [0]->length);
      previous = &sdna_simcon->surface->b_id;
      next = &sdna_simcon->dna->seg [0]->b_id;
      joint_create (i, simcon->world_id, previous, next, dna->joint [i], joint_position, 0);
      dna->joint [i]->next_body_size = dna->seg [0]->length;
    }
    /* the case of the DNA-bead joint */
    else if (i==nsegments) {
      /* calculates the position of the joint */
      segment_next_joint (joint_position, dna->seg [nsegments-1]);
      previous = &sdna_simcon->dna->seg [sdna_simcon->dna->nsegments-1]->b_id;
      next = &sdna_simcon->bead->b_id;
      joint_create (i, simcon->world_id, previous, next, dna->joint [i], joint_position, 1);
      dna->joint [i]->next_body_size = sdna_simcon->bead->radius*2.;
    }
    else {
      /* gets position and orientation of the two bodies to connect */
      segment_next_joint (joint_position, dna->seg [i-1]);
      previous = &sdna_simcon->dna->seg [i-1]->b_id;
      next = &sdna_simcon->dna->seg [i]->b_id;
      joint_create (i, simcon->world_id, previous, next, dna->joint [i], joint_position, 1);
      dna->joint [i]->next_body_size = dna->seg [i]->length;
    }
  }
}
