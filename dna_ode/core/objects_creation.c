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

#include "objects_creation.h"

/* 
 * SIMULATED OBJECT CREATION FUNCTIONS
 */

/* creates the surface */
surface_ode * surface_create (dSpaceID space_id, t_real a, t_real b, t_real c, t_real d) {
  surface_ode *surface = malloc (sizeof (surface_ode));

  /* assigns all the parameters of the surface */
  surface->a = a;
  surface->b = b;
  surface->c = c;
  surface->d = d;
  surface->g_id = dCreatePlane (space_id, surface->a, surface->b, surface->c, surface->d);
  surface->b_id = SURFACE_BODY_ID; /* this should be correct as this is a non-placeable geometry/body */

  /* returns a pointer to the newly created object */
  return surface;
}



/* returns a DNA segment that has "center" as center of mass, and (u, v, t) as principal axes */
dna_segment * dna_segment_create (unsigned int odepdb_id,
    FILE *odepdb_file,
    dWorldID world_id,
    dSpaceID space_id,
    t_real mass,
    t_real radius,
    t_real length,
    unsigned int axis,
    int finite_rotation,
    int gyroscopic_mode,
    t_real *center,
    t_real *u,
    t_real *v,
    t_real *t) {
  dMass m_id;
  t_matrix inertia;
  dna_segment *seg = (dna_segment *) malloc (sizeof (dna_segment));

  /* creates the body */
  seg->b_id = dBodyCreate (world_id);

  /* creates the mass and assigns the mass to be that of a cylinder */
  dMassSetCylinderTotal (&m_id, mass, axis, radius, length);

  /* assigns the mass to the body */
  dBodySetMass (seg->b_id, &m_id);

  /* sets the rotation mode and the gyroscopic mode to the body */
  dBodySetFiniteRotationMode (seg->b_id, finite_rotation);
  dBodySetGyroscopicMode (seg->b_id, gyroscopic_mode);

  /* creates the geometry and assigns the geometry to the body */
  seg->g_id = dCreateCylinder (space_id, radius, length);
  dGeomSetBody (seg->g_id, seg->b_id);

  /* assigns internal parameters to the segment */
  seg->mass = mass;
  seg->radius = radius;
  seg->length = length;

  /* assigns position to the segment */
  body_set_position (seg->b_id, center);

  /* assigns orientation to the segment and inertia matrix */
  body_set_axes (seg->b_id, u, v, t);
  body_inertia_tensor (inertia, seg->b_id, 1);
  matrix_copy (seg->inertia, inertia);

  /* assigns odepdb identity */
  seg->odepdb_id = odepdb_id;

  /* writes to odpdb file */
  file_message (odepdb_file, "%d cylinder %f %f %d\n",
      odepdb_id,
      radius,
      length,
      axis);

  return seg;
}



/* this function initializes and assigns all the properties to a joint that has already been
 * malloc'd externally to the function. This function is invoked from the dna_joints_create function */
void joint_create (unsigned int id,
    dWorldID world_id,
    dBodyID *prev,
    dBodyID *next,
    joint_ode *joint,
    const t_real *position,
    unsigned int assign_feedback) {
  /* assign the joint identity */
  joint->id = id;

  /* we attach it to the DNA and the surface */
  joint->j_id = dJointCreateBall (world_id, 0);
  dJointAttach (joint->j_id, *prev, *next);
  dJointSetBallAnchor (joint->j_id, position[0], position[1], position[2]);

  /* we assign the "previous" and "next" body pointers */
  joint->previous_body = prev;
  joint->next_body = next;

  /* we assign the joint feedback */
  if (assign_feedback) dJointSetFeedback (joint->j_id, &joint->jf_id);

  /* initializes twist */
  joint->twist = 0.;
}
