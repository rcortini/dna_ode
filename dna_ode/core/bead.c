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

#include "bead.h"

/* function to initialize the bead */
void bead_init (bead_ode *bead,
    unsigned int odepdb_id,
    FILE *odepdb_file,
    dWorldID world_id,
    dSpaceID space_id,
    int finite_rotation,
    int gyroscopic_mode,
    t_real scale_l,
    t_real scale_e,
    t_real scale_f,
    config_t *simconfig) {
  dMass m_id;
  t_matrix inertia;
  int mode;

  /* init all bead variables */
  bead->bead_turns = 0.;
  bead->bead_angular_k = 0.;
  bead->bead_pos_k = 0.;
  bead->rho = 0.;
  bead->trap_x0 = 0.;
  bead->trap_y0 = 0.;
  bead->trap_z0 = 0.;
  bead->bead_force = 0.;
  bead->bead_torque = 0.;

  /* get mode */
  mode = bead->mode = mycfg_read_int (simconfig, "mode");

  /* if we are in "fixed n mode", we must get the number of turns. Otherwise we must get the torque */
  if (mode==FIXED_N) {
    bead->bead_angular_k = mycfg_read_double (simconfig, "bead_angular_k") * 1e-14 / scale_e; /* convert from pN nm to dyn cm */
    bead->bead_align_k = mycfg_read_double (simconfig, "bead_align_k") * 1e-14 / scale_e; /* convert from pN nm to dyn cm */
    bead->bead_force = mycfg_read_double (simconfig, "bead_force") * 1e-7 / scale_f; /* convert from pN to dyn */
    bead->bead_turns = mycfg_read_double (simconfig, "bead_turns");
  }
  else if (mode==FIXED_TORQUE) {
    bead->bead_force = mycfg_read_double (simconfig, "bead_force") * 1e-7 / scale_f; /* convert from pN to dyn */
    bead->bead_torque = mycfg_read_double (simconfig, "bead_torque") * 1e-14 / scale_e; /* convert from pN nm to dyn cm */
  }
  else if (mode==FIXED_EXT) {
    bead->rho = mycfg_read_double (simconfig, "rho");
    bead->bead_pos_k = mycfg_read_double (simconfig, "bead_pos_k") * scale_l/scale_f; /* 1 pN/nm = 1 dyn/cm  (wow...)*/
  }
  else {
    err_message ("Invalid \"mode\" parameter: is %d\n", mode);
    exit (EXIT_FAILURE);
  }

  /* assign parameters of the bead: these parameters are assigned in local units and don't need scaling */
  bead->mass = mycfg_read_double (simconfig, "bead_mass");
  bead->radius = mycfg_read_double (simconfig, "bead_radius");
  bead->radius_effective = mycfg_read_double (simconfig, "bead_radius_effective");

  /* assigns all the parameters of the bead */
  bead->b_id = dBodyCreate (world_id);
  dMassSetSphereTotal (&m_id, bead->mass, bead->radius_effective);
  dBodySetMass (bead->b_id, &m_id);
  dBodySetFiniteRotationMode (bead->b_id, finite_rotation);
  dBodySetGyroscopicMode (bead->b_id, gyroscopic_mode);
  bead->g_id = dCreateSphere (space_id, bead->radius);
  dGeomSetBody (bead->g_id, bead->b_id);

  /* set inverse inertia tensor */
  body_inertia_tensor (inertia, bead->b_id, 1);
  matrix_copy (bead->inertia, inertia);

  /* assigns odepdb identity */
  bead->odepdb_id = odepdb_id;

  /* set force and torque to zero */
  dBodySetForce (bead->b_id, 0., 0., 0.);
  dBodySetTorque (bead->b_id, 0., 0., 0.);

  /* writes to odpdb file */
  file_message (odepdb_file, "%d sphere %f\n",
      odepdb_id,
      bead->radius);
}


/* This function is called whenever a bead is present in the simulations,
 * before a call to WorldStep. It adds force and torque to a bead, depending
 * on which simulation mode was chosen */
void bead_force_and_torque (bead_ode *bead) {
  t_real x, y, z;
  t_vector bead_position, bead_t, z_axis;
  t_vector bead_force, bead_torque;
  const bead_mode mode = bead->mode;

  /* init z axis */
  vector_set (z_axis, 0., 0., 1.);

  switch (mode) {
    /* in FIXED_N mode we assign the force on the bead, and the torque as
     * t = -k (n-n_0) */
    case FIXED_N :
      body_axis (bead_t, bead->b_id, 2);

      /* calculate force and set it */
      vector_scale (bead_force, z_axis, bead->bead_force);
      body_add_force (bead->b_id, bead_force);

      /* calculate torque and set it */
      /* magnetic trap component */
      vector_scale (bead_torque, z_axis, -bead->bead_angular_k*(*bead->Lk - bead->bead_turns));
      body_add_rel_torque (bead->b_id, bead_torque);

      /* align component */
      vector_cross (bead_torque, bead_t, z_axis);
      dScaleVector3 (bead_torque, bead->bead_align_k);
      body_add_torque (bead->b_id, bead_torque);
      break;

    /* in FIXED_TORQUE mode we assign the force and the torque on the bead */
    case FIXED_TORQUE :
      /* calculate force and set it */
      vector_scale (bead_force, z_axis, bead->bead_force);
      body_add_force (bead->b_id, bead_force);

      /* calculate torque and set it */
      vector_scale (bead_torque, z_axis, bead->bead_torque);
      body_add_rel_torque (bead->b_id, bead_torque);
      break;

    /* in FIXED_EXT mode we don't assign torque and assign force from an isotropic
     * "optical trap": F = -k *(|r-r_0|) */
    case FIXED_EXT :
      /* get bead position */
      body_position (bead_position, bead->b_id);
      x = bead_position [0];
      y = bead_position [1];
      z = bead_position [2];

      /* calculate and assign force. The force is harmonic isotropic */
      vector_set (bead_force,
	  -bead->bead_pos_k * (x-bead->trap_x0),
	  -bead->bead_pos_k * (y-bead->trap_y0),
	  -bead->bead_pos_k * (z-bead->trap_z0));
      body_add_force (bead->b_id, bead_force);
      break;
    default :
      err_message ("No valid bead mode supplied\n");
      exit (EXIT_FAILURE);
      break;
  }
}
