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

#include "langevin.h"

/* initializes a sim_body giving the parameters that
 * will be used in the thermostat routine */
void sim_body_init (sim_body *body,
    unsigned int odepdb_id,
    char *name,
    t_real sigma_langevin,
    t_real ode_step,
    dBodyID *b_id) {
  unsigned int i;
  t_real mass, trans_friction, trans_sigma;
  t_vector inertia, rot_friction, rot_sigma;

  /* get mass and inertia vector (principal components) */
  mass = body_mass (*b_id);
  body_inertia (inertia, *b_id);

  /* set values of translational and rotational friction and sigma
   * based on the mass, the inertia, ode_step, and sigma_langevin */
  trans_friction = mass * sigma_langevin;
  trans_sigma = sqrt (2.* trans_friction/ode_step);
  vector_scale (rot_friction, inertia, sigma_langevin);
  for (i=0; i<3; i++) rot_sigma [i] = sqrt (2.* rot_friction [i]/ode_step);

  /* assign values to the sim_body structure */
  body->odepdb_id = odepdb_id;
  body->name = name;
  body->b_id = b_id;
  body->trans_friction = trans_friction;
  body->trans_sigma = trans_sigma;
  body->rot_friction_x = rot_friction [0];
  body->rot_friction_y = rot_friction [1];
  body->rot_friction_z = rot_friction [2];
  body->rot_sigma_x = rot_sigma [0];
  body->rot_sigma_y = rot_sigma [1];
  body->rot_sigma_z = rot_sigma [2];
}

/* initializes the thermostat. Before invoking this function, user needs to:
 * 1. set sigma_langevin
 * 2. set ode_step
 * 3. set beta
 * 4. set mE 
 * 5. set dna_seg_(rot,trans)_(sigma,friction) */
/* WARNING: this is a trick! Since in the programs where there is a bead the Langevin
 * parameters are assigned to be exactly the same as in the case of a dna segment, this
 * function will work in any case. But this is a trick */
void init_langevin_thermostat (thermostat *my_thermostat,
    char *therm_type,
    t_real *Kt,
    t_real sigma_langevin,
    t_real ode_step,
    t_real beta,
    t_real mE) {

  /* assign the thermostat type and parameters */
  if (strcmp (therm_type,"global")==0) {
    struct global_thermostat_params *p = (struct global_thermostat_params *) malloc (sizeof (struct global_thermostat_params));

    /* assign type and parameters */
    my_thermostat->type = GLOBAL;
    my_thermostat->params = p;

    /* assign parameters */
    p->sigma_langevin = sigma_langevin;
    p->ode_step = ode_step;
    p->beta = beta;
    p->Kt = Kt;
    p->mE = mE;
  }
  else if (strcmp (therm_type,"local")==0) {
    /* assign type and parameters */
    my_thermostat->type = LOCAL;
    my_thermostat->params = NULL;
  }
  else if (strcmp (therm_type,"off")==0) {
    my_thermostat->type=OFF;
    /* no parameters to assign in this case! */
  }
  else {
    err_message ("Invalid thermostat type: %s\n", therm_type);
    exit (EXIT_FAILURE);
  }
}



/* THERMOSTAT */
void thermostat_all_bodies (thermostat *my_thermostat, sim_body *body_list, const unsigned int nbodies) {

  if (my_thermostat->type!=OFF) {
    unsigned int i;

    /* -----------------------
     * GLOBAL THERMOSTAT 
     ----------------------- */
    if (my_thermostat->type==GLOBAL) {
      t_real scale_factor;
      t_vector axis, vel, omega, torque, inertia;
      t_real m, mE, Kt, beta, ode_step, sigma_langevin;

      /* dereference the thermostat parameters to "global" type */
      struct global_thermostat_params *p = (struct global_thermostat_params *) my_thermostat->params;

      /* get the parameters */
      mE = p->mE;
      Kt = *p->Kt;
      ode_step = p->ode_step;
      beta = p->beta;
      sigma_langevin = p->sigma_langevin;

      /* now we calculate the scale factor that is needed for the Langevin-Euler global thermostat step.
       * This calculation is performed once for every simulation step */
      scale_factor = sigma_langevin * ((1. - 1. / (2.*mE * beta)) * (mE/Kt)-1.) + sqrt (sigma_langevin / (beta * Kt * ode_step)) * ran_gaussian_01 ();

      /* cycle on all bodies */
      for (i=0; i<nbodies; i++) {
	unsigned int j;

	/* step the next body */
	dBodyID b = *body_list [i].b_id;

	/* get the mass */
	m = body_mass (b);

	/* rescale the linear velocity of the body */
	body_velocity (vel, b);
        dScaleVector3 (vel, scale_factor*m);
	body_add_force (b, vel);

	/* rescale the angular velocity of the body */
	body_angular_velocity (omega, b);
	body_inertia (inertia, b);

	/* calculate the torque to assign to the body */
	vector_set_zero (torque);
	for (j=0; j<3; j++) {
	  t_vector aj;
	  vector_set_zero (aj);
	  body_axis (axis, b, j);
          vector_scale (aj, axis, scale_factor*inertia[j]*vector_dot (omega, axis));
	  vector_sum (torque, torque, aj);
	}

	/* now add the calculated torque to the body */
	body_add_torque (b, torque);
      }
    }

    /* -----------------------
     * LOCAL THERMOSTAT 
     ----------------------- */
    else if (my_thermostat->type==LOCAL) {
      for (i=0; i<nbodies; i++) {
	sim_body b = body_list [i];
	langevin_euler_step (&b);
      }
    }
    else {
      err_message ("Error: invalid thermostat type\n");
      exit (EXIT_FAILURE);
    }
  }
}


/* this function adds force and torque to a body in the Langevin-Euler dynamics */
void langevin_euler_step (sim_body *body) {
  const dBodyID b = *body->b_id;
  const t_real trans_friction = body->trans_friction;
  const t_real trans_sigma = body->trans_sigma;
  t_vector vel, omega;

  /* get the linear and angular velocities of the body */
  body_velocity (vel, b);
  body_angular_velocity_principal_axes (omega, b);

  /* add the computed force to the body */
  dBodyAddForce (b, trans_sigma * ran_gaussian_01 () - trans_friction * vel[0],
                    trans_sigma * ran_gaussian_01 () - trans_friction * vel[1],
                    trans_sigma * ran_gaussian_01 () - trans_friction * vel[2]);

  /* add the computed torque to the body */
  dBodyAddRelTorque (b, body->rot_sigma_x * ran_gaussian_01 () - body->rot_friction_x * omega[0],
                        body->rot_sigma_y * ran_gaussian_01 () - body->rot_friction_y * omega[1],
                        body->rot_sigma_z * ran_gaussian_01 () - body->rot_friction_z * omega[2]);

}



/* this function performs a single step for the Langevin-Euler global thermostat */
void langevin_euler_global_step (dBodyID b, t_real scale_factor) {
  unsigned int j;
  t_vector axis, force, torque, inertia, omega;
  t_real m = body_mass (b);

  /* rescale the linear velocity of the body */
  body_velocity (force, b);
  dScaleVector3 (force, scale_factor*m);
  body_add_force (b, force);

  /* rescale the angular velocity of the body */
  body_angular_velocity (omega, b);
  body_inertia (inertia, b);

  /* calculate the torque to assign to the body */
  vector_set_zero (torque);
  for (j=0; j<3; j++) {
    t_vector aj;
    vector_set_zero (aj);
    body_axis (axis, b, j);
    vector_scale (aj, axis, scale_factor*inertia[j]*vector_dot (omega, axis));
    vector_sum (torque, torque, aj);
  }

  /* now add the calculated torque to the body */
  body_add_torque (b, torque);
}
