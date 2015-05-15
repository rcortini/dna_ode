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

#include "config.h"
#include "dna_ode_core.h"
#include "dna.h"
#include "dna_ode.h"



/* initializes the simulation context */
void init_simulation_context (simcontext *simcon, simparameters *simparams) {
  int hs_minlev, hs_maxlev;

  /* SIMCONTEXT INITIALIZATION */

  simcon->simparams = simparams;

  /* ODE INITIALIZATION */

  /* creation of space, world, and contact joint group */
  simcon->world_id = dWorldCreate ();
  dWorldSetGravity (simcon->world_id, 0., 0., 0.); /* gravity-less world */
  dWorldSetERP (simcon->world_id, simparams->erp_glob);
  dWorldSetCFM (simcon->world_id, simparams->cfm_glob);
  simcon->contactjoints = dJointGroupCreate (MAX_CONTACTS); /* note that it takes as an argument a max_size, but now it is not used */

  /* creates the ODE hash space, with its optimal levels */
  simcon->space_id = dHashSpaceCreate (0);
  hs_minlev = 0.;
  hs_maxlev = (int) (gsl_sf_log (simparams->dna_segment_length)/gsl_sf_log (2.)) + 1;
  dHashSpaceSetLevels (simcon->space_id, hs_minlev, hs_maxlev);

  /* assign the "successive over-relaxation" (SOR) algorithm parameters to ODE */
  dWorldSetQuickStepNumIterations (simcon->world_id, simparams->sor_numiter);
  dWorldSetQuickStepW (simcon->world_id, simparams->sor_omega);

  /* set the number of simulated bodies */
  simcon->nbodies = system_nbodies (simparams->system, simparams->dna_nsegments);

  /* sets to zero the number of contacts */
  simcon->numc = 0;

  /* GB and GT */

  dna_gb_gt (simparams->dna_kuhn_length, simparams->dna_segment_length, simparams->dna_ltp, simparams->tol, &simparams->gb, &simparams->gt);

  /* INTER-DNA POTENTIAL */

  /* check if the forces are on, and what potential was chosen by the user */
  if (simparams->forces_on) {
    if (strcmp (simparams->potential_type, "LJ")==0) {
      struct LJ_parameters *par = malloc (sizeof (struct LJ_parameters));
      
      /* assign the LJ parameters */
      LJ_parameters_calculate (simcon->simparams->dna_segment_length,
	  			 simcon->simparams->LJ_sigma,
	  			 simcon->simparams->LJ_eps,
				 par);

      /* assign the potential function */
      simcon->U.f = &LJ_pair_potential_and_force;
      simcon->U.p = par;
    }

    else if (strcmp (simparams->potential_type, "Todd")==0) {
      struct Todd_parameters *par = malloc (sizeof (struct Todd_parameters));

      /* calculate the parameters for the Todd potential from the user-supplied parameters */
      Todd_parameters_calculate (simcon->simparams->dna_segment_length,
	  			 simcon->simparams->CA,
				 simcon->simparams->r_eq,
				 simcon->simparams->lambda,
				 par);

      /* and assign the potential function */
      simcon->U.f = &Todd_pair_potential_and_force;
      simcon->U.p = par;
    }
    else if (strcmp (simparams->potential_type, "LJ_trunc")==0) {
      struct LJ_trunc_parameters *par = malloc (sizeof (struct LJ_trunc_parameters));
      
      /* assign the LJ parameters */
      LJ_trunc_parameters_calculate (simcon->simparams->dna_segment_length,
	  			 simcon->simparams->LJ_sigma,
	  			 simcon->simparams->LJ_eps,
	  			 simcon->simparams->r_trunc,
				 par);

      /* assign the potential function */
      simcon->U.f = &LJ_trunc_pair_potential_and_force;
      simcon->U.p = par;
    }
  }

  /* EXTERNAL FORCES */

  if (simparams->external_forces_on) {
    if (strcmp (simparams->external_potential_type, "harmonic_z")==0) {
      struct harmonic_z_parameters *par = malloc (sizeof (struct harmonic_z_parameters));

      /* assign the parameters to the function */
      par->k = simparams->k_harmonic;
      par->z0 = simparams->z0_harmonic;

      /* and assign the external potential function */
      simcon->external_U.f = &harmonic_z_potential;
      simcon->external_U.p = par;
    }
    else if (strcmp (simparams->external_potential_type, "morse_z")==0) {
      struct morse_z_parameters *par = malloc (sizeof (struct morse_z_parameters));

      /* assign the parameters to the function */
      par->eps = simparams->eps_morse;
      par->z0 = simparams->z0_morse;
      par->lambda = simparams->lambda_morse;

      /* and assign the external potential function */
      simcon->external_U.f = &morse_z_potential;
      simcon->external_U.p = par;
    }
    else if (strcmp (simparams->external_potential_type, "gravity")==0) {
      dWorldSetGravity (simcon->world_id, 0., 0., simparams->gravity_val);
      simcon->surface = surface_create (simcon->space_id, 0., 0., 1., simparams->gravity_z0);
      simcon->external_U.f = NULL;
      simcon->external_U.p = NULL;
    }
    else if (strcmp (simparams->external_potential_type, "morse_z_surface")==0) {
      struct morse_z_parameters *par = malloc (sizeof (struct morse_z_parameters));

      /* assign the parameters to the function */
      par->eps = simparams->eps_morse;
      par->z0 = simparams->z0_morse;
      par->lambda = simparams->lambda_morse;
      simcon->surface = surface_create (simcon->space_id, 0., 0., 1., simparams->z0_surface);

      /* and assign the external potential function */
      simcon->external_U.f = &morse_z_potential;
      simcon->external_U.p = par;
    }
  }

  /* MECHANICAL MODEL */

  if (strcmp (simparams->mechanical_model, "kinkable")==0) {
    t_real *pars = (t_real *) malloc (2 * sizeof (t_real));
    pars [0] = simparams->gb;
    pars [1] = simparams->gt;

    /* assign mechanical model */
    simcon->mechanical_model.f = body_twist_and_torque_kinkable;
    simcon->mechanical_model.p = pars;
  }
  else if (strcmp (simparams->mechanical_model, "harmonic")==0) {
    t_real *pars = (t_real *) malloc (2 * sizeof (t_real));
    pars [0] = simparams->gb;
    pars [1] = simparams->gt;

    /* assign mechanical model */
    simcon->mechanical_model.f = body_twist_and_torque_harmonic;
    simcon->mechanical_model.p = pars;
  }
  else if (strcmp (simparams->mechanical_model, "harmonic4")==0) {
    t_real *pars = (t_real *) malloc (3 * sizeof (t_real));
    pars [0] = simparams->gb;
    pars [1] = simparams->gt;
    pars [2] = simparams->k4;

    /* assign mechanical model */
    simcon->mechanical_model.f = body_twist_and_torque_harmonic4;
    simcon->mechanical_model.p = pars;
  }
  else if (strcmp (simparams->mechanical_model, "nicked_kinkable")==0) {
    simcon->mechanical_model.f = body_twist_and_torque_nicked_kinkable;
    simcon->mechanical_model.p = &simparams->gb;
  }
  else if (strcmp (simparams->mechanical_model, "nicked_harmonic")==0) {
    simcon->mechanical_model.f = body_twist_and_torque_nicked_harmonic;
    simcon->mechanical_model.p = &simparams->gb;
  }
  else if (strcmp (simparams->mechanical_model, "nicked_harmonic4")==0) {
    t_real *pars = (t_real *) malloc (2 * sizeof (t_real));
    pars [0] = simparams->gb;
    pars [1] = simparams->k4;

    /* assign mechanical model */
    simcon->mechanical_model.f = body_twist_and_torque_nicked_harmonic4;
    simcon->mechanical_model.p = &simparams->gb;
  }
  else {
    err_message ("Invalid mechanical model \"%s\"", simparams->mechanical_model);
    exit (EXIT_FAILURE);
  }

  /* OUTPUT FILES HANDLING */

  /* open the output files */
  simcon->odepdb_file = my_fopen (simcon->simparams->odepdb_file_name, "w");
  simcon->trajectory_file = my_fopen (simcon->simparams->trajectory_file_name, "w");
  simcon->data_file = my_fopen (simcon->simparams->data_file_name, "w");

  /* CLOCK START */

  simcon->elapsed_time = 0.;
  simcon->sim_tick = 0;
  simcon->sim_clock = 0;

  /* PRINT CONFIGURATION TO TERMINAL */
  {
    char hostname [STR_BUFF_SIZE];
    time_t now;
    struct tm tstruct;
    char time_string [STR_BUFF_SIZE];

    /* get the current time and date and convert it to a human-readable string */
    now = time (0);
    tstruct = *localtime (&now);
    strftime (time_string, STR_BUFF_SIZE, "%Y-%m-%d %X", &tstruct);

    /* start time and host name */
    gethostname (hostname, STR_BUFF_SIZE);
    message ("Job ran on %s at %s\n", hostname, time_string);

    /* program name and version */
    message ("Package: " PACKAGE_STRING "\n");
    message ("Program: %s\n", simcon->simparams->system);

    /* random number generator seed */
    message ("Random number generator seed: %d\n", simcon->simparams->seed);

    /* ODE parameters */
#ifdef dODE_VERSION
    message ("ODE version: %s\n", dODE_VERSION);
#endif
    message ("ODE configuration: %s\n", dGetConfiguration ());

    /* scaling factors */
    message ("scale_t: %e s\n", simparams->scale_t);
    message ("scale_m: %e g\n", simparams->scale_m);
    message ("scale_l: %e cm\n", simparams->scale_l);
    message ("scale_f: %e dyn\n", simparams->scale_f);
    message ("scale_e: %e erg\n", simparams->scale_e);

    /* environment parameters */
    message ("Ion concentration: %.2f mM\n", simcon->simparams->ion_concentration);

    /* DNA parameters */
    message ("erp: %.1f, cfm: %.1e\n", simcon->simparams->erp_glob, simcon->simparams->cfm_glob);

    /* force field parameters */
    if (simcon->simparams->forces_on) {
      message ("potential_type = %s\n", simcon->simparams->potential_type);
      if (strcmp (simcon->simparams->potential_type, "LJ")==0) {
	message ("LJ_eps = %f\n", simcon->simparams->LJ_eps);
	message ("LJ_sigma = %f\n", simcon->simparams->LJ_sigma);
      }
      else if (strcmp (simcon->simparams->potential_type, "Todd")==0) {
	message ("lambda = %f\n", simcon->simparams->lambda);
	message ("r_eq = %f\n", simcon->simparams->r_eq);
	message ("CA = %f\n", simcon->simparams->CA);
      }
      else if (strcmp (simcon->simparams->potential_type, "LJ_trunc")==0) {
	message ("LJ_eps = %f\n", simcon->simparams->LJ_eps);
	message ("LJ_sigma = %f\n", simcon->simparams->LJ_sigma);
	message ("r_trunc = %f\n", simcon->simparams->r_trunc);
      }
    }

    /* integrator parameters */
    message ("Langevin sigma = %.3f (%.3e s^-1)\n",
	simparams->sigma_langevin,
	simparams->sigma_langevin / simparams->scale_t);
    message ("ode_step: %.5f\n", simcon->simparams->ode_step);
  }
}
