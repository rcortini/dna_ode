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

#include "config_read.h"
#include "dna_ode_core.h"
#include "dna_ode-run.h"
#include "simparameters.h"
 


/* reads all the necessary parameters from the user-specified configuration file:
 * Reads the following common parameters:
 * 1. ODE
 * 2. environment
 * 3. DNA
 * 4. thermostat
 * 5. inter-DNA forces
 * 6. external forces
 * 7. file names
 * 8. time-related parameters
 */
void init_simulation_parameters (simparameters *simparams, config_t *simconfig) {
  /***********************************************
   *                                             *
   * PHASE 1: INITIALIZATION                     *
   *                                             *
   * ********************************************/

  /* init simconfig */
  simparams->simconfig = simconfig;

  /***********************************************
   *                                             *
   * PHASE 2: READ PARAMETERS FROM CONFIG FILE   *
   *                                             *
   * ********************************************/

  /* get system name and identity*/
  simparams->system = mycfg_read_string (simconfig, "system");
  simparams->sysid = which_system (simparams->system);

  /*
   * ODE PARAMETERS
   */
  simparams->erp_glob = mycfg_read_double (simconfig, "erp_glob");
  simparams->cfm_glob = mycfg_read_double (simconfig, "cfm_glob");
  simparams->sor_numiter = mycfg_read_int (simconfig, "sor_numiter");
  simparams->sor_omega = mycfg_read_double (simconfig, "sor_omega");
  simparams->finite_rotation = mycfg_read_int (simconfig, "finite_rotation");
  simparams->gyroscopic_mode = mycfg_read_int (simconfig, "gyroscopic_mode");
  simparams->contact_mu = mycfg_read_double (simconfig, "contact_mu");
  simparams->contact_mu2 = mycfg_read_double (simconfig, "contact_mu2");
  simparams->contact_bounce = mycfg_read_double (simconfig, "contact_bounce");
  simparams->contact_bounce_vel = mycfg_read_double (simconfig, "contact_bounce_vel");

  /*
   * ENVIRONMENT PARAMETERS
   */
  simparams->water_dielectric = mycfg_read_double (simconfig, "water_dielectric");
  simparams->temperature = mycfg_read_double (simconfig, "temperature");
  simparams->ion_concentration = mycfg_read_double (simconfig, "ion_concentration");
  simparams->l_bjerrum = E_CHARGE*E_CHARGE/(simparams->water_dielectric * K_B * simparams->temperature);
  simparams->beta = 1./(K_B * simparams->temperature);

  /*
   * DNA PARAMETERS
   */
  simparams->dna_rad = mycfg_read_double (simconfig, "dna_rad") * 1e-8; /* convert to cm */
  simparams->dna_axialrise = mycfg_read_double (simconfig, "dna_axialrise") * 1e-8; /* convert to cm */
  simparams->dna_nsegments = mycfg_read_int (simconfig, "dna_nsegments");
  simparams->dna_nbasepairs_in_segment = mycfg_read_int (simconfig, "dna_nbasepairs_in_segment");
  simparams->dna_basepair_mass = mycfg_read_double (simconfig, "dna_basepair_mass");
  simparams->dna_basepairs_per_turn = mycfg_read_double (simconfig, "dna_basepairs_per_turn");
  simparams->dna_lbp = mycfg_read_double (simconfig, "dna_lbp") * 1e-8; /* convert to cm */
  simparams->dna_ltp = mycfg_read_double (simconfig, "dna_ltp") * 1e-8; /* convert to cm */
  simparams->dna_kuhn_length = simparams->dna_lbp * 2.;
  simparams->dna_nbasepairs = simparams->dna_nsegments * simparams->dna_nbasepairs_in_segment;
  simparams->dna_segment_length = simparams->dna_axialrise * simparams->dna_nbasepairs_in_segment;
  simparams->dna_length = simparams->dna_nsegments * simparams->dna_segment_length;
  simparams->dna_effective_rad = simparams->dna_rad +
    1./sqrt (8*M_PI*simparams->l_bjerrum * N_AVO * simparams->ion_concentration*1e-3); /* r_eff = r_DNA + l_Debye */
  simparams->dna_segment_mass = simparams->dna_basepair_mass * simparams->dna_nbasepairs_in_segment;
  simparams->tol = mycfg_read_double (simconfig, "tol");

  /*
   * THERMOSTAT PARAMETERS
   */
  simparams->thermostat_type = mycfg_read_string (simconfig, "thermostat_type");
  simparams->sigma_langevin = mycfg_read_double (simconfig, "sigma_langevin");

  /*
   * INTER-DNA FORCES PARAMETERS
   */
  simparams->forces_on = mycfg_read_int (simconfig, "forces_on");
  if (simparams->forces_on) {
    simparams->potential_type = mycfg_read_string (simconfig, "potential_type");
    if (strcmp (simparams->potential_type, "LJ")==0) {
      simparams->LJ_sigma = mycfg_read_double (simconfig, "LJ_sigma") * 1.e-8; /* convert from Angstroms to cm */
      simparams->LJ_eps = mycfg_read_double (simconfig, "LJ_eps") / 1.e-7; /* convert from kT/nm to kT/cm */
    }
    else if (strcmp (simparams->potential_type, "Todd")==0) {
      simparams->lambda = mycfg_read_double (simconfig, "lambda") * 1.e-8; /* convert from Angstroms to cm */
      simparams->CA = mycfg_read_double (simconfig, "CA") * 1.e7; /* convert from pN/nm^2 to dyn/cm^2 */
      simparams->r_eq = mycfg_read_double (simconfig, "r_eq") * 1.e-8; /* convert from Angstroms to cm */
    }
    else if (strcmp (simparams->potential_type, "LJ_trunc")==0) {
      simparams->LJ_sigma = mycfg_read_double (simconfig, "LJ_sigma") * 1.e-8; /* convert from Angstroms to cm */
      simparams->LJ_eps = mycfg_read_double (simconfig, "LJ_eps") / 1.e-7; /* convert from kT/nm to kT/cm */
      simparams->r_trunc = mycfg_read_double (simconfig, "r_trunc");
    }
    else {
      err_message ("Invalid potential function %s\n", simparams->potential_type);
      exit (EXIT_FAILURE);
    }
  }

  /*
   * EXTERNAL FORCES PARAMETERS
   */
  simparams->external_forces_on = mycfg_read_int (simconfig, "external_forces_on");
  if (simparams->external_forces_on) {
    simparams->external_potential_type = mycfg_read_string (simconfig, "external_potential_type");
    if (strcmp (simparams->external_potential_type, "harmonic_z")==0) {
      simparams->k_harmonic = mycfg_read_double (simconfig, "k_harmonic");
      simparams->z0_harmonic = mycfg_read_double (simconfig, "z0_harmonic");
    }
    else if (strcmp (simparams->external_potential_type, "morse_z")==0) {
      simparams->eps_morse = mycfg_read_double (simconfig, "eps_morse");
      simparams->z0_morse = mycfg_read_double (simconfig, "z0_morse");
      simparams->lambda_morse = mycfg_read_double (simconfig, "lambda_morse");
    }
    else if (strcmp (simparams->external_potential_type, "gravity")==0) {
      simparams->gravity_val = mycfg_read_double (simconfig, "gravity_val");
      simparams->gravity_z0 = mycfg_read_double (simconfig, "gravity_z0");
    }
    else if (strcmp (simparams->external_potential_type, "morse_z_surface")==0) {
      simparams->eps_morse = mycfg_read_double (simconfig, "eps_morse");
      simparams->z0_morse = mycfg_read_double (simconfig, "z0_morse");
      simparams->z0_surface = mycfg_read_double (simconfig, "z0_surface");
      simparams->lambda_morse = mycfg_read_double (simconfig, "lambda_morse");
    }
    else {
      err_message ("Invalid external potential function %s\n", simparams->external_potential_type);
      exit (EXIT_FAILURE);
    }
  }

  /*
   * MECHANICAL MODEL PARAMETERS
   */
  simparams->mechanical_model = mycfg_read_string (simconfig, "mechanical_model");
  if (strcmp (simparams->mechanical_model, "harmonic4")==0 ||
      strcmp (simparams->mechanical_model, "nicked_harmonic4")==0)
    simparams->k4 = mycfg_read_double (simconfig, "k4");

  /*
   * FILE NAMES
   */
  simparams->initial_configuration_file_name = mycfg_read_string (simconfig, "initial_configuration_file_name");
  simparams->restart_file_name = mycfg_read_string (simconfig, "restart_file_name");
  simparams->trajectory_file_name = mycfg_read_string (simconfig, "trajectory_file_name");
  simparams->data_file_name = mycfg_read_string (simconfig, "data_file_name");
  simparams->odepdb_file_name = mycfg_read_string (simconfig, "odepdb");

  /*
   * TIME-RELATED PARAMETERS
   */
  simparams->ode_step = mycfg_read_double (simconfig, "ode_step");
  simparams->total_time = mycfg_read_double (simconfig, "total_time");
  simparams->output_interval = mycfg_read_int64 (simconfig, "output_interval");
  simparams->print_interval = mycfg_read_int64 (simconfig, "print_interval");


  /***********************************************
   *                                             *
   * PHASE 3: SCALING                            *
   *                                             *
   * ********************************************/


  /* scaling factors for mass, distance, energy, force, and time */
  simparams->scale_m = simparams->dna_segment_mass;  /* mass */
  simparams->scale_l = simparams->dna_segment_length;  /* distance */
  simparams->scale_e = K_B*simparams->temperature;  /* energy */
  simparams->scale_f = simparams->scale_e/simparams->scale_l;  /* force */
  simparams->scale_t = sqrt (simparams->scale_m * simparams->scale_l * simparams->scale_l/simparams->scale_e);  /* time */

  /* distances scaling */
  simparams->dna_segment_length /= simparams->scale_l;
  simparams->dna_kuhn_length /= simparams->scale_l;
  simparams->dna_length /= simparams->scale_l;
  simparams->dna_effective_rad /= simparams->scale_l;
  simparams->dna_lbp /= simparams->scale_l;
  simparams->dna_ltp /= simparams->scale_l;

  /* scaling of the parameters associated with inter-DNA interaction potentials */
  if (simparams->forces_on) {
    if (strcmp (simparams->potential_type, "LJ")==0) {
      simparams->LJ_sigma /= simparams->scale_l;

      /* LJ_eps is in kT/nm so we scale only by length */
      simparams->LJ_eps *= simparams->scale_l;
    }
    if (strcmp (simparams->potential_type, "Todd")==0) {
      simparams->lambda /= simparams->scale_l;
      simparams->r_eq /= simparams->scale_l;

      /* CA is in pN/nm^2 so we scale twice by length and once by force */
      simparams->CA *= simparams->scale_l;
      simparams->CA *= simparams->scale_l;
      simparams->CA /= simparams->scale_f;
    }
    if (strcmp (simparams->potential_type, "LJ_trunc")==0) {
      simparams->LJ_sigma /= simparams->scale_l;

      /* LJ_eps is in kT/nm so we scale only by length */
      simparams->LJ_eps *= simparams->scale_l;
    }
  }

  /* masses scaling */
  simparams->dna_segment_mass /= simparams->scale_m;

  /* energy scaling */
  simparams->beta *= simparams->scale_e;

  /* force scaling */

  /* time scaling */
  simparams->sigma_langevin *= simparams->scale_t;
}
