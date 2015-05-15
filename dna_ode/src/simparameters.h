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

#ifndef __SIMPARAMETERS_H__
#define __SIMPARAMETERS_H__


#include "config_read.h"
#include "dna_ode_core.h"


/* this structure contains all the simulation parameters */
struct simparameters {
  char *system;
  system_id sysid;

  /* configuration object */
  config_t *simconfig;

  /* random number generator seed */
  unsigned long seed;

  /* ODE world params */
  t_real erp_glob;
  t_real cfm_glob;

  /* ODE dynamics params */
  unsigned int finite_rotation;
  unsigned int gyroscopic_mode;

  /* environment parameters */
  t_real water_dielectric;
  t_real temperature;
  t_real beta;
  t_real l_bjerrum;
  t_real ion_concentration;

  /* DNA parameters */
  unsigned int dna_nbasepairs;
  unsigned int dna_nbasepairs_in_segment;
  unsigned int dna_nsegments;
  t_real dna_rad;
  t_real dna_axialrise;
  t_real dna_length;
  t_real dna_effective_rad;
  t_real dna_basepair_mass;
  t_real dna_segment_mass;
  t_real dna_segment_length;
  t_real dna_basepairs_per_turn;
  t_real dna_lbp; /* bending persistence length */
  t_real dna_ltp; /* twisting persistence length */
  t_real dna_kuhn_length;

  /* contact point parameters */
  t_real contact_mu;
  t_real contact_mu2;
  t_real contact_bounce;
  t_real contact_bounce_vel;

  /* parameters for the calculation of twisting and bending energies */
  t_real tol;
  t_real mcos;
  t_real gb;
  t_real gt;

  /* inter-segment DNA forces parameters */
  unsigned int forces_on;
  char * potential_type;

  /* for Todd potential */
  t_real lambda;
  t_real CA;
  t_real r_eq;

  /* for LJ potential */
  t_real LJ_sigma;
  t_real LJ_eps;
  t_real r_trunc; /* this is for LJ_trunc */

  /* inter-segment DNA forces parameters */
  unsigned int external_forces_on;
  char * external_potential_type;

  /* harmonic_z external potential parameters */
  t_real k_harmonic;
  t_real z0_harmonic;

  /* morse_z external potential parameters */
  t_real eps_morse;
  t_real z0_morse;
  t_real lambda_morse;

  /* gravity external potential parameters */
  t_real gravity_val;
  t_real gravity_z0;

  /* morse_surface external potential parameters */
  t_real z0_surface;

  /* Langevin dynamics parameters */
  char * thermostat_type;
  t_real sigma_langevin;
  t_real mE;

  /* mechanical model */
  char * mechanical_model;
  t_real k4;

  /* ODE "successive over-relaxation params */
  int sor_numiter;
  t_real sor_omega;

  /* time-related stuff */
  t_real ode_step;
  unsigned long output_interval;
  unsigned long print_interval;
  t_real total_time;

  /* initial configuration file */
  char * restart_file_name;
  char * initial_configuration_file_name;
  char * trajectory_file_name;
  char * data_file_name;
  char * odepdb_file_name;

  /* scaling parameters */
  t_real scale_t; /* time */
  t_real scale_e; /* energy */
  t_real scale_l; /* distance */
  t_real scale_m; /* mass */
  t_real scale_f; /* force */
};
typedef struct simparameters simparameters;



/* FUNCTION PROTOTYPES */

void init_simulation_parameters (simparameters *simparams, config_t *simconfig);

#endif
