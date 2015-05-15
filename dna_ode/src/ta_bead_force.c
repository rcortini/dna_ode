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
#include "dna_ode.h"
#include "dna_ode-analyze.h"




/* calculates force acting on bead for a given frame */
int calculate_bead_force (trajectory_frame *frame, void *p) {
  t_real *par = (t_real *) p;
  t_real z, z0, trap_k, bead_force, scale_l;
  const unsigned int nbodies = frame->nbodies;
  t_real *frame_data;

  /* calculate parameters */
  scale_l = par [0];
  z0 = par [1];
  trap_k = par [2];
  z = frame->body_i_state [nbodies-1].r[2] * scale_l;

  /* calculate bead force */
  bead_force = -trap_k * (z-z0);

  /* assign calculated values to frame data */
  frame_data = (t_real *) malloc (sizeof (t_real));
  frame_data [0] = bead_force;
  frame->data = frame_data;

  /* print data */
  message ("%f %f\n", frame->sim_clock, bead_force);

  return LOAD_TRAJ_FRAME_SUCCESS;
}



int ta_bead_force (int argc, char *argv [], config_t *config) {
  int load_trajectory_ret_code;
  unsigned int dna_nsegments, nbodies, dna_nbasepairs_in_segment;
  char *trajectory_file_name, *system;
  trajectory traj;
  t_real dna_length, bead_radius, dna_axialrise, trap_k, scale_l, rho;
  t_real par [3];
  (void) argv;

  /* check if user gave the right number of parameters in input */
  if (argc != 4) {
    err_message ("Incorrect number of parameters\n");
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* read the values that we will use from the configuration file */
  system = mycfg_read_string (config, "system");
  trajectory_file_name = mycfg_read_string (config, "trajectory_file_name");
  dna_nsegments = mycfg_read_int (config, "dna_nsegments");
  dna_axialrise = mycfg_read_double (config, "dna_axialrise")/10.; /* convert to nanometers */
  dna_nbasepairs_in_segment = mycfg_read_int (config, "dna_nbasepairs_in_segment");
  dna_length = dna_nsegments * dna_axialrise * dna_nbasepairs_in_segment;
  bead_radius = dna_length/(2*M_PI);
  trap_k = mycfg_read_double (config, "trap_k");
  rho = mycfg_read_double (config, "rho");
  scale_l = dna_nbasepairs_in_segment * dna_axialrise;

  /* define nbodies and parameters */
  nbodies = system_nbodies (system, dna_nsegments);
  par [0] = scale_l;
  par [1] = dna_length * rho + bead_radius; /* this is z0 */
  par [2] = trap_k;

  /* load trajectory */
  load_trajectory_ret_code = load_trajectory (nbodies, trajectory_file_name, &traj, calculate_bead_force, par);

  /* check if success */
  if (load_trajectory_ret_code==LOAD_TRAJECTORY_FAIL) {
    err_message ("Trajectory loading failed. Exiting\n");
    return 1;
  }

  /* free memory, close files, and exit */
  free_trajectory (&traj);
  return 0;
}
