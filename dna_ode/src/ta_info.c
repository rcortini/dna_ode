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



int ta_info (int argc, char *argv [], config_t *config) {
  int load_trajectory_ret_code;
  unsigned int dna_nsegments, nbodies;
  char *trajectory_file_name, *system;
  trajectory traj;
  (void) argv;

  /* check if user gave the right number of parameters in input */
  if (argc != 4) {
    err_message ("Incorrect number of parameters\n");
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* read the values that we will use from the configuration file */
  trajectory_file_name = mycfg_read_string (config, "trajectory_file_name");
  system = mycfg_read_string (config, "system");
  dna_nsegments = mycfg_read_int (config, "dna_nsegments");
  nbodies = system_nbodies (system, dna_nsegments);

  /* load trajectory */
  load_trajectory_ret_code = load_trajectory (nbodies, trajectory_file_name, &traj, NULL, NULL);

  /* check if success */
  if (load_trajectory_ret_code==LOAD_TRAJECTORY_FAIL) {
    err_message ("Trajectory loading failed. Exiting\n");
    return load_trajectory_ret_code;
  }

  /* print info */
  message ("Trajectory file name: %s\n", trajectory_file_name);
  message ("nbodies: %d\n", nbodies);
  message ("nsteps: %d\n", traj.nsteps);

  /* free memory, close files, and exit */
  free_trajectory (&traj);
  return 0;
}
