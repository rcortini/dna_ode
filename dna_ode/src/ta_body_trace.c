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




int ta_body_trace (int argc, char *argv [], config_t *config) {
  int load_trajectory_ret_code;
  unsigned int i, dna_nsegments, segment_to_trace, nbodies;
  char *trajectory_file_name, *system;
  FILE *trajectory_file;
  trajectory traj;

  /* check if user gave the right number of parameters in input */
  if (argc != 5) {
    err_message ("Incorrect number of parameters\n");
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* get parameters */
  segment_to_trace = atoi (argv [4]);

  /* read the values that we will use from the configuration file */
  trajectory_file_name = mycfg_read_string (config, "trajectory_file_name");
  dna_nsegments = mycfg_read_int (config, "dna_nsegments");
  system = mycfg_read_string (config, "system");
  nbodies = system_nbodies (system, dna_nsegments);

  /* check if the segment we want to trace is a segment that is valid */
  if (segment_to_trace>nbodies-1) {
    err_message ("Invalid body number to trace (must be smaller than %d)\n", nbodies-1);
    return 1;
  }

  /* open trajectory file */
  trajectory_file = my_fopen (trajectory_file_name, "r");

  /* check if trajectory file exists */
  if (trajectory_file == NULL) {
    err_message ("Trajectory file %s not found\n", trajectory_file_name);
    return 1;
  }

  /* load trajectory */
  load_trajectory_ret_code = load_trajectory (nbodies, trajectory_file_name, &traj, NULL, NULL);

  /* check if success */
  if (load_trajectory_ret_code==LOAD_TRAJECTORY_FAIL) {
    err_message ("Trajectory loading failed. Exiting\n");
    return load_trajectory_ret_code;
  }

  /* now print the trace of the requested segment */
  for (i=0; i<traj.nsteps; i++) {
    print_body_state (stdout, &traj.frame [i], segment_to_trace);
  }

  /* free memory */
  fclose (trajectory_file);
  return 0;
}
