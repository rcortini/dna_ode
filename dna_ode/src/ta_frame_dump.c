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




int ta_frame_dump (int argc, char *argv [], config_t *config) {
  int load_trajectory_ret_code, load_trajectory_frame_ret_code;
  unsigned int dna_nsegments, frame_number, nbodies, i;
  char *trajectory_file_name, *system;
  trajectory_frame frame;
  FILE *trajectory_file;

  /* check if user gave the right number of parameters in input */
  if (argc != 5) {
    err_message ("Incorrect number of parameters\n");
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* set the values of the parameters */
  frame_number = atoi (argv [4]);

  /* read the values that we will use from the configuration file */
  trajectory_file_name = mycfg_read_string (config, "trajectory_file_name");
  dna_nsegments = mycfg_read_int (config, "dna_nsegments");
  system = mycfg_read_string (config, "system");
  nbodies = system_nbodies (system, dna_nsegments);

  /* open trajectory file */
  trajectory_file = my_fopen (trajectory_file_name, "r");

  /* place the file pointer to the frame_number-th frame of the trajectory file */
  load_trajectory_ret_code = skip_to_frame_n (frame_number, nbodies, trajectory_file);
  if (load_trajectory_ret_code==LOAD_TRAJECTORY_FAIL) {
     err_message ("Trajectory loading failed. Exiting\n");
     return load_trajectory_ret_code;
  }

  /* now that our file pointer is at the correct point, load the trajectory frame */
  load_trajectory_frame_ret_code = load_trajectory_frame (nbodies, trajectory_file, &frame, NULL, NULL);
  if (load_trajectory_frame_ret_code==LOAD_TRAJ_FRAME_FAIL) {
     err_message ("Trajectory loading failed. Exiting\n");
     return load_trajectory_frame_ret_code;
  }

  /* write data */
  for (i=0; i<nbodies; i++) {
    print_body_state (stdout, &frame, i);
  }

  /* free memory */
  fclose (trajectory_file);
  return 0;
}
