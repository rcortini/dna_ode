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



int calculate_distance_map (trajectory_frame *frame, void *p) {
  unsigned int i, j;
  unsigned int *nbodies = (unsigned int *) p;
  distance_map *map = distance_map_alloc (*nbodies);

  /* assign distance map to frame data */
  frame->data = map;

  /* calculate the distance map */
  for (i=0; i<*nbodies; i++) {
    for (j=i+1; j<*nbodies; j++) {
      /* calculate the distance vector */
      t_vector diff;
      vector_diff (diff,
	  frame->body_i_state [i].r,
	  frame->body_i_state [j].r);

      /* set the distance map */
      set_distance_ij (i, j, vector_norm (diff), map);
    }
  }
  
  /* print distance map for each frame */
  print_distance_map (map);

  return LOAD_TRAJ_FRAME_SUCCESS;
}



int ta_distance_map (int argc, char *argv [], config_t *config) {
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
  dna_nsegments = mycfg_read_int (config, "dna_nsegments");
  system = mycfg_read_string (config, "system");
  nbodies = system_nbodies (system, dna_nsegments);

  /* open trajectory file and skip to the frame that we want to analyze */
  load_trajectory_ret_code = load_trajectory (nbodies, trajectory_file_name, &traj, calculate_distance_map, &dna_nsegments);

  /* check if success */
  if (load_trajectory_ret_code==LOAD_TRAJECTORY_FAIL) {
    err_message ("Trajectory loading failed. Exiting\n");
    return load_trajectory_ret_code;
  }

  /* free memory, close files, and exit */
  free_trajectory (&traj);
  return 0;
}
