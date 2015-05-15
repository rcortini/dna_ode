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

#include <gsl/gsl_histogram.h>
#include "dna_ode_core.h"
#include "dna_ode.h"
#include "dna_ode-analyze.h"




int calculate_bending_angles (trajectory_frame *frame, void *p) {
  unsigned int i;
  const unsigned int nbodies = frame->nbodies;
  t_real *angles;
  (void) p;

  /* alloc memory and assign to frame data */
  frame->data = (t_real *) malloc ((nbodies-1)*sizeof (t_real));
  angles = (t_real *) frame->data;

  /* get the angles */
  for (i=0; i<(nbodies-1); i++) {
    /* calculate angle between the vectors */
    angles [i] = fabs (angle_between_vectors (frame->body_i_state [i].t, frame->body_i_state [i+1].t));
  }

  return LOAD_TRAJ_FRAME_SUCCESS;
}



int ta_angle_hist (int argc, char *argv [], config_t *config) {
  int load_trajectory_ret_code;
  unsigned int dna_nsegments, i, j, nbodies, bin_number;
  char *trajectory_file_name, *output_file_name, *system;
  FILE *output_file;
  trajectory traj;
  gsl_histogram *angle_hist;
  t_real min_angle, max_angle, sum;

  /* check if user gave the right number of parameters in input */
  if (argc != 8) {
    err_message ("Incorrect number of parameters\n");
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* set the values of the parameters */
  bin_number = atoi (argv [4]);
  min_angle = atof (argv [5]);
  max_angle = atof (argv [6]);
  output_file_name = argv [7];

  /* read the values that we will use from the configuration file */
  trajectory_file_name = mycfg_read_string (config, "trajectory_file_name");
  dna_nsegments = mycfg_read_int (config, "dna_nsegments");
  system = mycfg_read_string (config, "system");
  nbodies = system_nbodies (system, dna_nsegments);

  /* load trajectory */
  load_trajectory_ret_code = load_trajectory (nbodies, trajectory_file_name, &traj, calculate_bending_angles, NULL);

  /* check if success */
  if (load_trajectory_ret_code==LOAD_TRAJECTORY_FAIL) {
    err_message ("Trajectory loading failed. Exiting\n");
    return 1;
  }

  /* calculate the histogram of angles */
  angle_hist = gsl_histogram_alloc (bin_number);
  gsl_histogram_set_ranges_uniform (angle_hist, min_angle, max_angle);
  for (i=0; i<traj.nsteps; i++) {
    t_real *angles = (t_real *) traj.frame [i].data;
    for (j=0; j<(nbodies-1); j++) {
      gsl_histogram_increment (angle_hist, angles [j]);
    }
  }

  /* get total number of values in the histogram */
  sum = gsl_histogram_sum (angle_hist);

  /* Output file. Print header */
  output_file = my_fopen (output_file_name, "w");
  file_message (output_file, "# <bin_center>     <frequency>\n");
  for (i=0; i<bin_number; i++) {
    t_real lower, upper, mean;
    
    /* get the range of the i-th bin in the histogram, and calculate mean */
    gsl_histogram_get_range (angle_hist, i, &lower, &upper);
    mean = (upper+lower)/2.;

    /* write line in file */
    file_message (output_file, "%f %f\n",
	mean,
	gsl_histogram_get (angle_hist, i)/sum);
  }

  /* free memory, close files, and exit */
  free_trajectory (&traj);
  gsl_histogram_free (angle_hist);
  fclose (output_file);

  return 0;
}
