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




int calculate_distances (trajectory_frame *frame, void *p) {
  unsigned int i, j, count;
  unsigned int *nbodies = (unsigned int *) p;
  const unsigned int N = (*nbodies)*(*nbodies-1)/2;
  t_real *distances;

  /* alloc memory and assign to frame data */
  frame->data = (t_real *) malloc (N*sizeof (t_real));
  distances = (t_real *) frame->data;

  /* get the distance map */
  count = 0;
  for (i=0; i<*nbodies; i++) {
    for (j=i+1; j<*nbodies; j++) {
      /* calculate the distance vector */
      t_vector diff;
      vector_diff (diff,
	  frame->body_i_state [i].r,
	  frame->body_i_state [j].r);

      /* set the distance norm */
      distances [count++] = vector_norm (diff);
    }
  }

  return LOAD_TRAJ_FRAME_SUCCESS;
}



int ta_distance_hist (int argc, char *argv [], config_t *config) {
  int load_trajectory_ret_code;
  unsigned int dna_nsegments, i, j, nbodies, bin_number;
  char *trajectory_file_name, *system;
  trajectory traj;
  gsl_histogram *distance_hist;
  t_real min_dist, max_dist, sum;

  /* check if user gave the right number of parameters in input */
  if (argc != 6) {
    err_message ("Incorrect number of parameters\n");
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* set the values of the parameters */
  bin_number = atoi (argv [4]);
  min_dist = 0.;
  max_dist = atof (argv [5]);

  /* read the values that we will use from the configuration file */
  trajectory_file_name = mycfg_read_string (config, "trajectory_file_name");
  dna_nsegments = mycfg_read_int (config, "dna_nsegments");
  system = mycfg_read_string (config, "system");
  nbodies = system_nbodies (system, dna_nsegments);

  /* load trajectory */
  load_trajectory_ret_code = load_trajectory (nbodies, trajectory_file_name, &traj, calculate_distances, &dna_nsegments);

  /* check if success */
  if (load_trajectory_ret_code==LOAD_TRAJECTORY_FAIL) {
    err_message ("Trajectory loading failed. Exiting\n");
    return load_trajectory_ret_code;
  }

  /* calculate the histogram of angles */
  distance_hist = gsl_histogram_alloc (bin_number);
  gsl_histogram_set_ranges_uniform (distance_hist, min_dist, max_dist);
  for (i=0; i<traj.nsteps; i++) {
    t_real *distances = (t_real *) traj.frame [i].data;
    for (j=0; j<nbodies*(nbodies-1)/2; j++) {
      gsl_histogram_increment (distance_hist, distances [j]);
    }
  }

  /* get total number of values in the histogram */
  sum = gsl_histogram_sum (distance_hist);

  /* Output file. Print header */
  for (i=0; i<bin_number; i++) {
    t_real lower, upper, mean;
    
    /* get the range of the i-th bin in the histogram, and calculate mean */
    gsl_histogram_get_range (distance_hist, i, &lower, &upper);
    mean = (upper+lower)/2.;

    /* write line in file */
    message ("%f %f\n",
	mean,
	gsl_histogram_get (distance_hist, i)/sum);
  }

  /* free memory, close files, and exit */
  free_trajectory (&traj);
  gsl_histogram_free (distance_hist);
  return 0;
}
