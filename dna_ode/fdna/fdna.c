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
#include "fdna.h"

int fdna_main (simcontext *simcon) {
  unsigned int initial_configuration_file_exists;

  /* initialize system-specific simulation context */
  fdna_simcontext *fdna_simcon = (fdna_simcontext *) malloc (sizeof (*fdna_simcon));
  fdna_init_simulation_context (simcon, fdna_simcon);
 
  /*
   * ASSIGN THE INITIAL CONFIGURATION OF THE DNA AND THE BEAD. CHECK IF USER PROVIDED A FILE FOR THE INITIAL CONFIGURATION,
   * AND IN SUCH CASE LOAD THE FILE
   */

  /* check if file name was given and it exists */
  initial_configuration_file_exists = 0;
  if (simcon->simparams->initial_configuration_file_name [0] != '\0') {
    simcon->initial_configuration_file = my_fopen (simcon->simparams->initial_configuration_file_name, "r");
    if (simcon->initial_configuration_file != NULL)
      initial_configuration_file_exists = 1;
  }

  /* if the configuration file exists, load configuration, otherwise create a new one */
  if (initial_configuration_file_exists) {
    message ("Loading initial configuration from %s\n", simcon->simparams->initial_configuration_file_name);
    fdna_load_initial_configuration_from_file (simcon, fdna_simcon);
    fdna_joints_create (simcon, fdna_simcon->dna);
  }
  else {
    /* in this case we create a random new initial configuration, and the clock starts from 0 */
    simcon->sim_clock_0 = 0.;
    fdna_joints_create (simcon, fdna_simcon->dna);
    fdna_init_molecule_positions_orientations (simcon, fdna_simcon);
    fdna_init_molecule_velocities (simcon, fdna_simcon);

    /* writes initial configuration to output file */
    fdna_output_data (simcon, fdna_simcon);
  }

  /*
   * SIMULATION
   */

  /* start the clock: sim_clock_0 was initialized either in "load_initial_configuration_from_file" or at the end
   * of the generation of the random initial configuration (see above). elapsed_time is always initialized to 0 at
   * the end of "init_simulation_context". */
  simcon->sim_clock = simcon->sim_clock_0 + simcon->elapsed_time;

  /* start simulation loop: lasts until elapsed_time is equal to total_time */
  while (simcon->elapsed_time < simcon->simparams->total_time) {
    /* performs a single simulation step */
    fdna_simstep (simcon, fdna_simcon);

    /* data output is written if the current step is an integer multiple of output_interval */
    if (simcon->sim_tick % simcon->simparams->output_interval == 0) {
      fdna_output_data (simcon, fdna_simcon);
    }

    /* print output is written if the current step is an integer multiple of print_interval */
    if (simcon->sim_tick % simcon->simparams->print_interval == 0) {
      fdna_print_data (simcon, fdna_simcon);
      /* print_joints_data (simcon); */
    }
  }

  /*
   * FINALIZATION
   */

  /* free memory */
  fdna_free_memory (fdna_simcon);

  /* exit */
  return 0;
}
