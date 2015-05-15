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
#include "dna_ode-run.h"

void print_usage_run (const char *program_name) {
  err_message ("Usage: %s run <conf_file> <seed>\n", program_name);
}

int dna_ode_run (int argc, char **argv, config_t *simconfig) {
  int retcode;
  simparameters *simparams;
  simcontext *simcon;

  /*
   * INITIALIZATION
   */

  /* check that seed was given */
  if (argc < 4) {
    print_usage_run (argv [0]);
    exit (EXIT_FAILURE);
  }

  /* memory allocation for the simulation parameters, statistics, and context */
  simparams = (simparameters *) malloc (sizeof (simparameters));
  simcon = (simcontext *) malloc (sizeof (simcontext));

  /* get RNG seed */
  simparams->seed = atoi (argv [3]);

  /* initialization of ODE */
  dInitODE ();

  /*
   * INITIALIZE SIMULATION PARAMETERS, RANDOM NUMBER GENERATION, AND CONTEXT
   */
  init_simulation_parameters (simparams, simconfig);
  my_init_RNG (simparams->seed);
  init_simulation_context (simcon, simparams);

  /*
   * NOW WE SWITCH TO THE SYSTEM CHOSEN
   */
  switch (simparams->sysid) {
    case SDNA :
      retcode = sdna_main (simcon);
      break;
    default :
      err_message ("Invalid system ID\n");
      retcode = 1;
  }

  /* final stuff to terminal */
  finish_output (simcon);

  /* close output files */
  close_output_files (simcon);

  /* free memory */
  free_memory (simcon);

  /* close ODE */
  dCloseODE();

  return retcode;
}
