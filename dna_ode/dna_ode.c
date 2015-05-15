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

void print_usage (const char *program_name) {
  err_message ("Usage: %s <command> <conf_file> <options>\n", program_name);
  err_message ("          command = run, analyze, getp\n");
}

void free_simconfig (config_t *simconfig) {
  config_destroy (simconfig);
  free (simconfig);
}

int main (int argc, char **argv) {
  unsigned int config_file_retcode, prog_retcode;
  config_t * simconfig;
  const char *command = argv [1];
  const char *config_file_name = argv [2];

  /* verify that at least two arguments were given */
  if (argc<3) {
    err_message ("Incorrect usage\n");
    print_usage (DNA_ODE_NAME);
    exit (EXIT_FAILURE);
  }

  /* initialize the configuration file reader (from libconfig) */
  simconfig = (config_t *) malloc (sizeof (config_t));
  config_init (simconfig);

  /* read config file */
  config_file_retcode = config_read_file (simconfig, config_file_name);
  if (config_file_retcode != CONFIG_TRUE) {
    err_message ("Error (%d): could not read config file %s\n", config_file_retcode, config_file_name);
    free_simconfig (simconfig);
    exit (EXIT_FAILURE);
  }

  /* now execute the command */
  if (strcmp (command, "run")==0)
    prog_retcode = dna_ode_run (argc, argv, simconfig);
  else if (strcmp (command, "analyze")==0)
    prog_retcode = dna_ode_analyze (argc, argv, simconfig);
  else if (strcmp (command, "getp")==0)
    prog_retcode = dna_ode_getp (argc, argv, simconfig);
  else {
    err_message ("Unrecognized command %s\n", command);
    print_usage (DNA_ODE_NAME);
    exit (EXIT_FAILURE);
  }

  /* free the memory and return */
  free_simconfig (simconfig);

  return prog_retcode;
}
