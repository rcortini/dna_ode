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

#include "dna_ode.h"


/* usage message */
void print_usage_getp (const char *program_name) {
  err_message ("Usage: %s <config_file> <parameter_name>\n", program_name);
}



int dna_ode_getp (int argc, char **argv, config_t *simconfig) {
  int retcode;
  char *parameter_name;

  /* check if user supplied a valid program call */
  if (argc != 4) {
    err_message ("Invalid program call\n");
    print_usage_getp (DNA_ODE_NAME);
    exit (EXIT_FAILURE);
  }

  /* get requested parameter name */
  parameter_name = argv [3];

  /* Cycle on all possible parameter type. If it finds the right parameter, print and exit.
   * Otherwise exit with error */

  /* int */
  {
    int *res = (int *) malloc (sizeof (int)); 
    retcode = config_lookup_int (simconfig, parameter_name, res);
    if (retcode == CONFIG_TRUE) {
      message ("%d\n", *res);
      free (res);
      return 0;
    }
    free (res);
  }

  /* int 64 */
  {
    long long *res = (long long *) malloc (sizeof (long long)); 
    retcode = config_lookup_int64 (simconfig, parameter_name, res);
    if (retcode == CONFIG_TRUE) {
      printf ("%lld\n", *res);
      free (res);
      return 0;
    }
    free (res);
  }

  /* double */
  {
    double *res = (double *) malloc (sizeof (double)); 
    retcode = config_lookup_float (simconfig, parameter_name, res);
    if (retcode == CONFIG_TRUE) {
      printf ("%.3f\n", *res);
      free (res);
      return 0;
    }
    free (res);
  }

  /* string */
  {
    const char *res;
    retcode = config_lookup_string (simconfig, parameter_name, &res);
    if (retcode == CONFIG_TRUE) {
      printf ("%s\n", res);
      return 0;
    }
  }

  /* bool */
  {
    int *res = (int *) malloc (sizeof (int)); 
    retcode = config_lookup_bool (simconfig, parameter_name, res);
    if (retcode == CONFIG_TRUE) {
      printf ("%d\n", *res);
      free (res);
      return 0;
    }
    free (res);
  }

  /* if control arrives here, CONFIG_TRUE was never returned. It means that we could not read
   * specified parameter */
  err_message ("Could not read parameter %s\n", parameter_name);
  return 1;
}
