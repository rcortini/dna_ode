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

#include <stdio.h>
#include <stdlib.h>
#include "config_read.h"



/*
 *
 * CONFIGURATION FILE PARSING FUNCTIONS
 * Author: Bertrand Care
 * email: care@lptl.jussieu.fr
 *
 */

/*
  generic function for reading key=value pairs from a file
  using libconfig.
  It's "safe" because the program exits if a parameter could not be read.

  In the file, the parameters are stored under the form :
  """
  paramname = res;
  """
  paramtype specifies the type of res
  0 : int 
  1 : int 64
  2 : double
  3 : string
  4 : bool
  
  Convenience functions are available below for easier use (one for each parameter type)
  These functions call mycfg_safe_read with the right paramtype.
 */
void * mycfg_safe_read (config_t * conf, char * paramname, unsigned int paramtype) {
  int retcode;
  void * res;

  switch (paramtype) {
    case 0:
      /* int */
      res = (int *) malloc (sizeof (int)); 
      retcode = config_lookup_int (conf, paramname, res);
      break;
    case 1:
      /* int 64 */
      res = (long long *) malloc (sizeof (long long)); 
      retcode = config_lookup_int64 (conf, paramname, res);
      break;
    case 2:
      /* double */
      res = (double *) malloc (sizeof (double)); 
      retcode = config_lookup_float (conf, paramname, res);
      break;
    case 3:
      /* string */
      res = (char **) malloc (sizeof (char **)); 
      retcode = config_lookup_string (conf, paramname, res);
      break;
    case 4:
      /* bool */
      res = (unsigned int *) malloc (sizeof (unsigned int)); 
      retcode = config_lookup_bool (conf, paramname, res);
      break;
    default : 
      fprintf (stderr, "Error : wrong parameter type specified\n");
      exit (EXIT_FAILURE);
      break;
    }
  
  /* check if the return code of the config_lookup functions is CONFIG_TRUE */
  if (retcode != CONFIG_TRUE) {
      fprintf (stderr, "Error (%d): could not read setting %s\n", retcode, paramname);
      exit (EXIT_FAILURE);
  }

  /* return the pointer to the result */
  return res;
}

/*
  ------------------------------
  convenience functions for loading parameters
  they call mycfg_safe_read with the right paramtype value
  and cast the result before returning it.
  Use them in the big parameter-loading function below.
  types :
  - int : 0
  - int64 : 1
  - double : 2
  - string : 3
  - bool : 4
  ------------------------------
*/
int mycfg_read_int (config_t * conf, char * paramname) {
  int * par_pointer = mycfg_safe_read (conf, paramname, 0);
  int par = *par_pointer;
  free (par_pointer);
  return par;
}

long mycfg_read_int64 (config_t * conf, char * paramname) {
  long long * par_pointer = mycfg_safe_read (conf, paramname, 1);
  long long par = *par_pointer;
  free (par_pointer);
  return par;
}

double mycfg_read_double(config_t * conf, char * paramname) {
  double * par_pointer = mycfg_safe_read (conf, paramname, 2);
  double par = *par_pointer;
  free (par_pointer);
  return par;
}

char * mycfg_read_string (config_t * conf, char * paramname) {
  char ** par_pointer = mycfg_safe_read (conf, paramname, 3);
  char *par = *par_pointer;
  free (par_pointer);
  return par;
}

unsigned int mycfg_read_bool (config_t * conf, char * paramname) {
  unsigned int * par_pointer = mycfg_safe_read (conf, paramname, 4);
  unsigned int par = *par_pointer;
  free (par_pointer);
  return par;
}
