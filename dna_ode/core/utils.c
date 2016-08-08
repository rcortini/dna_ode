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

#include "utils.h"

FILE *my_fopen (const char *path, const char *mode) {
  FILE *fp = fopen (path, mode);
  if (fp == NULL) {
    err_message ("Could not open %s\n", path);
    exit (EXIT_FAILURE);
  }
  else
    return fp;
}



/* a function to display a message */
void message (char *text, ...) {
  va_list args;
  va_start (args, text);
  vfprintf (stdout, text, args);
  va_end (args);
  fflush (stdout);
}



/* a function to print to file */
void file_message (FILE *file, char *text, ...) {
  va_list args;
  va_start (args, text);
  vfprintf (file, text, args);
  va_end (args);
  fflush (file);
}



/* a function to display an error message */
void err_message (char *text, ...) {
  va_list args;
  va_start (args, text);
  vfprintf (stderr, text, args);
  va_end (args);
  fflush (stderr);
}


/* prints a vector3d */
void print_vector3d (const t_real *vec) {
  message ("%.16e %.16e %.16e", vec[0], vec[1], vec[2]);
}



/* prints a matrix */
void print_matrix (const t_real *m) {
  unsigned int i, j;
  for (i=0; i<3; i++) {
    for (j=0; j<4; j++) {
      message ("%.16e ", m[4*i+j]);
    }
    message ("\n");
  }
}



/* prints a gsl matrix */
void print_gsl_matrix (unsigned int nrows, unsigned int ncolumns, gsl_matrix *m) {
  unsigned int i, j;
  for (i=0; i<nrows; i++) {
    for (j=0; j<ncolumns; j++) {
      message ("%.3f ", gsl_matrix_get (m, i, j));
    }
    message ("\n");
  }
}



/* prints a gsl vector */
void print_gsl_vector (unsigned int nrows, gsl_vector *m) {
  unsigned int i;
  for (i=0; i<nrows; i++) {
    message ("%.3f ", gsl_vector_get (m, i));
  }
  message ("\n");
}



/* returns the system identity based on the string "system" */
system_id which_system (char *system) {
  if (strcmp (system, "fdna")==0)
    return FDNA;
  if (strcmp (system, "sdna")==0)
    return SDNA;
  if (strcmp (system, "cdna")==0)
    return CDNA;
  else {
    err_message ("Invalid system type: %s\n", system);
    exit (EXIT_FAILURE);
  }
}



/* returns the number of bodies in a system depending on what system it is */
unsigned int system_nbodies (char *system, unsigned int dna_nsegments) {
  system_id sysid = which_system (system);
  switch (sysid) {
    case (FDNA) :
      return dna_nsegments;
      break;
    case (SDNA) :
      return dna_nsegments + 1;
      break;
    case (CDNA) :
      return dna_nsegments;
      break;
    default :
      err_message ("Invalid system type: %s\n", system);
      exit (EXIT_FAILURE);
      break;
  };
}

unsigned int system_njoints (char *system, unsigned int dna_nsegments) {
  system_id sysid = which_system (system);
  switch (sysid) {
    case (FDNA) :
      return dna_nsegments-1;
      break;
    case (SDNA) :
      return dna_nsegments + 1;
      break;
    case (CDNA) :
      return dna_nsegments;
      break;
    default :
      err_message ("Invalid system type: %s\n", system);
      exit (EXIT_FAILURE);
      break;
  }
}

/* Metropolis Monte Carlo method */
int metropolis (t_real beta, t_real E_initial, t_real E_final) {
  if (E_final <= E_initial)
    return 1;
  else {
    if (ran_uniform (0.,1.) < exp (beta*(E_initial-E_final)))
      return 1;
    else
      return 0;
  }
}
