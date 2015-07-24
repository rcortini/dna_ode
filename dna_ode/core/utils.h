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

#ifndef __UTILS_H__
#define __UTILS_H__

#include "dna_ode_core.h"
#include "random.h"

FILE *my_fopen (const char *path, const char *mode);

void message (char *text, ...);

void file_message (FILE *file, char *text, ...);

void err_message (char *text, ...);

void print_vector3d (const t_real *vec);

void print_matrix (const t_real *m);

void print_gsl_matrix (unsigned int nrows, unsigned int ncolumns, gsl_matrix *m);

void print_gsl_vector (unsigned int nrows, gsl_vector *m);

system_id which_system (char *system);

unsigned int system_nbodies (char *system, unsigned int dna_nsegments);

unsigned int system_njoints (char *system, unsigned int dna_nsegments);

int metropolis (t_real beta, t_real E_initial, t_real E_final);

#endif
