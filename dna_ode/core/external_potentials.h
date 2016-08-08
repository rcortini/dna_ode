/* dna_ode: simulating DNA using the ODE physics engine
 *
 * Copyright (C) 2014, 2015, 2016  Ruggero Cortini

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

#ifndef __EXTERNAL_POTENTIALS_H__
#define __EXTERNAL_POTENTIALS_H__

#include "dna_ode_core.h"
#include "ode_body_functions.h"

struct external_potential_f {
  void (*f) (dBodyID b, void *p, t_real *u, t_real *f);
  void *p;
};
typedef struct external_potential_f external_potential_f;

t_real external_forces_add (external_potential_f *external_U, dna_ode *dna);

void external_potential (const dBodyID b, external_potential_f *U, t_real *u, t_real *f);

struct harmonic_z_parameters {
  t_real k;
  t_real z0;
};

void harmonic_z_potential (dBodyID b, void *p, t_real *u, t_real *f);

struct morse_z_parameters {
  t_real lambda;
  t_real eps;
  t_real z0;
};

void morse_z_potential (dBodyID b, void *p, t_real *u, t_real *f);

#endif
