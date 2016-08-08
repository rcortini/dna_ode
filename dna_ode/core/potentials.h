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

#ifndef __POTENTIALS_H__
#define __POTENTIALS_H__

#include "dna_ode_core.h"
#include "ode_body_functions.h"


struct potential_f {
  void (*f) (dBodyID b1, dBodyID b2, void *p, t_real *u12, t_real *f12);
  void *p;
};
typedef struct potential_f potential_f;

t_real force_loop (potential_f *U, dna_ode *dna);

void potential (const dBodyID b1, const dBodyID b2, potential_f *U, t_real *u12, t_real *f12);

/* Lennard-Jones potential */

/* the LJ parameters container */
struct LJ_parameters {
  t_real sigma;
  t_real eps;
};

void LJ_pair_potential_and_force (const dBodyID b1, const dBodyID b2, void *p, t_real *u12, t_real *f12);

void LJ_parameters_calculate (t_real l, t_real LJ_sigma, t_real LJ_eps, struct LJ_parameters *p);

/* Todd potential */

/* the Todd parameters container */
struct Todd_parameters {
  t_real CA;
  t_real r_eq;
  t_real lambda;
  t_real A;
  t_real B;
};

void Todd_pair_potential_and_force (const dBodyID b1, const dBodyID b2, void *p, t_real *u12, t_real *f12);

void Todd_parameters_calculate (t_real l, t_real CA, t_real r_eq, t_real lambda, struct Todd_parameters *p);

/* Lennard-Jones truncated potential */

/* the LJ truncated parameters container */
struct LJ_trunc_parameters {
  t_real sigma;
  t_real eps;
  t_real V0;
  t_real rc2;
};

void LJ_trunc_pair_potential_and_force (const dBodyID b1, const dBodyID b2, void *p, t_real *u12, t_real *f12);

void LJ_trunc_parameters_calculate (t_real l, t_real LJ_sigma, t_real LJ_eps, t_real r_trunc, struct LJ_trunc_parameters *p);



#endif
