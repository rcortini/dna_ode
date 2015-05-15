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

#ifndef __CURVES_H__
#define __CURVES_H__

#include "dna_ode_core.h"
#include "dna_ode_math.h"

/* CURVE FUNCTION PROTOTYPE DEFINTIONS */

struct curve_f {
  void (*r) (t_real *r, t_real s, void *p);
  void *p;

  void (*t) (t_real *r, t_real s, void *p);
  void (*u) (t_real *r, t_real s, void *p);
  void (*v) (t_real *r, t_real s, void *p);

  unsigned int has_t;
  unsigned int has_u;
  unsigned int has_v;
};

typedef struct curve_f curve_f;

#define CURVE_EVAL(vec,s,f) ((f)->r((vec),(s),((f)->p)))

curve_f * curve_init (void (*r) (t_real *r, t_real s, void *p), void *p);

void curve_set_t (curve_f *f, void (*t) (t_real *r, t_real s, void *p));

void curve_set_u (curve_f *f, void (*u) (t_real *r, t_real s, void *p));

void curve_set_v (curve_f *f, void (*v) (t_real *r, t_real s, void *p));

void curve_t_eval (t_real *, t_real s, curve_f *f);

void curve_u_eval (t_real *, t_real s, curve_f *f);

void curve_v_eval (t_real *, t_real s, curve_f *f);

t_real s2_s1_d (t_real s1, t_real d, curve_f *f);

void curve_free (curve_f *f);



/**********************************/
/* CURVE OPERATIONS               */
/**********************************/



struct translated_curve_parameters {
  curve_f *f;
  t_real *R0;
};

void translated_curve (t_real *, t_real s, void *p);

struct rotated_curve_parameters {
  curve_f *f;
  t_vector rotation_axis;
  t_real rotation_angle;
};

void rotated_curve (t_real *, t_real s, void *p);

struct rot_trans_curve_parameters {
  curve_f *f;
  t_vector R0;
  t_vector rotation_axis;
  t_vector rotation_origin;
  t_real rotation_angle;
};

void rot_trans_curve (t_real *, t_real s, void *p);



/**********************************/
/* DEFINITION OF SOME CURVES      */
/**********************************/



/* HELIX */



struct helix_parameters {
  t_real radius;
  t_real pitch;
};

void helix (t_real *, t_real s, void *p);

void helix_t (t_real *, t_real s, void *p);

void helix_u (t_real *, t_real s, void *p);

#endif
