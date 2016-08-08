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

#ifndef __MECHANICAL_MODELS_H__
#define __MECHANICAL_MODELS_H__

#include "dna_ode_core.h"
#include "dna_ode_math.h"

/* these macros are useful not to rewrite the same things throughout
 * this file */

#define T1_T2(R1,R2) dCalcVectorDot3_44 (R1+2, R2+2)
#define T1_x_T2(bn,R1,R2) dCalcVectorCross3_144 (bn, R1+2, R2+2)
#define U1_U2(R1,R2) dCalcVectorDot3_44 (R1, R2)
#define V1_V2(R1,R2) dCalcVectorDot3_44 (R1+1, R2+1)
#define V1_U2(R1,R2) dCalcVectorDot3_44 (R1+1, R2)
#define U1_V2(R1,R2) dCalcVectorDot3_44 (R1, R2+1)

/* calculates the twist between two bodies, once
 * 1 + t1 * t2 is known */
static inline void twist_12 (const t_real *R1, const t_real *R2,
    t_real one_plus_t1_dot_t2, t_real *phi) {
  t_real cTw, sTw;

  /* calculate cosine and sine of the twist */
  cTw = (U1_U2 (R1,R2) + V1_V2 (R1,R2))/one_plus_t1_dot_t2;
  sTw = (V1_U2 (R1,R2) - U1_V2 (R1,R2))/one_plus_t1_dot_t2;

  /* calculate the twist, with correct sign */
  *phi = safe_acos (cTw) * signum (sTw);
}

static inline t_real theta_angle (const t_real *R1, const t_real *R2) {
  return fabs (safe_acos (T1_T2 (R1,R2)));
}

static inline void t1_plus_t2_vec (t_real *res, const t_real *R1, const t_real *R2) {
  vector_set (res,
      R1 [2] + R2 [2],
      R1 [6] + R2 [6],
      R1 [10] + R2 [10]);
}

/* calculates the twist angle (times 2pi) for a given pair of bodies */
static inline t_real body_twist (const t_real *R1, const t_real *R2) {
  t_real one_plus_t1_dot_t2, Tw;

  /* calculate the twist */
  one_plus_t1_dot_t2 = 1. + T1_T2 (R1, R2);
  twist_12 (R1, R2, one_plus_t1_dot_t2, &Tw);

  /* return the calculated value of the twist */
  return Tw;
}

struct mechanical_model_f {
  void (*f) (t_real *, const t_real *, const t_real *, t_real *, void *);
  void *p;
};
typedef struct mechanical_model_f mechanical_model_f;

void body_twist_and_torque_kinkable (t_real *torque, const t_real *, const t_real *, t_real *Tw, void *);

void body_twist_and_torque_harmonic (t_real *torque, const t_real *, const t_real *, t_real *Tw, void *);

void body_twist_and_torque_nicked_kinkable (t_real *torque, const t_real *, const t_real *, t_real *Tw, void *);

void body_twist_and_torque_nicked_harmonic (t_real *torque, const t_real *, const t_real *, t_real *Tw, void *);

void body_twist_and_torque_harmonic4 (t_real *torque, const t_real *, const t_real *, t_real *Tw, void *);

void body_twist_and_torque_nicked_harmonic4 (t_real *torque, const t_real *, const t_real *, t_real *Tw, void *);

typedef struct {
  unsigned int njoints;
  t_real *gb;
  t_real *gt;
  unsigned int count;
} inhom_par;

void body_twist_and_torque_harmonic_inhom (t_real *torque, const t_real *R1, const t_real *R2, t_real *Tw, void *p);

#endif
