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

#include "curves.h"

/* INITIALIZATION FUNCTIONS */

curve_f * curve_init (void (*r) (t_real *r, t_real s, void *p), void *p) {
  /* alloc memory */
  curve_f *f = malloc (sizeof (curve_f));
  if (f==NULL) {
    err_message ("Curve alloc failed. Out of memory\n");
    exit (EXIT_FAILURE);
  }

  /* assign curve structure parameters */
  f->r = r;
  f->p = p;

  /* set internal values */
  f->has_t = 0;
  f->has_u = 0;
  f->has_v = 0;

  /* return pointer to alloc'd struct */
  return f;
}

void curve_set_t (curve_f *f, void (*t) (t_real *r, t_real s, void *p)) {
  f->has_t = 1;
  f->t = t;
}

void curve_set_u (curve_f *f, void (*u) (t_real *r, t_real s, void *p)) {
  f->has_u = 1;
  f->u = u;
}

void curve_set_v (curve_f *f, void (*v) (t_real *r, t_real s, void *p)) {
  f->has_v = 1;
  f->v = v;
}



/* WRAPPERS. Note that these functions are not declared in headers
 * as they are to be considered only for internal use, by the numerical
 * utilities */



/* rx */
double curve_rx (double s, void *p) {
  curve_f *f = (curve_f *) p;
  t_vector r;
  CURVE_EVAL (r, s, f);
  return r[0];
}



/* ry */
double curve_ry (double s, void *p) {
  curve_f *f = (curve_f *) p;
  t_vector r;
  CURVE_EVAL (r, s, f);
  return r[1];
}



/* rz */
double curve_rz (double s, void *p) {
  curve_f *f = (curve_f *) p;
  t_vector r;
  CURVE_EVAL (r, s, f);
  return r[2];
}



/* tx */
double curve_tx (double s, void *p) {
  curve_f *f = (curve_f *) p;
  t_vector t;
  CURVE_EVAL (t, s, f);
  return t[0];
}



/* ty */
double curve_ty (double s, void *p) {
  curve_f *f = (curve_f *) p;
  t_vector t;
  CURVE_EVAL (t, s, f);
  return t[1];
}



/* tz */
double curve_tz (double s, void *p) {
  curve_f *f = (curve_f *) p;
  t_vector t;
  CURVE_EVAL (t, s, f);
  return t[2];
}



/* CALCULATE THE T, U, AND V VECTORS */



/* calculates t (s) = d r (s) / d s */
void curve_t_eval (t_real *t, t_real s, curve_f *f) {
  /* check if we have an analytical expression for t (s) */
  if (f->has_t) {
    f->t (t, s, f->p);
  }
  else {
    /* in this case, we calculate the numerical derivative of the function */
    struct f_derivative_params f_derivative_p;
    t_real tx, ty, tz;

    /* tx */
    f_derivative_p.func = curve_rx;
    f_derivative_p.p = f;
    tx = f_derivative (s, &f_derivative_p);

    /* ty */
    f_derivative_p.func = curve_ry;
    f_derivative_p.p = f;
    ty = f_derivative (s, &f_derivative_p);

    /* tz */
    f_derivative_p.func = curve_rz;
    f_derivative_p.p = f;
    tz = f_derivative (s, &f_derivative_p);

    vector_set (t, tx, ty, tz);
  }

  dNormalize3 (t); /* return normalized t */
}



/* calculates u (s) = d^2 r (s)/d s^2 */
void curve_u_eval (t_real *u, t_real s, curve_f *f) {

  /* check if we have an analytical expression for t (s) */
  if (f->has_u) {
    f->u (u, s, f->p);
  }
  else {
    /* in this case, we calculate the numerical u vector. It has
     * to be the second derivative of r(s) wrt s.
     * WARNING: if we need to calculate the numerical second derivative,
     * it becomes very unstable. Much better to provide an analytical t (s)
     * and then calculate the derivative numerically. */
    struct f_derivative_params f_derivative_p;
    t_real ux, uy, uz;

    /* ux */
    f_derivative_p.func = curve_tx;
    f_derivative_p.p = f;
    ux = f_derivative (s, &f_derivative_p);

    /* uy */
    f_derivative_p.func = curve_ty;
    f_derivative_p.p = f;
    uy = f_derivative (s, &f_derivative_p);

    /* uz */
    f_derivative_p.func = curve_tz;
    f_derivative_p.p = f;
    uz = f_derivative (s, &f_derivative_p);

    vector_set (u, ux, uy, uz);
  }

  dNormalize3 (u);
}



/* calculates v (s) = t (s) x u (s) */
void curve_v_eval (t_real *v, t_real s, curve_f *f) {
  /* check if we have an analytical expression for t (s) */
  if (f->has_v) {
    f->v (v, s, f->p);
  }
  else {
    /* In this case, we calculate v = t x u */
    t_vector t, u;
    curve_t_eval (t, s, f);
    curve_u_eval (u, s, f);
    vector_cross (v, t, u);
  }

  dNormalize3 (v);
}



/* NUMERICAL UTILITIES */



/* curve distance between s1 and s2 */
struct curve_distance_s1_s2_parameters {
  curve_f *f;
  double s1;
  double d;
};



/* This function finds the distance between r(s1) and r(s2).
 * Is used for the function s2_s1_d. */
double curve_distance_s1_s2 (double s2, void *p) {
  struct curve_distance_s1_s2_parameters *par = (struct curve_distance_s1_s2_parameters *) p;
  curve_f *f = par->f;
  double s1 = par->s1;
  t_vector r1, r2, diff;

  /* calculate r (s1) - r (s2) */
  CURVE_EVAL (r1, s1, f);
  CURVE_EVAL (r2, s2, f);
  vector_diff (diff, r1, r2);

  return vector_norm (diff) - par->d;
}



/* Calculates s2_d (s1) function. This function returns the value of the curve
 * parameter s2 corresponding to the point at which the distance between r (s1)
 * and r (s2) is equal to d */
t_real s2_s1_d (t_real s1, t_real d, curve_f *f) {
  t_real s2, upper_bound;
  int f_root_exit_status;
  struct curve_distance_s1_s2_parameters curve_distance_s1_s2_p;
  struct f_root_params f_root_p;

  /* define the parameters of the function that we need to study */
  curve_distance_s1_s2_p.f = f;
  curve_distance_s1_s2_p.s1 = s1;
  curve_distance_s1_s2_p.d = d;

  /* define the parameters of the root-finder */
  f_root_p.f = curve_distance_s1_s2;
  f_root_p.p = &curve_distance_s1_s2_p;
  f_root_p.type = gsl_root_fsolver_brent;
  f_root_p.eps_abs = 1e-7;
  f_root_p.eps_rel = 0.;
  f_root_p.max_iter = 100;

  /* find an appropriate upper bound. Note that this is only a tentative guess.
   * The distance function is not necessarily monotonic and we could find multiple solutions. */
  upper_bound = s1 + d;
  while (curve_distance_s1_s2 (upper_bound, &curve_distance_s1_s2_p)<0.) {
    upper_bound += d;
  };

  /* find root and check if we found it */
  f_root_exit_status = f_root (s1, upper_bound, &s2, &f_root_p);
  if (f_root_exit_status != GSL_SUCCESS) {
    err_message ("Error: impossible to find s2\n");
    exit (EXIT_FAILURE);
  }

  return s2;
}



/* FINALIZATION FUNCTIONS */



void curve_free (curve_f *f) {
  free (f);
}



/**********************************/
/* CURVE OPERATIONS               */
/**********************************/



/* builds the curve r_translated (s) = r (s) + R0 
 * R0 is the translation vector*/
void translated_curve (t_real *translated_vec, t_real s, void *p) {
  struct translated_curve_parameters *par = (struct translated_curve_parameters *) p;

  CURVE_EVAL (translated_vec, s, par->f);
  vector_inplace_sum (translated_vec, par->R0);
}



/* builds the curve r_rotated (s) = R r(s)
 * R is the rotation matrix which is built from rotation axis and angle */
void rotated_curve (t_real *rotated_vec, t_real s, void *p) {
  struct rotated_curve_parameters *par = (struct rotated_curve_parameters *) p;
  t_vector buff;

  CURVE_EVAL (buff, s, par->f);
  vector_rotate_around_axis (rotated_vec, par->rotation_axis, par->rotation_angle, buff);
}



/* builds the curve r_tr (s) = R r(s) + R0
 * R0 is the translation vector
 * R is the rotation matrix which is built from rotation axis and angle */
void rot_trans_curve (t_real *translated_vec, t_real s, void *p) {
  struct rot_trans_curve_parameters *par = (struct rot_trans_curve_parameters *) p;
  t_vector r0, rs, rotated_vec, rs_r0;

  dCopyVector3 (r0, par->rotation_origin);
  CURVE_EVAL (rs, s, par->f);
  vector_diff (rs_r0, rs, r0);

  vector_rotate_around_axis (rotated_vec, par->rotation_axis, par->rotation_angle, rs_r0);
  vector_sum (translated_vec, rotated_vec, par->R0);
  vector_inplace_sum (translated_vec, r0);
}



/**********************************/
/* DEFINITION OF SOME CURVES      */
/**********************************/



/* HELIX */



void helix (t_real *r, t_real s, void *p) {
  t_real x, y, z;
  struct helix_parameters *par = (struct helix_parameters *) p;

  /* calculate r */
  x = par->radius * cos (2*M_PI*s/par->pitch);
  y = par->radius * sin (2*M_PI*s/par->pitch);
  z = s;
  vector_set (r, x, y, z);
}



void helix_t (t_real *t, t_real s, void *p) {
  struct helix_parameters *par = (struct helix_parameters *) p;
  t_real tx, ty, tz;

  /* calculate t */
  tx = - 2*M_PI/par->pitch*par->radius * sin (2*M_PI*s/par->pitch);
  ty = 2*M_PI/par->pitch*par->radius * cos (2*M_PI*s/par->pitch);
  tz = 1.;
  vector_set (t, tx, ty, tz);
}



void helix_u (t_real *u, t_real s, void *p) {
  struct helix_parameters *par = (struct helix_parameters *) p;
  t_real ux, uy, uz;

  /* calculate u */
  ux = - (2*M_PI/par->pitch*par->radius)*(2*M_PI/par->pitch*par->radius) * cos (2*M_PI*s/par->pitch);
  uy = - (2*M_PI/par->pitch*par->radius)*(2*M_PI/par->pitch*par->radius) * sin (2*M_PI*s/par->pitch);
  uz = 1.;
  vector_set (u, ux, uy, uz);
}
