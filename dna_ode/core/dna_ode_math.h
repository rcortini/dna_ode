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

#ifndef __DNA_ODE_MATH_H__
#define __DNA_ODE_MATH_H__


#include "dna_ode_core.h"
#include "utils.h"


t_real signum (t_real x);

t_real safe_acos (t_real x);

t_real safe_asin (t_real x);

t_real average (t_real *v, unsigned int N);

t_real devst (t_real *v, t_real mean, unsigned int N);

void average_and_devst (t_real *v, unsigned int N, t_real *average, t_real *devst);

/* VECTOR MATH                      */

void vector_set_zero (t_real *v);

void vector_set (t_real *v, t_real x, t_real y, t_real z);

void vector_copy (t_real *res, const t_real *v);

t_real vector_norm (const t_real *v);

void vector_inplace_sum (t_real *res, const t_real *a);

void vector_sum (t_real *res, const t_real *a, const t_real *b);

void vector_inplace_diff (t_real *res, const t_real *a);

void vector_diff (t_real *res, const t_real *a, const t_real *b);

void vector_scale (t_real *res, const t_real *v, const t_real K);

void vector_lincomb (t_real *res, const t_real *u, const t_real *v, t_real lambda, t_real mu);

void vector_rotate_around_axis (t_real *res, const t_real *rotation_axis, t_real angle, const t_real *old_vector);

t_real vector_dot (const t_real *u, const t_real *v);

void vector_cross (t_real *res, const t_real *u, const t_real *v);

void axis_from_quaternion (t_real *res, const t_real *q, int i);

void quaternion_from_axes (t_real *q, const t_real *u, const t_real *v, const t_real *t);

t_real angle_between_vectors (const t_real *v1, const t_real *v2);

void vector_unit (t_real *res);

void vector_project_2d_unit (t_real *res, const t_real *v);

/* MATRIX MATH                      */

void matrix_set (t_real *res, const t_real *c1, const t_real *c2, const t_real *c3);

void matrix_copy (t_real *m, const t_real *m1);

void matrix_set_identity (t_real *m);

void matrix_negate (t_real *m);

void matrix_transpose (t_real *m, const t_real *m1);

void matrix_inplace_sum (t_real *m, const t_real *m1);

void matrix_sum (t_real *res, const t_real *m1, const t_real *m2);

void matrix_inplace_diff (t_real *m, const t_real *m1);

void matrix_diff (t_real *res, const t_real *m1, const t_real *m2);

void matrix_scale (t_real *res, const t_real *m, const t_real r);

void matrix_vector_cross (t_real *res, const t_real *u);

t_real matrix_trace (const t_real *m);

void matrix_vector_product (t_real *res, const t_real *m, const t_real *v);

void matrix_product (t_real *res, const t_real *m1, const t_real *m2);

void matrix_ABAt (t_real *m, const t_real *A, const t_real *B);

t_real matrix_determinant (const t_real *m);

void matrix_inverse_CH (t_real *res, const t_real *m);

void matrix_inverse (t_real *res, const t_real *m);

t_real langevin_function (t_real x);

/* GSL MATH WRAPPER FUNCTIONS */

struct f_root_params {
  double (*f) (double x, void *p);
  void *p;
  const gsl_root_fsolver_type *type;
  double eps_abs;
  double eps_rel;
  unsigned int max_iter;
};

int f_root (double x_min, double x_max, double *root, struct f_root_params *par);

struct f_multimin_par {
  unsigned int dim;
  t_real size_tol;
  unsigned int max_iter;
  gsl_vector *step_size;
  const gsl_multimin_fminimizer_type *type;
  t_real (*func) (const gsl_vector *X, void *p);
  void *func_par;
};

int f_multimin (gsl_vector *x_start, gsl_vector *x_min, struct f_multimin_par *par);

struct f_derivative_params {
  double (*func) (double x, void *p);
  void *p;
};

double f_derivative (double x, struct f_derivative_params *p);

/* this structure contains all the information necessary to perform the
 * non-linear least-square fit */
struct fdf_fit_parameters {
  size_t n;
  double *x;
  double *y;
  double *sigma;
  double (*model_f) (double x, const gsl_vector *par);
  double (*model_df) (unsigned int i, double x, const gsl_vector *par);
  size_t p;
  struct data *d;
  double eps_abs;
  double eps_rel;
  size_t max_iter;
  const gsl_multifit_fdfsolver_type *type;
};

int chi_f (const gsl_vector *X, void *p, gsl_vector *f);

int chi_df (const gsl_vector *X, void *par, gsl_matrix * J);

int chi_fdf (const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix * J);

double chi2_from_fit (gsl_vector *fit, struct fdf_fit_parameters *fit_pars);

/* DISTANCE MAP FUNCTIONS */

/* distance map structure definition */
struct distance_map {
  unsigned int nbodies;
  unsigned int dim;
  t_real *distance;
};
typedef struct distance_map distance_map;

distance_map * distance_map_alloc (unsigned int nbodies);

distance_map * distance_map_calloc (unsigned int nbodies);

int get_map_index_ij (unsigned int i, unsigned int j);

void set_distance_ij (unsigned int i, unsigned int j, t_real distance, distance_map *map);

t_real get_distance_ij (unsigned int i, unsigned int j, distance_map *map);

void print_distance_map (distance_map *map);

void file_print_distance_map (distance_map *map, FILE *file);

void distance_map_sum (distance_map *dest, distance_map *map);

void distance_map_scale (distance_map *map, t_real k);

void distance_map_free (distance_map *map);

void set_gsl_matrix_from_matrix (unsigned int i, unsigned int j, const t_real *v, gsl_matrix *gsl_m);

void set_gsl_vector_from_vector (unsigned int i, const t_real *v, gsl_vector *gsl_v);

void check_LDL_decomposition (unsigned int N,
    t_matrix *Lij,
    t_matrix *Dii,
    t_matrix *Hij,
    t_matrix *Hii);

void check_forward_substitution (unsigned int N,
    t_matrix *Lij,
    t_vector *rhs,
    t_vector *buffer);

void check_back_substitution (unsigned int N,
    t_matrix *Lij,
    t_matrix *Dii,
    t_vector *lambda,
    t_vector *buffer);

void check_LDL_solution (unsigned int N,
  t_matrix *Hii,
  t_matrix *Hij,
  t_vector *lambda,
  t_vector *rhs);

#endif
