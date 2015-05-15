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

#include "dna_ode_math.h"

/*
 * MATH FUNCTIONS
 */

/* the signum function */
t_real signum (t_real x) {
  if (x>0)
    return 1.;
  else if (x<0)
    return -1.;
  else
    return 0.;
}



/* this is a "safe" call to the arccosine function. It takes care of the case in which the argument of acos is slightly larger than +/- 1 */
t_real safe_acos (t_real x) {
  if (x < -1.0) x = -1.0;
  else if (x > 1.0) x = 1.0;
  return acos (x);
}



/* this is a "safe" call to the arcsine function. It takes care of the case in which the argument of acos is slightly larger than +/- 1 */
t_real safe_asin (t_real x) {
  if (x < -1.0) x = -1.0;
  else if (x > 1.0) x = 1.0;
  return asin (x);
}



/* computes average of an array of real numbers */
t_real average (t_real *v, unsigned int N) {
  unsigned int i;
  double sum, av;

  /* check that we have at least one value */
  if (N <= 1) {
    err_message ("Insufficient data\n");
    exit (EXIT_FAILURE);
  }

  /* cycle until the end of the input stream is reached */
  sum = 0.;

  /* read input */
  for (i=0; i<N; i++) sum += v [i];

  /* calculate average and return */
  av = sum/N;
  return av;
}



/* computes standard deviation of an array of real numbers, given the mean in input */
t_real devst (t_real *v, t_real mean, unsigned int N) {
  unsigned int i;
  double sum, ds;

  /* check that we have at least one value */
  if (N <= 2) {
    err_message ("Insufficient data\n");
    exit (EXIT_FAILURE);
  }

  /* cycle until the end of the input stream is reached */
  sum = 0.;

  /* read input */
  for (i=0; i<N; i++) sum += (v [i]-mean)*(v [i]-mean);

  /* calculate average and return */
  ds = sqrt (sum/(N-1.));
  return ds;
}



/* computes average and standard deviation of an array of real numbers. It does
 * it in one go, one for cycle (more efficient) */
void average_and_devst (t_real *v, unsigned int N, t_real *av, t_real *ds) {
  unsigned int i;
  double value, sum, sum2;

  /* check that we have at least one value */
  if (N <= 2) {
    err_message ("Insufficient data\n");
    exit (EXIT_FAILURE);
  }

  /* read input */
  sum=sum2=0.;
  for (i=0; i<N; i++) {
    /* calculate the sum and the sum of the squares of the input data */
    value = v [i];
    sum+=value;
    sum2+=value*value;
  }

  /* calculate average and standard deviation */
  *av = sum/N;
  
  /* calculate standard deviation */
  *ds = sqrt (N/(N-1.)*(sum2/N - (*av)* (*av)));
  *ds /= sqrt (N);
}



/************************************
 * VECTOR MATH                      *
 ***********************************/



/* set to null vector */
void vector_set_zero (t_real *v) {
  dSetZero (v,3);
}



/* assigns the components of a vector */
void vector_set (t_real *v, t_real x, t_real y, t_real z) {
  v [0] = x;
  v [1] = y;
  v [2] = z;
}



/* creates an identical copy of a vector */
void vector_copy (t_real *res, const t_real *v) {
  res [0] = v [0];
  res [1] = v [1];
  res [2] = v [2];
}



/* returns the norm of a vector */
t_real vector_norm (const t_real *v) {
  return dCalcVectorLength3 (v);
}



/* calculates res = res + a */
void vector_inplace_sum (t_real *res, const t_real *a) {
  res [0] += a [0];
  res [1] += a [1];
  res [2] += a [2];
}



/* assigns c = a + b, where a, b, and c are vectors */
void vector_sum (t_real *res, const t_real *a, const t_real *b) {
  dAddVectors3 (res, a, b);
}



/* calculates res = res - a */
void vector_inplace_diff (t_real *res, const t_real *a) {
  res [0] -= a [0];
  res [1] -= a [1];
  res [2] -= a [2];
}



/* assigns c = a - b, where a, b, and c are vectors */
void vector_diff (t_real *res, const t_real *a, const t_real *b) {
  dSubtractVectors3 (res, a, b);
}



/* assigns c = K v, with c and v vectors and K scalar */
void vector_scale (t_real *res, const t_real *v, t_real K) {
  dCopyScaledVector3 (res, v, K);
}



/* this assigns c = lambda*v + mu*u */
void vector_lincomb (t_real *res, const t_real *u, const t_real *v, t_real lambda, t_real mu) {
  dAddScaledVectors3 (res, u, v, lambda, mu);
}



/* calculates the dot product between vectors u and v */
t_real vector_dot (const t_real *u, const t_real *v) {
  return dDot (u, v, 3);
}



/* returns a vector that is the cross product of v and u */
void vector_cross (t_real *res, const t_real *u, const t_real *v) {
  dCalcVectorCross3 (res, u, v);
}



/* rotates old_vector around rotation_axis by angle, and assigns the result to new_vector */
void vector_rotate_around_axis (t_real *res, const t_real *rotation_axis, t_real angle, const t_real *old_vector) {
  dMatrix3 mr;

  /* check if the length of the rotation_axis vector is not null */
  if (dCalcVectorLength3 (rotation_axis)!=0.)
    dRFromAxisAndAngle (mr, rotation_axis[0], rotation_axis[1], rotation_axis[2], angle);
  else
    dRSetIdentity(mr);

  /* assign the result */
  dMultiply0 (res, mr, old_vector, 3,3,1);
}



/* returns the axis i, given a quaternion.
 * i = 0 -> u
 * i = 1 -> v
 * i = 2 -> t
 * CODE DIRECTLY TRANSLATED FROM ODE's rotation.cpp */
void axis_from_quaternion (t_real *axis, const t_real *q, int i) { 
  switch (i) {
    case 0:
      axis[0] = 1 - 2 * (q[2]*q[2] + q[3]*q[3]);
      axis[1] = 2 * (q[1]*q[2] + q[0]*q[3]);
      axis[2] = 2 * (q[1]*q[3] - q[0]*q[2]);
      break;
    case 1:
      axis[0] = 2 * (q[1]*q[2] - q[0]*q[3]);
      axis[1] = 1 - 2 * (q[1]*q[1] + q[3]*q[3]);
      axis[2] = 2 * (q[2]*q[3] + q[0]*q[1]);
      break;
    case 2:
      axis[0] = 2 * (q[1]*q[3] + q[0]*q[2]);
      axis[1] = 2 * (q[2]*q[3] - q[0]*q[1]);
      axis[2] = 1 - 2 * (q[1]*q[1] + q[2]*q[2]);
      break;
    default:
      err_message ("Axis index out of range\n");
      exit (EXIT_FAILURE);
  }
}



/* calculates the quaternion, given a set of three body vectors u, v, t.
 * CODE DIRECTLY TRANSLATED FROM ODE's rotation.cpp */
void quaternion_from_axes (t_real *q, const t_real *u, const t_real *v, const t_real *t) {
  t_real tr, s;
  tr = u[0] + v[1] + t[2];
  if (tr >= 0) {
    s = sqrt (tr + 1.);
    q [0] = 0.5 * s;
    s = 0.5/s;
    q [1] = (v[2] - t[1]) * s;
    q [2] = (t[0] - u[2]) * s;
    q [3] = (u[1] - v[0]) * s;
  }
  else {
    /* find the largest diagonal element and jump to the appropriate case */
    t_real max = GSL_MAX (u[0], GSL_MAX (v[1], t[2]));
    if (max == u[0]) {
      s = sqrt (1. + u[0] - v[1] - t[2]);
      q [1] = 0.5 * s;
      s = 0.5/s;
      q [2] = (u[1] + v[0]) * s;
      q [3] = (t[0] + u[2]) * s;
      q [0] = (v[2] - t[1]) * s;
    }
    else if (max == v[1]) {
      s = sqrt (1. + v[1] - t[2] - u[0]);
      q [2] = 0.5 * s;
      s = 0.5/s;
      q [3] = (v[2] + t[1]) * s;
      q [1] = (u[1] + v[0]) * s;
      q [0] = (t[0] - u[2]) * s;
    }
    else {
      s = sqrt (1. + t[2] - u[0] - v[1]);
      q [3] = 0.5 * s;
      s = 0.5/s;
      q [1] = (t[0] + u[2]) * s;
      q [2] = (v[2] + t[1]) * s;
      q [0] = (u[1] - v[0]) * s;
    }
  }
}



/* returns the angle between two vectors, in radians */
t_real angle_between_vectors (const t_real *v1, const t_real *v2) {
  return safe_acos (vector_dot (v1, v2)/(vector_norm (v1) * vector_norm (v2)));
}



/* returns the unit vector corresponding to the given vector */
void vector_unit (t_real *res) {
  t_real norm = vector_norm (res);
  if (norm==0.)
    vector_set (res, 0., 0., 0.);
  else
    dNormalize3 (res);
}



/* returns the unit vector corresponding to the given vector */
void vector_project_2d_unit (t_real *res, const t_real *v) {
  t_real norm = sqrt (v[0]*v[0] + v[1]*v[1]);
  res [0] = v[0]/norm;
  res [1] = v[1]/norm;
  res [2] = 0.;
}


/************************************
 * MATRIX MATH                      *
 ***********************************/



/* creates a matrix out of three vectors */
void matrix_set (t_real *m, const t_real *c1, const t_real *c2, const t_real *c3) {
  /* row 1 */
  m [0] = c1 [0];
  m [1] = c1 [1];
  m [2] = c1 [2];

  /* row 2 */
  m [4] = c2 [0];
  m [5] = c2 [1];
  m [6] = c2 [2];

  /* row 3 */
  m [8] = c3 [0];
  m [9] = c3 [1];
  m [10] = c3 [2];

  /* extra elements */
  m [3] = m [7] = m [11] = 0.;
}



/* copies a matrix */
void matrix_copy (t_real *m, const t_real *m1) {
  /* row 1 */
  m [0] = m1 [0];
  m [1] = m1 [1];
  m [2] = m1 [2];

  /* row 2 */
  m [4] = m1 [4];
  m [5] = m1 [5];
  m [6] = m1 [6];

  /* row 3 */
  m [8] = m1 [8];
  m [9] = m1 [9];
  m [10] = m1 [10];

  /* extra elements */
  m [3] = m [7] = m [11] = 0.;
}



/* set identity matrix */
void matrix_set_identity (t_real *m) {
  t_vector c1, c2, c3;
  vector_set (c1, 1., 0., 0.);
  vector_set (c2, 0., 1., 0.);
  vector_set (c3, 0., 0., 1.);
  matrix_set (m, c1, c2, c3);
}



/* the transpose of a matrix */
void matrix_transpose (t_real *m, const t_real *m1) {
  /* row 1 */
  m [0] = m1 [0];
  m [1] = m1 [4];
  m [2] = m1 [8];

  /* row 2 */
  m [4] = m1 [1];
  m [5] = m1 [5];
  m [6] = m1 [9];

  /* row 3 */
  m [8] = m1 [2];
  m [9] = m1 [6];
  m [10] = m1 [10];

  /* extra elements */
  m [3] = m [7] = m [11] = 0.;
}


void matrix_negate (t_real *m) {
  dNegateVector3 (m  );
  dNegateVector3 (m+4);
  dNegateVector3 (m+8);
}


/* m = m + m1 */
void matrix_inplace_sum (t_real *m, const t_real *m1) {
  vector_inplace_sum (m    , m1);
  vector_inplace_sum (m + 4, m1 + 4);
  vector_inplace_sum (m + 8, m1 + 8);
}



/* m = m - m1 */
void matrix_inplace_diff (t_real *m, const t_real *m1) {
  vector_inplace_diff (m    , m1);
  vector_inplace_diff (m + 4, m1 + 4);
  vector_inplace_diff (m + 8, m1 + 8);
}



/* sum of two matrices */
void matrix_sum (t_real *m, const t_real *m1, const t_real *m2) {
  vector_sum (m    , m1    , m2);
  vector_sum (m + 4, m1 + 4, m2 + 4);
  vector_sum (m + 8, m1 + 8, m2 + 8);
  m [3] = m [7] = m [11] = 0.;
}



/* difference of two matrices */
void matrix_diff (t_real *m, const t_real *m1, const t_real *m2) {
  vector_diff (m    , m1    , m2);
  vector_diff (m + 4, m1 + 4, m2 + 4);
  vector_diff (m + 8, m1 + 8, m2 + 8);
  m [3] = m [7] = m [11] = 0.;
}


/* matrix scale by a real number */
void matrix_scale (t_real *n, const t_real *m, t_real r) {
  vector_scale (n,     m    , r);
  vector_scale (n + 4, m + 4, r);
  vector_scale (n + 8, m + 8, r);
}



/* creates a matrix associated with cross product */
void matrix_vector_cross (t_real *m, const t_real *u) {
  dSetZero (m, 12);
  dSetCrossMatrixMinus (m, u, 4);
}



/* returns the trace of a matrix */
t_real matrix_trace (const t_real *m) {
  return m[0] + m[5] + m[10];
}



/* matrix-vector product */
void matrix_vector_product (t_real *res, const t_real *m, const t_real *v) {
  dMultiply0 (res, m, v, 3, 3, 1);
}



/* matrix-matrix product */
void matrix_product (t_real *m, const t_real *m1, const t_real *m2) {
  dMultiply0 (m, m1, m2, 3, 3, 3);
  m[3] = m[7] = m[11] = 0.;
}


/* calculates the matrix m = A B A^t */
void matrix_ABAt (t_real *m, const t_real *A, const t_real *B) {
  t_matrix BAt, At;
  matrix_transpose (At, A);
  matrix_product (BAt, B, At);
  matrix_product (m, A, BAt);
}



/* matrix determinant */
t_real matrix_determinant (const t_real *m) {
  return dCalcMatrix3Det (m);
}



/* matrix inverse using Caley - Hamilton method */
/* void matrix_inverse_CH (t_real *res, const t_real *m) {
  t_real tr, dtr;
  t_matrix m1, m2;

  tr = matrix_trace (m);
  matrix_product (m2, m, m);
  dtr = .5*(tr*tr - matrix_trace (m2));

  m1 = matrix_sum (matrix_scale (m, -tr), m2);
  m1.c1.x += dtr;
  m1.c2.y += dtr;
  m1.c3.z += dtr;
  m1 = matrix_scale (m1, 1./matrix_determinant (m));

  return m1;
} */



/* invert matrix and return its determinant */
void matrix_inverse (t_real *res, const t_real *m) {
  dInvertMatrix3 (res, m);
  res [3] = res [7] = res [11] = 0.;
}



/*****************************************
 * GENERAL UTILITY FUNCTIONS             *
 * ***************************************/



/* Langevin function */
t_real langevin_function (t_real x) {
  return 1./tanh (x) - 1./x;
}



/****************************************
 * METHODS INHERITED FROM GSL FUNCTIONS *
 * **************************************/



/* this function returns the solution to f (x) = 0 */
int f_root (double x_min, double x_max, double *root, struct f_root_params *par) {
  unsigned int iter;
  int status;
  double x_lo, x_hi;
  gsl_root_fsolver *s;

  /* defines the function to pass to the solver */
  gsl_function func;
  func.function = par->f;
  func.params = par->p;

  /* allocates and initializes the minimizer */
  s = gsl_root_fsolver_alloc (par->type);
  gsl_root_fsolver_set (s, &func, x_min, x_max);

  /* start the iteration to find the root */
  iter = 0;
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    *root = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, par->eps_abs, par->eps_rel);

    DPRINT ("%d: x = %f\n", iter, *root);
  }
  while (status==GSL_CONTINUE && iter<par->max_iter);

  /* free the memory and return the found root */
  gsl_root_fsolver_free (s);
  return status;
}



/* calculates the numerical derivative of a function at x */
double f_derivative (double x, struct f_derivative_params *p) {
  double dx, df, err;
  gsl_function func;
  double h = GSL_SQRT_DBL_EPSILON;

  /* set the fields for the function */
  func.function = p->func;
  func.params = p->p;

  /* set the increment of the function */
  dx = x==0. ? h : x*h;

  /* invokes the derivative calculator and returns */
  gsl_deriv_central (&func, x, dx, &df, &err);
  /* TODO: check err */
  return df;
}



/* This function finds the minimum of a real-valued vector-argument function,
 * given x_start as starting point. Uses algorithm of minimization without
 * derivatives. The parameters are passed through the f_multimin_par structure,
 * which also contains the pointer to the function to minimize, and the parameters
 * of the function to minimize. */
int f_multimin (gsl_vector *x_start, gsl_vector *x_min, struct f_multimin_par *par) {
  unsigned int iter;
  int status;
  t_real size;
  gsl_multimin_fminimizer *s;
  gsl_multimin_function my_func;

  /* initializes the function to minimize */
  my_func.n = par->dim;
  my_func.f = par->func;
  my_func.params = par->func_par;

  /* initializes the minimizer */
  s = gsl_multimin_fminimizer_alloc (par->type, par->dim);
  gsl_multimin_fminimizer_set (s, &my_func, x_start, par->step_size);

  /* iterate the minimizer */
  iter = 0;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);

    if (status) {
      err_message ("f_multimin: warning: %s\n", gsl_strerror (status));
      break;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, par->size_tol);
  }
  while (status==GSL_CONTINUE && iter<par->max_iter);

  /* set the minimum, free memory, and exit */
  gsl_vector_memcpy (x_min, s->x);
  gsl_multimin_fminimizer_free (s);
  return GSL_SUCCESS;
}



/* this function calculates chi = F_i (x_i, y_i, p) = f(x_i, p) - y_i/sigma_i */
int chi_f (const gsl_vector *X, void *p, gsl_vector *f) {
  struct fdf_fit_parameters *fit_p = (struct fdf_fit_parameters *) p;
  size_t n = fit_p->n;
  double *x = fit_p->x;
  double *y = fit_p->y;
  double *sigma = fit_p->sigma;
  size_t i;

  for (i=0; i<n; i++) {
    double Yi = fit_p->model_f (x [i], X);
    gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
  }

  return GSL_SUCCESS;
}



/* calculates the derivative components of the chi function */
int chi_df (const gsl_vector *X, void *par, gsl_matrix * J) {
  struct fdf_fit_parameters *fit_p = (struct fdf_fit_parameters *) par;
  size_t n = fit_p->n, p = fit_p->p;
  double *sigma = fit_p->sigma;
  double *x = fit_p->x;
  size_t i, j;

  for (i=0; i<n; i++) {
    /* Jacobian matrix J(i,j) = dfi / dxj, */
    /* where fi = (Yi - yi)/sigma[i]       */
    double s = sigma [i];
    for (j=0; j<p; j++) {
      gsl_matrix_set (J, i, j, fit_p->model_df (j, x[i], X)/s); 
    }
  }
  return GSL_SUCCESS;
}



/* calculates chi function and its derivatives at the same time */
int chi_fdf (const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix * J) {
  chi_f (x, p, f);
  chi_df (x, p, J);

  return GSL_SUCCESS;
}



/* performs a non-linear best-fit of parameters from a set of weighted data and a model, and
 * its derivatives, based on the algorithm passed through the fit_p pointer, using a least-square minimization
 * method */
int fdf_fit (const gsl_vector *x_start, struct fdf_fit_parameters *fit_p, gsl_vector *fit, gsl_matrix *covar) {
  unsigned int iter = 0, p = fit_p->p, n = fit_p->n;
  const gsl_multifit_fdfsolver_type *T = fit_p->type;
  gsl_multifit_fdfsolver *s;
  int status;
  gsl_multifit_function_fdf f;

  /* init function to fit */
  f.f = &chi_f;
  f.df = &chi_df;
  f.fdf = &chi_fdf;
  f.n = n;
  f.p = fit_p->p;
  f.params = fit_p;

  /* initialize the fitter */
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, x_start);

  /* iterate */
  do
  {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s);

    /* printf ("status = %s\n", gsl_strerror (status)); */

    if (status)
      break;

    status = gsl_multifit_test_delta (s->dx, s->x, fit_p->eps_abs, fit_p->eps_rel);
  }
  while (status == GSL_CONTINUE && iter < fit_p->max_iter);

  /* assign the fit vector and the covariance matrix */
  gsl_multifit_covar (s->J, 0.0, covar);
  gsl_vector_memcpy (fit, s->x);

  /* printf ("status = %s\n", gsl_strerror (status)); */

  gsl_multifit_fdfsolver_free (s);

  return status;
}



/* calculates the value of the chi^2, as a function of the best-fit vector of 
 * parameters, and the parameters of the function fitter */
double chi2_from_fit (gsl_vector *fit, struct fdf_fit_parameters *fit_pars) {
  unsigned int i;
  double chi2 = 0.;
  gsl_vector *f = gsl_vector_alloc (fit_pars->n);
  chi_f (fit, fit_pars, f);

  for (i=0; i<fit_pars->n; i++) {
    double fi = gsl_vector_get (f, i);
    chi2 += fi*fi;
  }

  gsl_vector_free (f);
  return chi2;
}



/* DISTANCE MAP FUNCTIONS */



/* alloc a distance map */
distance_map * distance_map_alloc (unsigned int nbodies) {
  unsigned int dim = nbodies * (nbodies-1)/2;
  t_real * distance = (t_real *) malloc (dim * sizeof (t_real));
  distance_map * map = (distance_map *) malloc (sizeof (*map));
  if (map==NULL) {
    err_message ("Insufficient memory\n");
    exit (EXIT_FAILURE);
  }
  map->nbodies = nbodies;
  map->dim = dim;
  map->distance = distance;
  return map;
}



/* calloc a distance map */
distance_map * distance_map_calloc (unsigned int nbodies) {
  unsigned int i, dim = nbodies * (nbodies-1)/2;
  t_real * distance = (t_real *) malloc (dim * sizeof (t_real));
  distance_map * map = (distance_map *) malloc (sizeof (*map));
  map->nbodies = nbodies;
  map->dim = dim;
  map->distance = distance;
  for (i=0; i<dim; i++) map->distance [i] = 0.0;
  return map;
}



/* get the index of the distance vector that corresponds to the distance between i and j */
int get_map_index_ij (unsigned int i, unsigned int j) {
  unsigned int max_index = i>j ? i : j;
  unsigned int min_index = i<j ? i : j;
  return (max_index)*(max_index-1)/2 + min_index;
}



/* sets the distance between body i and j in the distance map */
void set_distance_ij (unsigned int i, unsigned int j, t_real distance, distance_map *map) {
  if (i==j)
    return;
  else if (i>=map->nbodies || j>=map->nbodies) {
    err_message ("Error. Index exceeding nbodies\n");
    distance_map_free (map);
    exit (EXIT_FAILURE);
  }
  else {
    int map_index = get_map_index_ij (i, j);
    map->distance [map_index] = distance;
  }
}



/* gets the distance between body i and j in the distance map */
t_real get_distance_ij (unsigned int i, unsigned int j, distance_map *map) {
  if (i==j)
    return 0.0;
  else if (i>=map->nbodies || j>=map->nbodies) {
    err_message ("Error. Index exceeding nbodies\n");
    distance_map_free (map);
    exit (EXIT_FAILURE);
  }
  else {
    unsigned int map_index = get_map_index_ij (i, j);
    return map->distance [map_index];
  }
}



/* print the distance map */
void print_distance_map (distance_map *map) {
  unsigned int i, j;
  const unsigned int nbodies = map->nbodies;

  /* print the map */
  for (i=1; i<nbodies; i++)
    for (j=0; j<i; j++)
      message ("%d %d %f\n", i, j, get_distance_ij (i, j, map));
}



/* print the distance map */
void file_print_distance_map (distance_map *map, FILE *file) {
  unsigned int i, j;
  const unsigned int nbodies = map->nbodies;

  /* print the map */
  for (i=1; i<nbodies; i++)
    for (j=0; j<i; j++)
      file_message (file, "%d %d %f\n", i, j, get_distance_ij (i, j, map));
}



/* sum two distance maps element by element */
void distance_map_sum (distance_map *dest, distance_map *map) {
  const unsigned int dim = map->dim;
  unsigned int i;

  /* now set the values of the contact map */
  for (i=0; i<dim; i++) {
    dest->distance [i] = dest->distance [i] + map->distance [i];
  }
}



void distance_map_scale (distance_map *map, t_real k) {
  const unsigned int dim = map->dim;
  unsigned int i;

  /* now scale the values */
  for (i=0; i<dim; i++) {
    map->distance [i] /= k;
  }
}



/* free a distance map */
void distance_map_free (distance_map *map) {
  free (map->distance);
  free (map);
}



/* set the elements of the gsl_vector gsl_v according to the vector v. This sets a block in the
 * gsl_vector */
void set_gsl_vector_from_vector (unsigned int i, const t_real *v, gsl_vector *gsl_v) {
  unsigned int n;
  const unsigned int ind = 3*i;
  for (n=0; n<3; n++) 
    gsl_vector_set (gsl_v, ind+n, v [n]);
}



/* set the elements of the gsl_matrix gsl_m according to the matrix M. This sets a block in the
 * gsl_matrix */
void set_gsl_matrix_from_matrix (unsigned int i, unsigned int j, const t_real *M, gsl_matrix *gsl_m) {
  unsigned int n, m;
  const unsigned int row = 3*i;
  const unsigned int col = 3*j;
  for (n=0; n<3; n++) 
    for (m=0; m<3; m++) 
      gsl_matrix_set (gsl_m, row+n, col+m, M [4*n+m]);
}



void check_LDL_decomposition (unsigned int N,
    t_matrix *Lij,
    t_matrix *Dii,
    t_matrix *Hij,
    t_matrix *Hii) {
  t_matrix id;
  unsigned int i;
  /* allocate the memory for the gsl matrices */
  gsl_matrix *H = gsl_matrix_alloc (3*N, 3*N);
  gsl_matrix *LDLt = gsl_matrix_alloc (3*N, 3*N);
  gsl_matrix *DLt = gsl_matrix_alloc (3*N, 3*N);
  gsl_matrix *D = gsl_matrix_alloc (3*N, 3*N);
  gsl_matrix *L = gsl_matrix_alloc (3*N, 3*N);
  gsl_matrix *diff = gsl_matrix_alloc (3*N, 3*N);

  matrix_set_identity (id);

  /* build the gsl matrices */
  for (i=0; i<N; i++) {
    /* H */
    set_gsl_matrix_from_matrix (i,i, Hii [i], H);
    if (i>0) {
      t_matrix transp;
      set_gsl_matrix_from_matrix (i,i-1, Hij [i-1], H);
      matrix_transpose (transp, Hij [i-1]);
      set_gsl_matrix_from_matrix (i-1,i, transp, H);
    }

    /* D */
    set_gsl_matrix_from_matrix (i,i, Dii [i], D);

    /* L */
    set_gsl_matrix_from_matrix (i,i, id, L);
    if (i>0)
      set_gsl_matrix_from_matrix (i,i-1, Lij [i-1], L);
  }

  /* calculate D L^t */
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1., D, L, 0., DLt);

  /* calculate L D L^t */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1., L, DLt, 0., LDLt);

  /* check LDL^t = H */
  gsl_matrix_memcpy (diff, H);
  gsl_matrix_sub (diff, LDLt);

  /* print results */
  printf ("LDLt:\n");
  print_gsl_matrix (3*N, 3*N, LDLt);
  printf ("H:\n");
  print_gsl_matrix (3*N, 3*N, H);
  printf ("diff:\n");
  print_gsl_matrix (3*N, 3*N, diff);

  /* free memory */
  gsl_matrix_free (H);
  gsl_matrix_free (LDLt);
  gsl_matrix_free (DLt);
  gsl_matrix_free (D);
  gsl_matrix_free (L);
}



/* check the equation L Y = X, as calculated using
 * L = Lij, Y = buffer, X = rhs */
void check_forward_substitution (unsigned int N,
    t_matrix *Lij,
    t_vector *rhs,
    t_vector *buffer) {

  t_matrix id;
  unsigned int i;

  /* allocate the memory for the gsl matrices */
  gsl_matrix *L = gsl_matrix_alloc (3*N, 3*N);
  gsl_vector *Y = gsl_vector_alloc (3*N);
  gsl_vector *X = gsl_vector_alloc (3*N);
  gsl_vector *LY = gsl_vector_alloc (3*N);
  gsl_vector *diff = gsl_vector_alloc (3*N);

  matrix_set_identity (id);

  /* build the gsl matrix and vectors */
  for (i=0; i<N; i++) {
    /* L */
    set_gsl_matrix_from_matrix (i,i, id, L);
    if (i>0)
      set_gsl_matrix_from_matrix (i,i-1, Lij [i-1], L);

    /* X and Y */
    set_gsl_vector_from_vector (i, rhs [i], X);
    set_gsl_vector_from_vector (i, buffer [i], Y);
  }

  /* compute LY */
  gsl_blas_dgemv (CblasNoTrans, 1., L, Y, 0., LY);

  /* check the equation L Y = X */
  gsl_vector_memcpy (diff, LY);
  gsl_vector_sub (diff, X);

  /* print results */
  printf ("L Y = \n");
  print_gsl_vector (3*N, LY);
  printf ("X = \n");
  print_gsl_vector (3*N, X);
  printf ("diff = \n");
  print_gsl_vector (3*N, diff);

  /* free memory */
  gsl_matrix_free (L);
  gsl_vector_free (Y);
  gsl_vector_free (X);
  gsl_vector_free (LY);
  gsl_vector_free (diff);
}



/* check the equation D L^t l = Y, where
 * D = Dii, L = Lij, l = lambda, Y = buffer */
void check_back_substitution (unsigned int N,
    t_matrix *Lij,
    t_matrix *Dii,
    t_vector *lambda,
    t_vector *buffer) {

  t_matrix id;
  unsigned int i;

  /* allocate the memory for the gsl matrices */
  gsl_matrix *L = gsl_matrix_alloc (3*N, 3*N);
  gsl_matrix *D = gsl_matrix_alloc (3*N, 3*N);
  gsl_matrix *DLt = gsl_matrix_alloc (3*N, 3*N);
  gsl_vector *Y = gsl_vector_alloc (3*N);
  gsl_vector *l = gsl_vector_alloc (3*N);
  gsl_vector *DLtl = gsl_vector_alloc (3*N);
  gsl_vector *diff = gsl_vector_alloc (3*N);

  matrix_set_identity (id);

  /* build the gsl matrix and vectors */
  for (i=0; i<N; i++) {
    /* D */
    set_gsl_matrix_from_matrix (i,i, Dii [i], D);

    /* L */
    set_gsl_matrix_from_matrix (i,i, id, L);
    if (i>0)
      set_gsl_matrix_from_matrix (i,i-1, Lij [i-1], L);

    /* lambda and Y */
    set_gsl_vector_from_vector (i, lambda [i], l);
    set_gsl_vector_from_vector (i, buffer [i], Y);
  }

  /* calculate D L^t */
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1., D, L, 0., DLt);

  /* compute D L^t l */
  gsl_blas_dgemv (CblasNoTrans, 1., DLt, l, 0., DLtl);

  /* check the equation D L^t l = Y */
  gsl_vector_memcpy (diff, DLtl);
  gsl_vector_sub (diff, Y);

  /* print results */
  printf ("D L^t l = \n");
  print_gsl_vector (3*N, DLtl);
  printf ("Y = \n");
  print_gsl_vector (3*N, Y);
  printf ("diff = \n");
  print_gsl_vector (3*N, diff);

  /* free memory */
  gsl_matrix_free (L);
  gsl_matrix_free (D);
  gsl_matrix_free (DLt);
  gsl_vector_free (Y);
  gsl_vector_free (l);
  gsl_vector_free (DLtl);
  gsl_vector_free (diff);
}



/* check H l = X, where
 * H = Hii,Hij  l = lambda,  X = rhs */
void check_LDL_solution (unsigned int N,
  t_matrix *Hii,
  t_matrix *Hij,
  t_vector *lambda,
  t_vector *rhs) {
  unsigned int i;

  /* allocate the memory for the gsl matrices */
  gsl_matrix *H = gsl_matrix_alloc (3*N, 3*N);
  gsl_vector *l = gsl_vector_alloc (3*N);
  gsl_vector *X = gsl_vector_alloc (3*N);
  gsl_vector *Hl = gsl_vector_alloc (3*N);
  gsl_vector *diff = gsl_vector_alloc (3*N);

  /* build the gsl matrix and vectors */
  for (i=0; i<N; i++) {
    /* H */
    set_gsl_matrix_from_matrix (i,i, Hii [i], H);
    if (i>0) {
      t_matrix transp;
      set_gsl_matrix_from_matrix (i,i-1, Hij [i-1], H);
      matrix_transpose (transp, Hij [i-1]);
      set_gsl_matrix_from_matrix (i-1,i, transp, H);
    }

    /* X and l */
    set_gsl_vector_from_vector (i, lambda [i], l);
    set_gsl_vector_from_vector (i, rhs [i], X);
  }

  /* calculate H l */
  gsl_blas_dgemv (CblasNoTrans, 1., H, l, 0., Hl);

  /* check Hl = X */
  gsl_vector_memcpy (diff, Hl);
  gsl_vector_sub (diff, X);

  /* print results */
  printf ("Hl:\n");
  print_gsl_vector (3*N, Hl);
  printf ("X:\n");
  print_gsl_vector (3*N, X);
  printf ("diff:\n");
  print_gsl_vector (3*N, diff);

  /* free memory */
  gsl_matrix_free (H);
  gsl_vector_free (Hl);
  gsl_vector_free (l);
  gsl_vector_free (X);
  gsl_vector_free (diff);
}
