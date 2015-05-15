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

#include "potentials.h"

/* this function calculates the force between all the DNA segments, and returns the total potential energy of the system */
t_real force_loop (potential_f *U, dna_ode *dna) {
  unsigned int i, j;
  t_real uij, utot = 0.;
  t_vector fij;
  const unsigned int N = dna->nsegments;

  /* this is the famous force loop */
  for (i=2; i<N; i++) {
    /* we do not calculate interaction between connected DNA segments */
    for (j=0; j<i-1; j++) {
      const dBodyID *seg_i = &dna->seg [i]->b_id;
      const dBodyID *seg_j = &dna->seg [j]->b_id;
      potential (*seg_i, *seg_j, U, &uij, fij);
      utot += uij;
      body_add_force (*seg_i, fij);
      dNegateVector3 (fij);
      body_add_force (*seg_j, fij);
    }
  }

  return utot;
}



/* a generic function that returns the value of the potential energy of two interacting bodies */
void potential (const dBodyID b1, const dBodyID b2, potential_f *U, t_real *u12, t_real *f12) {
  U->f (b1, b2, U->p, u12, f12);
}





/* this function returns the value of the Lennard-Jones pair interaction potential. The 
 * parameters of the potential are eps and sigma. In the same function, the value of the radial
 * force is also computed. The formulas are:
 * u (r) = 4 epsilon *[(sigma/r)^12 - (sigma/r)^6]
 * F (r) = 24 epsilon / r^2 *[2 (sigma/r)^12 - (sigma/r)^6]  */
void LJ_pair_potential_and_force (const dBodyID b1, const dBodyID b2, void *p, t_real *u12, t_real *f12) {
  t_real R2, force, r, temp, x2, x6, x12;
  const t_real *r1, *r2;
  t_vector r12;
  struct LJ_parameters *par = (struct LJ_parameters *) p;

  /* body positions */
  r1 = dBodyGetPosition (b1);
  r2 = dBodyGetPosition (b2);
  vector_diff (r12, r1, r2);

  /* distance between r1 and r2, squared */
  R2 = dCalcVectorLengthSquare3 (r12);
  
  /* fast calculation of each term */
  r = sqrt (R2);
  temp = par->sigma/r;
  x2 = temp*temp;
  x6 = x2*x2*x2;
  x12 = x6*x6;

  /* this is u (r) */
  *u12 = 4*par->eps*(x12-x6);

  /* this is the force, calculated as f = -(1/r)(d/dr) u (r).
   * The force then acting on the particles will be _F_ = f _r_ */
  force = 24*par->eps/R2 * (2*x12-x6);

  /* set the values of the force vector. This is the force acting on b1 due to the presence of b2 */
  vector_scale (f12, r12, force);
}



/* takes user-supplied data on the LJ parameters, and fills the LJ_parameters structure */
void LJ_parameters_calculate (t_real l, t_real LJ_sigma, t_real LJ_eps, struct LJ_parameters *p) {
  p->sigma = LJ_sigma;
  p->eps = l * LJ_eps;
}



/* this function returns the value of the Todd pair interaction potential, from
 * Todd, Parsegian, ... Rau "Attractive forces between cation condensed DNA 
 * double helices" Biophys. J. 94, 4775--4782 (2008).
 * The potential was re-adapted by Argudo and Purohit to describe the interaction of 
 * only a _pair_ of DNA molecules, in
 * Argudo and Purohit, "Competition between supercoils and toroids in single molecule 
 * DNA condensation", Biophys. J. 103 (1), 1--11 (2012).
 * The potential has the form:
 * U (R)/l = sqrt(3) lambda/(e^(2R/lambda)) [CR/4 (2R+lambda) - CA (R + lambda) e^(R/lambda)]
 * Here, r is the inter-DNA distance, L indicates the length of the DNA, lambda is 
 * the decay length, measured to be 4.8 A. Imposing that the osmotic pressure is zero at the
 * equilibrium distance r_eq (which is measured experimentally), we have that
 * CR = CA exp (r_eq/lambda), so that in this potential there is only one free parameter, which
 * is CA. Defining A = l sqrt(3) lambda CA, and B = exp (r_eq/lambda), the potential can
 * be written as
 * U (R) = A [ B (2R + lambda) e^(-2R/lambda) - (R+lambda) e^(-R/lambda) ].
 * These are the parameters A and B defined in the Todd_parameters structure. */
void Todd_pair_potential_and_force (const dBodyID b1, const dBodyID b2, void *p, t_real *u12, t_real *f12) {
  t_real R, e, e2, force;
  const t_real *r1, *r2;
  t_vector r12;
  struct Todd_parameters *par = (struct Todd_parameters *) p;

  /* body positions */
  r1 = dBodyGetPosition (b1);
  r2 = dBodyGetPosition (b2);
  vector_diff (r12, r1, r2);

  /* distance between r1 and r2 */
  R = dCalcVectorLength3 (r12);
  
  /* calculate the exponential */ 
  e = exp (-R/par->lambda);
  e2 = e*e;

  /* this is u (r) */
  *u12 = par->A * (par->B / 4. * (2*R+par->lambda)*e2 - (R+par->lambda)*e);

  /* this is the force, calculated as f = -(1/r)(d/dr) u (r).
   * The force then acting on the particles will be _F_ = f _r_ */
  force = par->A/par->lambda * (par->B*e2 - e);

  /* set the values of the force vector. This is the force acting on b1 due to the presence of b2 */
  vector_scale (f12, r12, force);
}


/* takes user-supplied data on the Todd parameters, and fills the Todd_parameters structure */
void Todd_parameters_calculate (t_real l, t_real CA, t_real r_eq, t_real lambda, struct Todd_parameters *p) {
  p->A = l*M_SQRT3*lambda*CA;
  p->B = exp (r_eq/lambda);
  p->lambda = lambda;
}



/* this function returns the value of the truncated Lennard-Jones pair interaction potential.
 * This is obtained by setting a value of the distance at which we suppose that the interaction
 * energy is zero. Such value is r* = r_trunc * LJ_sigma, where r_trunc is typically chosen to
 * be 2.5. To avoid discontinuities in the energy, the potential is
 * V_trunc (r) = V_LJ (r) - V_LJ (r*). The form of the force is exactly the same. */
void LJ_trunc_pair_potential_and_force (const dBodyID b1, const dBodyID b2, void *p, t_real *u12, t_real *f12) {
 t_real R2;
 const t_real *r1, *r2;
 t_vector r12;
 struct LJ_trunc_parameters *par = (struct LJ_trunc_parameters *) p;

 /* body positions */
 r1 = dBodyGetPosition (b1);
 r2 = dBodyGetPosition (b2);
 vector_diff (r12, r1, r2);

 /* distance between r1 and r2, squared */
 R2 = dCalcVectorLengthSquare3 (r12);

 if (R2>par->rc2) {
   *u12 = 0.;
   vector_set_zero (f12);
   return;
 }
 else {
   t_real force, r, temp, x2, x6, x12;

   /* fast calculation of each term */
   r = sqrt (R2);
   temp = par->sigma/r;
   x2 = temp*temp;
   x6 = x2*x2*x2;
   x12 = x6*x6;

   /* this is u (r) */
   *u12 = 4*par->eps*(x12-x6) - par->V0;

   /* this is the force, calculated as f = -(1/r)(d/dr) u (r).
    * The force then acting on the particles will be _F_ = f _r_ */
   force = 24*par->eps/R2 * (2*x12-x6);

   /* set the values of the force vector. This is the force acting on b1 due to the presence of b2 */
   vector_scale (f12, r12, force);
 }
}



/* takes user-supplied data on the LJ parameters, and fills the LJ_parameters structure */
void LJ_trunc_parameters_calculate (t_real l, t_real LJ_sigma, t_real LJ_eps, t_real r_trunc, struct LJ_trunc_parameters *p) {
  t_real rc;
 
  /* set parameters */
  p->sigma = LJ_sigma;
  p->eps = l * LJ_eps;

  /* calculate rc = r_trunc * LJ_sigma */
  rc = r_trunc*p->sigma;
  p->rc2 = rc*rc;

  /* calculate V0 = V_LJ (rc) */
  p->V0 = 4*p->eps*(pow (p->sigma/rc,12)-pow (p->sigma/rc, 6));
}
