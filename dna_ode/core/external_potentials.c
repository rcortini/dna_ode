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

#include "external_potentials.h"

/* this function adds the external forces to a DNA molecule */
t_real external_forces_add (external_potential_f *external_U, dna_ode *dna) {
  if (external_U->f == NULL)
    return 0.;
  else {
    unsigned int i;
    t_real ui, utot = 0.;
    t_vector fi;
    const unsigned int N = dna->nsegments;

    for (i=0; i<N; i++) {
      const dBodyID *seg_i = &dna->seg [i]->b_id;
      external_potential (*seg_i, external_U, &ui, fi);
      utot += ui;
      body_add_force (*seg_i, fi);
    }

    return utot;
  }
}





/* a generic function that returns the value of the potential energy of a body interacting with an external field */
void external_potential (const dBodyID b, external_potential_f *U, t_real *u12, t_real *f12) {
  U->f (b, U->p, u12, f12);
}



/* an external potential which couples the body to the z=z0 plane */
void harmonic_z_potential (dBodyID b, void *p, t_real *u, t_real *f) {
  struct harmonic_z_parameters *par = (struct harmonic_z_parameters *) p;
  t_real z = dBodyGetPosition (b)[2];

  /* potential */
  *u = par->k * (z-par->z0)*(z-par->z0);

  /* force */
  f[0] = 0.;
  f[1] = 0.;
  f[2] = -par->k * (z-par->z0);
}



/* Morse potential acting between the plane z=z0 and any body */
void morse_z_potential (dBodyID b, void *p, t_real *u, t_real *f) {
  struct morse_z_parameters *par = (struct morse_z_parameters *) p;
  t_real z = dBodyGetPosition (b)[2];
  t_real e = exp (-(z-par->z0)/par->lambda);
  t_real k = (1.-e);

  /* potential */
  *u = par->eps * (k*k - 1.);

  /* force */
  f[0] = 0.;
  f[1] = 0.;
  f[2] = -2 * par->eps / par->lambda * e * k;
}
