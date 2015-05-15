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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define MT_RNG_WARMUP_ROUNDS 50




/*
 * RANDOM NUMBER GENERATION FUNCTIONS 
 */




/* global RNG */
gsl_rng * RNG = (gsl_rng *)0;



/* initializes the random number generator, and should be invoked before the init_simulation_context function
 * and after the init_simulation_parameters function. A seed should be provided by input from the user */
void my_init_RNG (unsigned long seed) {
  int i;
  RNG = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set (RNG, seed);

  /*
    do a few "warm up" rounds 
    because MT is affected by small seeds
    (whose internal reprensetation as bits sequences
    contain a lot of successive 0's or 1's)
   */
  for (i=0; i<MT_RNG_WARMUP_ROUNDS; i++) gsl_rng_set (RNG, gsl_rng_get (RNG));
}

/* returns a random number uniformly distributed in the interval [min, max) */
double ran_uniform (double min, double max) {
  return gsl_ran_flat (RNG, min, max);
}

/* returns a random number distributed in the gaussian N (0,1) distribution */
double ran_gaussian_01 (void) {
  /* apparently the Ziggurat method is the fastest one available (source: GSL manual) */
  return gsl_ran_gaussian_ziggurat (RNG, 1.);
}

/* returns a random number distributed in the gaussian N (mean,sigma) distribution */
double ran_gaussian (double mean, double sigma) {
  return mean + gsl_ran_gaussian_ziggurat (RNG, sigma);
}

/* free random number generator */
void my_free_RNG () {
  gsl_rng_free (RNG);
}
