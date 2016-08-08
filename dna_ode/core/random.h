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

#ifndef __RANDOM_H__
#define __RANDOM_H__




void my_init_RNG (unsigned long seed);

double ran_uniform (double min, double max);

double ran_gaussian_01 (void);

double ran_gaussian (double mean, double sigma);

void my_free_RNG ();


#endif
