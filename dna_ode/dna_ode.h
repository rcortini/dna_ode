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

#ifndef __DNA_ODE_H__
#define __DNA_ODE_H__

#include <libconfig.h>
#include "dna_ode-run.h"

#define DNA_ODE_NAME "dna_ode"

int dna_ode_run (int argc, char **argv, config_t *simconfig);
int dna_ode_analyze (int argc, char **argv, config_t *simconfig);
int dna_ode_getp (int argc, char **argv, config_t *simconfig);

#endif
