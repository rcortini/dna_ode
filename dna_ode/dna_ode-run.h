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

#ifndef __DNA_ODE_RUN_H__
#define __DNA_ODE_RUN_H__

#include "dna_ode_core.h"
#include "simparameters.h"
#include "simcontext.h"
#include "output_functions.h"
#include "print_functions.h"
#include "collision_callback.h"
#include "fin_functions.h"

void print_usage_run (const char *program_name);
int sdna_main (simcontext *simcon);
int fdna_main (simcontext *simcon);
int cdna_main (simcontext *simcon);

#endif
