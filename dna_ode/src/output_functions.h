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

#ifndef __OUTPUT_FUNCTIONS_H__
#define __OUTPUT_FUNCTIONS_H__

#include "dna_ode_core.h"
#include "dna_ode.h"

void string_geometry_data (char * line, dGeomID g);

void string_body_data (char * line, dBodyID g);

void write_frame (const unsigned int nbodies, t_real sim_clock, FILE *trajectory_file, FILE *restart_file, sim_body *body_list);

void finish_output (simcontext * simcon);

#endif
