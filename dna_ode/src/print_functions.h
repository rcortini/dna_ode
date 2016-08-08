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

#ifndef __PRINT_FUNCTIONS_H__
#define __PRINT_FUNCTIONS_H__

#include "dna_ode_core.h"
#include "bead.h"
#include "dna_ode.h"

char * which_body (const unsigned int nbodies, sim_body *body_list, unsigned int body_id);

void print_molecule_velocities (const unsigned int nbodies, sim_body *body_list);

void print_molecule_forces (const unsigned int nbodies, sim_body *body_list);

void print_molecule_torques (const unsigned int nbodies, sim_body *body_list);

void print_molecule_kinetic_energies (const unsigned int nbodies, sim_body *body_list);

void print_body_info (sim_body b);

void print_bead_info (simcontext *simcon, bead_ode *bead, sim_body *bead_body);

void print_dna_info (simcontext *simcon, dna_segment *seg, sim_body *seg_body);

#endif
