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

#ifndef __DNA_ODE_ANALYZE_H__
#define __DNA_ODE_ANALYZE_H__


#include "config_read.h"
#include "trajectory_tools.h"



void print_usage_analyze (const char *program_name);



/* prototypes of the functions we can invoke as analysis types */

int ta_info (int argc, char *argv [], config_t *config);

int ta_writhe (int argc, char *argv [], config_t *config);

int ta_angle_hist (int argc, char *argv [], config_t *config);

int ta_distance_hist (int argc, char *argv [], config_t *config);

int ta_distance_map (int argc, char *argv [], config_t *config);

int ta_frame_dump (int argc, char *argv [], config_t *config);

int ta_body_trace (int argc, char *argv [], config_t *config);

int ta_bead_force (int argc, char *argv [], config_t *config);

int ta_last_segment_position (int argc, char *argv [], config_t *config);

#endif
